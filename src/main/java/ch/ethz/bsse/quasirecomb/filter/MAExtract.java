/**
 * Copyright (c) 2011-2012 Armin Töpfer
 *
 * This file is part of QuasiRecomb.
 *
 * QuasiRecomb is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * QuasiRecomb is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * QuasiRecomb. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.quasirecomb.filter;

import java.io.*;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class MAExtract {

    private static boolean verbose;
    private static int update;
    private static int window = 200;
    private static int increment = 1;
    private static double gapc;
    private static boolean filter_only = false;

    public static void filter(StringBuilder sb, String head, BufferedWriter out) throws IOException {
        String s = sb.toString();
        int begin = 0;
        int end = 0;
        for (int j = 0; j < s.length(); j++) {
            if (s.charAt(j) != '-') {
               begin = j; 
               break;
            }
        }
        for (int j = s.length()-1; j >=begin; j--) {
            if (s.charAt(j) != '-') {
               end = j; 
               break;
            }
        }
        
        if (acceptableRegion(s.substring(begin, end))) {
            out.write(head + "\n" + s + "\n");
        }
        
//        for (int j = 0; j < s.length(); j += increment) {
//            int x = j + window;
//            if (acceptableRegion(s.substring(j, x < s.length() ? x : s.length()))) {
//                while (x <= s.length()) {
//                    x++;
//                    if (!acceptableRegion(s.substring(j, x < s.length() ? x : s.length()))) {
//                        x--;
//                        break;
//                    }
//                }
//                out.write(head + "###" + j + "-" + x + "\n");
//                if (filter_only) {
//                    out.write(s + "\n");
//                } else {
//                    StringBuilder gb = new StringBuilder();
//                    for (int i = 0; i < j; i++) {
//                        gb.append("-");
//                    }
//                    gb.append(s.substring(j, x));
//                    for (int i = 0; i < x; i++) {
//                        gb.append("-");
//                    }
//                    out.write(gb.toString() + "\n");
//                }
//                break;
//            }
//        }
    }

    public static void calc(String input, String output, double p) {
        gapc = p;
        filter_only = true;
        parse(input, output, false);
    }

    /**
     * @param args the command line arguments
     */
//    public static void main(String[] args) throws JSAPException {
//
//        SimpleJSAP jsap = new SimpleJSAP("MAExtract", "",
//                new Parameter[]{
//                    new FlaggedOption("Input file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'i', "input", "Input multiple alignment fasta file."),
//                    new FlaggedOption("Output file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'o', "output", "Output multiple alignment fasta file."),
//                    new FlaggedOption("Begin index", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'b', "begin", "Begin index. Sequence starts with 1."),
//                    new FlaggedOption("End index", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'e', "end", "End index."),
//                    new FlaggedOption("Gap percentage", JSAP.DOUBLE_PARSER, "0.1", JSAP.NOT_REQUIRED, 'g', "gap", "Amount of gaps allowed."),
//                    new FlaggedOption("Window size", JSAP.INTEGER_PARSER, "300", JSAP.NOT_REQUIRED, 'w', "window", "Tests each read by given window size if filter succeeds."),
//                    new FlaggedOption("Increment", JSAP.INTEGER_PARSER, "10", JSAP.NOT_REQUIRED, 'c', "increment", "Slides the window with given increment."),
//                    new FlaggedOption("Update", JSAP.INTEGER_PARSER, "1", JSAP.REQUIRED, 'u', "Update interval for verbose mode."),
//                    new Switch("Verbose").setLongFlag("verbose").setHelp("Print debug information")
//                });
//        JSAPResult config = jsap.parse(args);
//        if (config.success()) {
//            if (config.contains("Begin index") && config.contains("End index")) {
//                begin = config.getInt("Begin index");
//                end = config.getInt("End index");
//            } else {
//                filter_only = true;
//            }
//            window = config.getInt("Window size");
//            increment = config.getInt("Increment");
//            String input = config.getString("Input file");
//            String output = config.getString("Output file");
//            verbose = config.getBoolean("Verbose");
//            update = config.getInt("Update");
//            gapc = config.getDouble("Gap percentage");
//            if (!new File(input).exists()) {
//                System.out.println("Input file does not exist.");
//                System.exit(0);
//            }
//            parse(input, output);
//        }
//    }
    public static void parse(String input, String output, boolean complete) {
        try {
            // Create file 
            FileWriter outstream = new FileWriter(output);
            try (BufferedWriter out = new BufferedWriter(outstream)) {
                FileInputStream fstream = new FileInputStream(input);
                StringBuilder sb;
                String head = null;
                int i = 0;
                try (DataInputStream in = new DataInputStream(fstream)) {
                    BufferedReader br = new BufferedReader(new InputStreamReader(in));
                    String strLine;
                    sb = new StringBuilder();

                    while ((strLine = br.readLine()) != null) {
                        if (strLine.startsWith(">")) {
                            if (verbose) {
                                i++;
                                if (i % update == 0) {
                                    System.out.println(i);
                                }
                            }
                            if (sb.length() > 0) {
                                if (complete) {
                                    if (acceptableRegion(sb.toString())) {
                                        out.write(head);
                                        out.write("\n");
                                        out.write(sb.toString());
                                        out.write("\n");
                                    }
                                } else {
                                    filter(sb, head, out);
                                }
                                sb.setLength(0);
                            }
                            head = strLine;
                        } else {
                            sb.append(strLine);
                        }
                    }

                }
                if (complete) {
                    if (acceptableRegion(sb.toString())) {
                        out.write(head);
                        out.write("\n");
                        out.write(sb.toString());
                        out.write("\n");
                    }
                } else {
                    filter(sb, head, out);
                }
            }
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
    }

    private static boolean acceptableRegion(String s) {
//        if (s.startsWith("-")) {
//            return false;
//        }
//        if (s.endsWith("-")) {
//            return false;
//        }
//        if (s.contains("N")) {
//            return false;
//        }
        int gaps = 0;
        for (char c : s.toCharArray()) {
            if (c == '-') {
                gaps++;
            }
        }
        if (gaps > (s.length() * gapc)) {
            return false;
        }
//        for (char c : s.toCharArray()) {
//            if (c != '-') {
//                return true;
//            }
//        }
//        return false;
        return true;
    }
}
