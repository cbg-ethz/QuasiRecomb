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
package ch.ethz.bsse.quasirecomb.quasiviz;

import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.*;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class QuasiViz {

    /**
     * @param args the command line arguments
     */
//    public static void main(String[] args) {
////        String path = "C:\\Users\\XLR\\Dropbox\\jHMM\\QuasiSimulation\\src\\sim\\K2\\simuHap.fasta";
//        String name = "K3";
//        String output = "C:\\Users\\XLR\\Dropbox\\QRData\\viz\\SRR002680\\"+name;
//        String input = "C:\\Users\\XLR\\Dropbox\\jHMM\\QuasiSimulation\\src\\sim\\"+name+"\\simuHap.fasta";
//        paint(input,output);
//    }
    public static void paint(String path, String output) throws NumberFormatException {
        paint(path, output, null);
    }

    public static void paint(String path, String output, String dot) throws NumberFormatException {
//        String path = "C:\\Users\\XLR\\Dropbox\\QRData\\viz\\SRR069887\\"+input+".fasta";
        Pair<String, String>[] reads = parseFarFile(path);
        Pair<String, char[]>[] chars = new Pair[reads.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < reads.length; i++) {
            chars[i] = Pair.with(reads[i].getValue0(), reads[i].getValue1().toCharArray());
        }
        sb.append("digraph G {\ngraph [rankdir=LR];\nnode [shape=box];\nnode [shape=box, color=white];\nedge [color=white]\n");
        List<Integer> l = new LinkedList<>();
        for (int i = 1; i < chars.length; i++) {
            for (int j = 0; j < chars[i].getValue1().length; j++) {
                if (!l.contains(j)) {
                    if (chars[i].getValue1()[j] != chars[0].getValue1()[j]) {
                        l.add(j);
                    }
                }
            }
        }
        snps = l.toArray(new Integer[l.size()]);
        for (int i = 0; i < chars.length; i++) {
            boolean isHaplotype = false;
            try {
                Integer.parseInt(chars[i].getValue0());
                isHaplotype = true;
            } catch (NumberFormatException ex) {
            }

            if (!isHaplotype) {
                sb.append("subgraph ").append(chars[i].getValue0()).append(" {\n");
                sb.append("\"").append(chars[i].getValue0()).append("\"").append(" -> ");
                for (int j = 0; j < chars[i].getValue1().length; j++) {
                    if (region(j)) {
                        node(sb, i, chars, j).append(" -> ");
                    }
                }
                sb.setLength(sb.length() - 3);
                sb.append(";\n}\n");
                for (int j = 0; j < chars[i].getValue1().length; j++) {
                    if (region(j)) {
                        if (chars[i].getValue1()[j] != chars[0].getValue1()[j]) {
                            node(sb, i, chars, j).append(" [style=filled, color=\"#" + color(chars[i].getValue1()[j]) + "\", fontcolor=\"black\"];\n");
                        } else {
                            node(sb, i, chars, j).append(" [fontcolor=white];\n");
                        }
                    }
                }
            } else if (isHaplotype && Double.parseDouble(chars[i].getValue0()) > 0 && Integer.parseInt(chars[i].getValue0()) < 99999) {
                sb.append("subgraph ").append(chars[i].getValue0()).append(" {\n");
                sb.append("\"").append(chars[i].getValue0()).append("\"").append(" -> ");
                for (int j = 0; j < chars[i].getValue1().length; j++) {
                    if (region(j)) {
                        node(sb, i, chars, j).append(" -> ");
                    }
                }
                sb.setLength(sb.length() - 3);
                sb.append(";\n}\n");
                for (int j = 0; j < chars[i].getValue1().length; j++) {
                    if (region(j)) {
                        if (chars[i].getValue1()[j] != chars[0].getValue1()[j]) {
                            node(sb, i, chars, j).append(" [style=filled, color=orange, fontcolor=black];\n");
                        } else {
                            node(sb, i, chars, j).append(" [fontcolor=white];\n");
                        }
                    }
                }
            }
        }
        sb.append("}");
        Utils.saveFile(output + ".gv", sb.toString());
        if (dot != null) {
            try {
                Runtime.getRuntime().exec(dot + " -Tpdf -o " + output + ".pdf " + output + ".gv");
//            Runtime.getRuntime().exec("C:\\Program Files (x86)\\Graphviz 2.28\\bin\\dot.exe -Tpdf -o "+output+".pdf "+output+".gv");
            } catch (IOException ex) {
                Logger.getLogger(QuasiViz.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    private static Integer[] snps;// = new int[]{41, 49, 77, 79, 80, 97, 112, 113, 114, 116, 117,};

    private static String color(char j) {
        if (j == 'A') {
            return "61AE24";
        } else if (j == 'C') {
            return "D61F3B";
        } else if (j == 'G') {
            return "FCA700";
        } else if (j == 'T') {
            return "4878A8";
        } else {
            return "9a9a9a";
        }
    }

    private static boolean region(int j) {
        for (int snp : snps) {
            if (j == snp) {
                return true;
            }
        }
        return false;
//        return (j==15 || j ==16 || j == 28 || j == 78 || j == 83 || j == 84 || j == 153 || j == 160 || j == 161 || j == 171 || j == 175 || j == 187 || j == 188 || j == 255 || j == 279 || j == 288);
    }

    public static StringBuilder node(StringBuilder sb, int i, Pair<String, char[]>[] chars, int j) {
        return sb.append("\"").append(i).append(" ").append(chars[i].getValue1()[j]).append(" ").append(j).append("\" ");
    }

    public static Pair[] parseFarFile(String location) {
        List<Pair> readList = new LinkedList<>();
        try {
            FileInputStream fstream = new FileInputStream(location);
            String header = null;
            StringBuilder sb;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (sb.length() > 0) {
                            readList.add(Pair.with(header.replaceAll(">", "").replaceAll("-", ""), sb.toString()));
                            sb.setLength(0);
                        }
                        header = strLine;
                    } else {
                        sb.append(strLine);
                    }
                }
            }
            readList.add(Pair.with(header.replaceAll(">", "").replaceAll("-", ""), sb.toString()));
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return readList.toArray(new Pair[readList.size()]);
    }
}
