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
package ch.ethz.bsse.quasirecomb.entropy;

import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ShannonEntropy {

    private static boolean verbose;

    /**
     * @param args the command line arguments
     */
    public static void entropy(String input) {

        verbose = Globals.DEBUG;
        String output;
        File file = new File(input);
        if (file.exists()) {
            Map<String, byte[]> far = parseFarFile(input);
            int N = far.size();
            int L = far.values().iterator().next().length;
            double[] positionEntropy = new double[L];
            int[] coverage = new int[L];
            
            for (int j = 0; j < L; j++) {
                int l = 0;
                Map<Byte, Integer> map = new HashMap<>();
                for (String p : far.keySet()) {
                    byte c = far.get(p)[j];
                    if (c != 'N') {
                        l++;
                        if (!map.containsKey(c)) {
                            map.put(c, 0);
                        }
                    }
                    map.put(c, map.get(c) + 1);
                }

                // calculate the entropy
                Double result = 0.0;
                for (byte base : map.keySet()) {
                    Double frequency = (double) map.get(base) / l;
                    result -= frequency * Math.log(frequency);
                }
                positionEntropy[j] = result;
                coverage[j] = l;
            }
            StringBuilder sb = new StringBuilder();
            for (double p : positionEntropy) {
                sb.append(p).append("\n");
            }
            StringBuilder cb = new StringBuilder();
            for (int i : coverage) {
                cb.append(i).append("\n");
            }
            String pwd = new File(".").getAbsolutePath();
            if (input.contains(File.separator)) {
                output = input.split(File.separator)[input.split(File.separator).length - 1];
            } else {
                output = input;
            }
            if (output.contains(".")) {
                System.out.println(Arrays.toString(output.split("\\.")));
                output = output.split("\\.")[0];
            }
            output = pwd.subSequence(0, pwd.length() - 1) + output;
            System.out.println(output);
            Utils.saveFile(output + "Data", sb.toString());
            Utils.saveFile(output + "Coverage", cb.toString());
            sb.setLength(0);
            sb.append("pdf(\"" + output + "Histogram.pdf\")\n");
            sb.append("data <- read.delim(\"" + output + "Data\", header=F)\n");
            sb.append("plot(seq(1,length(data$V1)),data$V1,type='l')\n");
            sb.append("dev.off()");
            Utils.saveFile(output + "Histogram.R", sb.toString());
        }
    }

    public static Map<String, byte[]> parseFarFile(String location) {
        Map<String, byte[]> map = new HashMap<>();
        try {
            FileInputStream fstream = new FileInputStream(location);
            StringBuilder sb;
            String head = null;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (sb.length() > 0) {
                            map.put(head, Utils.splitReadIntoByteArray(sb.toString()));
                        }
//                        String xy[] = strLine.split("###")[1].split("-");
//                        Pair p = Pair.with(Integer.parseInt(xy[0]), Integer.parseInt(xy[1]));
                        sb.setLength(0);
                        head = strLine;
                    } else {
                        sb.append(strLine);
                    }
                }
            }
            map.put(head, Utils.splitReadIntoByteArray(sb.toString()));
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return map;
    }
}
