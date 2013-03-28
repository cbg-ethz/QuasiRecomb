/**
 * Copyright (c) 2011-2013 Armin Töpfer
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
package ch.ethz.bsse.quasirecomb.utils;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import java.io.*;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Cutter {

    public static void cut(String input, String output, int begin, int end) {
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
                            if (Globals.getINSTANCE().isDEBUG()) {
                                i++;
                                if (i % 100 == 0) {
                                    System.out.println(i);
                                }
                            }
                            if (sb.length() > 0) {
                                out.write(head);
                                out.write("\n");
                                out.write(sb.toString().substring(begin, end));
                                out.write("\n");
                                sb.setLength(0);
                            }
                            head = strLine;
                        } else {
                            sb.append(strLine);
                        }
                    }

                }
                out.write(head);
                out.write("\n");
                out.write(sb.toString().substring(begin, end));
                out.write("\n");
            }
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
    }
}
