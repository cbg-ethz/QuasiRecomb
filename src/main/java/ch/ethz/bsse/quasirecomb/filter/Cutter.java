/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.filter;

import ch.ethz.bsse.quasirecomb.model.Globals;
import java.io.*;

/**
 *
 * @author toepfera
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
                            if (Globals.DEBUG) {
                                i++;
                                if (i % 100 == 0) {
                                    System.out.println(i);
                                }
                            }
                            if (sb.length() > 0) {
                                out.write(head);
                                out.write("\n");
                                out.write(sb.substring(begin - 1, end));
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
                out.write(sb.substring(begin - 1, end));
                out.write("\n");
            }
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
    }
}
