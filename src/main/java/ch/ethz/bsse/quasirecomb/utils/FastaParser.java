/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.utils;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 *
 * @author xlr
 */
public class FastaParser {

    /**
     *
     * @param location
     * @return
     */
    public static String[] parseFarFile(String location) {
        List<String> readList = new LinkedList<>();
        try {
            FileInputStream fstream = new FileInputStream(location);
            StringBuilder sb;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (sb.length() > 0) {
                            readList.add(sb.toString());
                            sb.setLength(0);
                        }
                    } else {
                        sb.append(strLine);
                    }
                }
                readList.add(sb.toString());
            }
//            readList.add(sb.toString());
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return readList.toArray(new String[readList.size()]);
    }

    /**
     *
     * @param location
     * @return
     */
    public static Map<String, Double> parseQuasispeciesFile(String location) {
        Map<String, Double> hapMap = new HashMap<>();
        try {
            FileInputStream fstream = new FileInputStream(location);
            StringBuilder sb;
            String head = null;
            int i = 0;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (sb.length() > 0) {
                            double freq;
                            try {
                                freq = Double.parseDouble(head);
                            } catch (NumberFormatException e) {
                                freq = Double.parseDouble(head.split("_")[1]);
                            }
                            hapMap.put(sb.toString(), freq);
                            sb.setLength(0);
                        }
                        head = strLine;
                    } else {
                        sb.append(strLine);
                    }
                }
                double freq;
                try {
                    freq = Double.parseDouble(head);
                } catch (NumberFormatException e) {
                    freq = Double.parseDouble(head.split("_")[1]);
                }
                hapMap.put(sb.toString(), freq);

            }
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return hapMap;
    }

    public static Map<String, String> parseHaplotypeFile(String location) {
        Map<String, String> hapMap = new HashMap<>();
        try {
            FileInputStream fstream = new FileInputStream(location);
            StringBuilder sb;
            String head = null;
            int i = 0;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (sb.length() > 0) {
                            hapMap.put(sb.toString(), head);
                            sb.setLength(0);
                        }
                        head = strLine;
                    } else {
                        sb.append(strLine);
                    }
                }
                hapMap.put(sb.toString(), head);

            }
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return hapMap;
    }
}
