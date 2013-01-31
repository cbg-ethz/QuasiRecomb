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
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
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
        Map<String, Double> hapMap = new ConcurrentHashMap<>();
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
                if (head != null) {
                    double freq;
                    try {
                        freq = Double.parseDouble(head);
                    } catch (NumberFormatException e) {
                        freq = Double.parseDouble(head.split("_")[1]);
                    }
                    hapMap.put(sb.toString(), freq);
                }
            }
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return hapMap;
    }

    public static Map<String, String> parseHaplotypeFile(String location) {
        Map<String, String> hapMap = new ConcurrentHashMap<>();
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

    public static Read[] parseFastaPairedEnd(String location) {
        Map<String, Read> pairedReads1 = new HashMap<>();
        Map<String, Read> pairedReads2 = new HashMap<>();
        Map<String, byte[]> haps = parseGlobalFarFile(location);
        for (Map.Entry<String, byte[]> head : haps.entrySet()) {
            String[] firstsplit = head.getKey().split("\\|");
            String[] split = firstsplit[0].split("_")[1].split("-");
            int begin = Integer.parseInt(split[0]);
            int end = Integer.parseInt(split[1]);
            byte[] seq = head.getValue();
            String[] secondsplit = firstsplit[1].split("/");
            String tag = secondsplit[0];
            int pairedNumber = Integer.parseInt(secondsplit[1]);
            if (pairedNumber == 1) {
                if (pairedReads1.containsKey(tag)) {
                    pairedReads1.get(tag).incCount();
                } else {
                    pairedReads1.put(tag, new Read(seq, begin, end));
                }
            } else {
                if (pairedReads2.containsKey(tag)) {
                    pairedReads2.get(tag).incCount();
                } else {
                    pairedReads2.put(tag, new Read(seq, begin, end));
                }
            }
        }

        if (Globals.getINSTANCE().isPAIRED()) {
            for (Map.Entry<String, Read> r : pairedReads1.entrySet()) {
                String tag = r.getKey();
                if (pairedReads2.containsKey(tag)) {
                    Read r2 = pairedReads2.get(tag);
                    r.getValue().setPairedEnd(r2.getSequence().clone(), r2.getBegin(), r2.getEnd());
                    r.getValue().rearrange();
                    r.getValue().merge();
                } else {
//                System.out.println("");
//                pairedReads1.remove(tag);
                }
            }
        } else {
            for (Map.Entry<String, Read> r : pairedReads2.entrySet()) {
                pairedReads1.put(r.getKey() + "-", r.getValue());
            }
        }
        Map<Integer, Read> hashed = new HashMap<>();
        for (Read r1 : pairedReads1.values()) {
            int hash = r1.hashCode();
            if (hashed.containsKey(hash)) {
                hashed.get(hash).incCount();
            } else {
                hashed.put(hash, r1);
            }
        }
        return hashed.values().toArray(new Read[hashed.size()]);
    }

    public static Map<String, byte[]> parseGlobalFarFile(String location) {
        Map<String, byte[]> hapMap = new ConcurrentHashMap<>();
        try {
            String head = null;

            StringBuilder sb = new StringBuilder();
            BufferedReader reader =
                    Files.newBufferedReader(Paths.get(location), StandardCharsets.UTF_8);
            String strLine = null;
            while ((strLine = reader.readLine()) != null) {
                if (strLine.startsWith(">")) {
                    if (sb.length() > 0) {
                        hapMap.put(head, BitMagic.splitReadIntoBytes(sb.toString()));
                        sb = new StringBuilder();
                    }
                    head = strLine;
                } else {
                    sb.append(strLine);
                }
                if (hapMap.size() % 50000 == 0) {
                    System.gc();
                }
            }
            System.gc();
            System.gc();
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return hapMap;
    }
}
