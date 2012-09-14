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
package ch.ethz.bsse.quasirecomb.utils;

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava3.sequencing.io.fastq.Fastq;
import org.biojava3.sequencing.io.fastq.FastqReader;
import org.biojava3.sequencing.io.fastq.IlluminaFastqReader;

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

    public static Read[] parseFastq(String location) {
        FastqReader fastqReader = new IlluminaFastqReader();
        Map<String, Read> pairedReads = new HashMap<>();
        List<Read> hashing = new LinkedList<>();
        try {
            Iterable<Fastq> reads = fastqReader.read(new File(location));
            for (Fastq f : reads) {
                //@generator-0_0.25_899_1068_0:0:0_0:0:0_0/2
                byte[] seq = Utils.splitReadIntoByteArray(f.getSequence());
                String description = f.getDescription();
                String tag = description.split(":")[0];
                final String[] firstSplit = description.split("_");
                //SAMPLED-0_100-300\1
                int pairedNumber = Integer.parseInt(description.split("/")[1]);
                switch (pairedNumber) {
                    case 1:
                        int begin = Integer.parseInt(firstSplit[1]);
                        int end = begin + seq.length;
                        pairedReads.put(tag, new Read(seq, begin, end));
                        break;
                    case 2:
                        Read mate = pairedReads.get(tag);
                        int end2 = Integer.parseInt(firstSplit[2]);
                        int begin2 = end2-seq.length;
                        mate.setPairedEnd(seq, begin2, end2);
                        boolean missing = true;
                        for (Read r : hashing) {
                            if (r.equals(mate)) {
                                r.incCount();
                                missing = false;
                                break;
                            }
                        }
                        if (missing) {
                            mate.setCount(1);
                            hashing.add(mate);
                        }
                        break;
                    default:
                        throw new IllegalStateException("Do not know paired end number " + pairedNumber + " of read " + firstSplit[0]);
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(FastaParser.class.getName()).log(Level.SEVERE, null, ex);
        }
        return hashing.toArray(new Read[hashing.size()]);
    }

    public static Map<String, String> parseGlobalFarFile(String location) {
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
                            hapMap.put(head, sb.toString());
                            sb.setLength(0);
                        }
                        head = strLine;
                    } else {
                        sb.append(strLine);
                    }
                }
                hapMap.put(head, sb.toString());

            }
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return hapMap;
    }
}
