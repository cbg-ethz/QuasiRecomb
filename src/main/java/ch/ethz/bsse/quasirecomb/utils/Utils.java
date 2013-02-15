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
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.informationholder.ReadTMP;
import com.google.common.collect.Lists;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Utils extends FastaParser {

    public final static String SAVEPATH = "";

    public static byte[] convertRead(Byte[] readSplit) {
        byte[] rs = new byte[readSplit.length];
        int length = readSplit.length;
        for (int i = 0; i < length; i++) {
            switch ((short) readSplit[i]) {
                case 65:
                    rs[i] = 0;
                    break;
                case 67:
                    rs[i] = 1;
                    break;
                case 71:
                    rs[i] = 2;
                    break;
                case 84:
                    rs[i] = 3;
                    break;
                case 45:
                case 78:
                    rs[i] = 4;
                    break;
                default:
                    System.out.println("Unknown " + (char) ((byte) readSplit[i]) + " " + readSplit[i]);
                    break;
            }
        }
        return rs;
    }

    private static boolean isFastaFormat(String path) {
        try {
            FileInputStream fstream = new FileInputStream(path);
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        return true;
                    } else {
                        return false;
                    }
                }
            }
        } catch (Exception e) {
            System.err.println("Error identifying format of input file: " + e.getMessage());
        }
        return false;
    }

    private static boolean isFastaGlobalFormat(String path) {
        try {
            FileInputStream fstream = new FileInputStream(path);
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">") && strLine.contains("_") && strLine.contains("-")) {
                        String[] split = strLine.split("_")[1].split("-");
                        try {
                            int begin = Integer.parseInt(split[0]);
                            int end = Integer.parseInt(split[1]);
                        } catch (NumberFormatException e) {
                            return false;
                        }
                        return true;
                    } else {
                        return false;
                    }
                }
            }
        } catch (Exception e) {
            System.err.println("Error identifying format of input file: " + e.getMessage());
        }
        return false;
    }

    private static boolean isFastaGlobalMatePairFormat(String path) {
        try {
            FileInputStream fstream = new FileInputStream(path);
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">") && strLine.contains("/")) {
                        return true;
                    } else {
                        return false;
                    }
                }
            }
        } catch (Exception e) {
            System.err.println("Error identifying format of input file: " + e.getMessage());
        }
        return false;
    }

    public static void mkdir(String save) {
        if (!new File(save).exists()) {
            if (!new File(save).mkdirs()) {
                throw new RuntimeException("Cannot create directory: " + save);
            }
        }
    }

    public static void saveOptimum(String save, OptimalResult or) {
        try {
            FileOutputStream fos = new FileOutputStream(save);
            try (ObjectOutputStream out = new ObjectOutputStream(fos)) {
                out.writeObject(or);
            }
        } catch (IOException ex) {
            System.out.println("Optimum Java saving\n" + ex.getMessage());
        }
    }

    public static void appendFile(String path, String sb) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(path, true);
            try (BufferedWriter out = new BufferedWriter(fstream)) {
                out.write(sb);
            }
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error append file: ");
            System.err.println(path);
            System.err.println(sb);
        }
    }

    public static void saveFile(String path, String sb) {
        try {
            // Create file 
            FileWriter fstream = new FileWriter(path);
            try (BufferedWriter out = new BufferedWriter(fstream)) {
                out.write(sb);
            }
        } catch (Exception e) {//Catch exception if any
            System.err.println("Error save file: ");
            System.err.println(path);
            System.err.println(sb);
        }
    }

    public static byte[][] convertSamples(int N, Map<String, Integer> map) {
        String[] reads = new String[N];
        int index = 0;
        for (String read : map.keySet()) {
            for (int i = 0; i < map.get(read); i++) {
                reads[index++] = read;
            }
        }
        return splitReadsIntoByteArrays(reads);
    }

    public static String reverse(byte[] packed, int length) {
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < length; j++) {
            sb.append(reverse(BitMagic.getPosition(packed, j)));
        }
        return sb.toString();
    }

    /**
     * Splits String consisting of {A,C,G,T,-} into corresponding byte array
     * with alphabet {0,1,2,3,4}.
     *
     * @param read String containing the read
     * @return byte array corresponding the read
     */
    public static byte[][] splitReadsIntoByteArrays(String[] reads) {
        byte[][] Rs = new byte[reads.length][reads[0].length()];
        for (int x = 0; x < reads.length; x++) {
            Rs[x] = BitMagic.splitReadIntoBytes(reads[x]);
        }
        return Rs;
    }

    public static byte[] splitReadIntoByteArray(String read) {
        byte[] Rs = new byte[read.length()];
        char[] readSplit = read.toCharArray();
        int length = readSplit.length;
        for (int i = 0; i < length; i++) {
            switch ((short) readSplit[i]) {
                case 65:
                    Rs[i] = 0;
                    break;
                case 67:
                    Rs[i] = 1;
                    break;
                case 71:
                    Rs[i] = 2;
                    break;
                case 84:
                    Rs[i] = 3;
                    break;
                case 45:
                    Rs[i] = 4;
                    break;
                default:
                    break;
            }
        }
        return Rs;
    }

    public static byte[] splitQualityIntoByteArray(String quality) {
        byte[] Rs = new byte[quality.length()];
        char[] readSplit = quality.toCharArray();
        int length = readSplit.length;
        for (int i = 0; i < length; i++) {
            Rs[i] = (byte) readSplit[i];
        }
        return Rs;
    }

    public static void saveClusteredReads(Map<byte[], Integer> clusterReads) {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<byte[], Integer> read : clusterReads.entrySet()) {
            sb.append(read.getValue()).append("\t");
            for (byte b : read.getKey()) {
                sb.append(Utils.reverse((int) b));
            }
            sb.append("\n");
        }
        Utils.saveFile(SAVEPATH + "support" + File.separator + "readsClustered.txt", sb.toString());
    }

    public static Read[] parseInput(String path) {
        if (isFastaGlobalMatePairFormat(path)) {
            return FastaParser.parseFastaPairedEnd(path);
        } else if (isFastaFormat(path)) {
            return parseFastaInput(path);
        } else {
            return parseBAMSAM(path);
        }
    }

    public static Read[] parseBAMSAM(String location) {
        Map<String, Read> readMap = new HashMap<>();
        File bam = new File(location);
        SAMFileReader sfr = new SAMFileReader(bam);
        double size = 0;
        for (final SAMRecord samRecord : sfr) {
            size++;
        }
        sfr.close();
        sfr = new SAMFileReader(bam);
        List<Future<ReadTMP>> readFutures = Lists.newArrayListWithExpectedSize((int) size);
        int counter = 0;
        for (final SAMRecord samRecord : sfr) {
            readFutures.add(Globals.getINSTANCE().getExecutor().submit(new SFRComputing(samRecord)));
            Globals.getINSTANCE().print("Parsing\t\t" + (Math.round((counter++ / size) * 100)) + "%");
        }
        Globals.getINSTANCE().print("Parsing\t\t100%");
        sfr.close();
        StringBuilder sb = new StringBuilder();
        for (Future<ReadTMP> future : readFutures) {
            try {
                ReadTMP read = future.get();
                if (read != null) {
                    String name = read.name;
                    int refStart = read.refStart;
                    byte[] readBases = read.readBases;
                    double[] quality = read.quality;
                    boolean hasQuality = read.hasQuality;
                    if (readMap.containsKey(name)) {
                        if (hasQuality) {
                            readMap.get(name).setPairedEnd(BitMagic.pack(readBases), refStart, refStart + readBases.length, quality);
                        } else {
                            readMap.get(name).setPairedEnd(BitMagic.pack(readBases), refStart, refStart + readBases.length);
                        }
                        Read r2 = readMap.get(name);
                        if (r2.isPaired()) {
                            sb.append(r2.getCrickBegin() - r2.getWatsonEnd()).append("\n");
                            if (r2.getCrickBegin() - r2.getWatsonEnd() < 0) {
                                System.out.println(name + "\t" + r2.getWatsonBegin() + "\t" + r2.getWatsonEnd() + "\t" + r2.getCrickBegin() + "\t" + r2.getCrickEnd());
//                            readMap.put(name, null);
                            }
                        }
                        if (Globals.getINSTANCE().isUNPAIRED()) {
                            Read r = readMap.get(name);
                            if (r.isPaired()) {
                                r.unpair();
                                if (hasQuality) {
                                    readMap.put(name + "_R", new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length, quality));
                                } else {
                                    readMap.put(name + "_R", new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length));
                                }
                            }
                        }
                        Globals.getINSTANCE().incPAIRED();
                    } else {
                        if (hasQuality) {
                            readMap.put(name, new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length, quality));
                        } else {
                            readMap.put(name, new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length));
                        }
                    }
                }
            } catch (InterruptedException | ExecutionException ex) {
                System.err.println(ex);
            }
        }
        readFutures.clear();

        Map<Integer, Read> hashed = new HashMap<>();
        for (Read r1 : readMap.values()) {
            if (r1 != null) {
                int hash = r1.hashCode();
                if (hashed.containsKey(hash)) {
                    hashed.get(hash).incCount();
                } else {
                    hashed.put(hash, r1);
                }
            }
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "indel.txt", sb.toString());
//        return null;
        return hashed.values().toArray(new Read[hashed.size()]);
    }

    public static Read[] parseFastaInput(String path) {
        List<Read> hashing = new ArrayList<>();
        if (isFastaGlobalFormat(path)) {
            Map<String, byte[]> haps = parseGlobalFarFile(path);
            for (Map.Entry<String, byte[]> head : haps.entrySet()) {
                String[] split = head.getKey().split("_")[1].split("-");
                int begin = Integer.parseInt(split[0]);
                int end = Integer.parseInt(split[1]);
                byte[] seq = head.getValue();
                boolean missing = true;
                for (Read r : hashing) {
                    if (Arrays.equals(r.getSequence(), seq)
                            && r.getBegin() == begin
                            && r.getEnd() == end) {
                        r.incCount();
                        missing = false;
                        break;
                    }
                }
                if (missing) {
                    hashing.add(new Read(seq, begin, end));
                }
            }
        } else {
            Map<Integer, Read> hashMap = new HashMap<>();
            String[] parseFarFile = parseFarFile(path);
            for (String s : parseFarFile) {
                byte[] packed = BitMagic.splitReadIntoBytes(s);
                Read r = new Read(packed, 0, s.length());
                if (hashMap.containsKey(r.hashCode())) {
                    hashMap.get(r.hashCode()).incCount();
                } else {
                    hashMap.put(r.hashCode(), new Read(packed, 0, s.length()));
                }
            }
            hashing.addAll(hashMap.values());
        }
        return hashing.toArray(new Read[hashing.size()]);
    }

    public static String reverse(int[] intArray) {
        StringBuilder sb = new StringBuilder();
        for (int i : intArray) {
            sb.append(reverse(i));
        }
        return sb.toString();
    }

    public static String reverse(byte[] bArray) {
        StringBuilder sb = new StringBuilder();
        for (byte b : bArray) {
            sb.append(reverse(b));
        }
        return sb.toString();
    }

    public static String reverse(byte i) {
        switch (i) {
            case 0:
                return "A";
            case 1:
                return "C";
            case 2:
                return "G";
            case 3:
                return "T";
            case 4:
                return "-";
        }
        throw new IllegalAccessError();
    }

    public static String reverse(int i) {
        switch (i) {
            case 0:
                return "A";
            case 1:
                return "C";
            case 2:
                return "G";
            case 3:
                return "T";
            case 4:
                return "-";
        }
        throw new IllegalAccessError("" + i);
    }

    public static char reverseChar(int v) {
        switch ((short) v) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            case 4:
                return '-';
            default:
                throw new IllegalStateException("cannot reverse " + v);
        }
    }

    public static void save(Map<String, Integer> map, String path) {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, Integer> read : map.entrySet()) {
            for (int i = 0; i < read.getValue(); i++) {
                sb.append(">").append(i).append("\n").append(read.getKey()).append("\n");
            }
        }
        saveFile(path, sb.toString());
    }

    public static void error() {
        System.out.println("    .o oOOOOOOOo                                            OOOo");
        System.out.println("    Ob.OOOOOOOo  OOOo.      oOOo.                      .adOOOOOOO");
        System.out.println("    OboO\"\"\"\"\"\"\"\"\"\"\"\".OOo. .oOOOOOo.    OOOo.oOOOOOo..\"\"\"\"\"\"\"\"\"'OO");
        System.out.println("    OOP.oOOOOOOOOOOO \"POOOOOOOOOOOo.   `\"OOOOOOOOOP,OOOOOOOOOOOB'");
        System.out.println("    `O'OOOO'     `OOOOo\"OOOOOOOOOOO` .adOOOOOOOOO\"oOOO'    `OOOOo");
        System.out.println("    .OOOO'            `OOOOOOOOOOOOOOOOOOOOOOOOOO'            `OO");
        System.out.println("    OOOOO                 '\"OOOOOOOOOOOOOOOO\"`                oOO");
        System.out.println("   oOOOOOba.                .adOOOOOOOOOOba               .adOOOOo.");
        System.out.println("  oOOOOOOOOOOOOOba.    .adOOOOOOOOOO@^OOOOOOOba.     .adOOOOOOOOOOOO");
        System.out.println(" OOOOOOOOOOOOOOOOO.OOOOOOOOOOOOOO\"`  '\"OOOOOOOOOOOOO.OOOOOOOOOOOOOO");
        System.out.println(" \"OOOO\"       \"YOoOOOOMOIONODOO\"`  .   '\"OOROAOPOEOOOoOY\"     \"OOO\"");
        System.out.println("    Y           'OOOOOOOOOOOOOO: .oOOo. :OOOOOOOOOOO?'         :`");
        System.out.println("    :            .oO%OOOOOOOOOOo.OOOOOO.oOOOOOOOOOOOO?         .");
        System.out.println("    .            oOOP\"%OOOOOOOOoOOOOOOO?oOOOOO?OOOO\"OOo");
        System.out.println("                 '%o  OOOO\"%OOOO%\"%OOOOO\"OOOOOO\"OOO':");
        System.out.println("                      `$\"  `OOOO' `O\"Y ' `OOOO'  o             .");
        System.out.println("    .                  .     OP\"          : o     .");
        System.out.println("                              :");
        System.out.println("                              .");
    }
}
