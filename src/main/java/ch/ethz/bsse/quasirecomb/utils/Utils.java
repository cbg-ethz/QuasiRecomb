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

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import java.io.*;
import java.util.*;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Utils extends FastaParser {

    public final static String SAVEPATH = "";

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
//            return null;
        }

    }

    public static Read[] parseBAMSAM(String location) {
        Map<String, Read> readMap = new HashMap<>();
        File bam = new File(location);
        SAMFileReader sfr = new SAMFileReader(bam);
        for (final SAMRecord samRecord : sfr) {
            List<AlignmentBlock> alignmentBlocks = samRecord.getAlignmentBlocks();
            if (alignmentBlocks.isEmpty()) {
                continue;
            }
            int refStart = alignmentBlocks.get(0).getReferenceStart() + alignmentBlocks.get(0).getReadStart() - 1;
            int readStart = alignmentBlocks.get(0).getReadStart() - 1;
//            boolean ignore = false;
            List<Byte> buildRead = new ArrayList<>();
            boolean d = false;
            for (CigarElement c : samRecord.getCigar().getCigarElements()) {
                switch (c.getOperator()) {
                    case M:
                        for (int i = 0; i < c.getLength(); i++) {
                            if (readStart >= samRecord.getReadBases().length) {
                                System.out.println("U1");
                                System.exit(9);
                            }
                            byte b = samRecord.getReadBases()[readStart++];
                            buildRead.add(b);
                        }
                        break;
                    case I:
//                        for (int i = 0; i < c.getLength(); i++) {
//                            buildRead.add((byte) "-".charAt(0));
                        readStart++;
//                        }
                        break;
                    case D:
                        d = true;
                        for (int i = 0; i < c.getLength(); i++) {
                            buildRead.add((byte) "-".charAt(0));
//                            readStart++;
                        }
                        break;
                    case S:
                        break;
                    case H:
                        System.out.println("H");
                        System.exit(9);
                        break;
                    case P:
                        System.out.println("P");
                        System.exit(9);
                        break;
                    case EQ:
                        System.out.println("EQ");
                        System.exit(9);
                        break;
                    case X:
                        System.out.println("X");
                        System.exit(9);
                        break;
                    case N:
                        System.out.println("N");
                        System.exit(9);
                        break;
                    default:
                        break;
                }
            }
            byte[] readBases = convertRead(buildRead.toArray(new Byte[buildRead.size()]));
            String name = samRecord.getReadName();
            if (readMap.containsKey(name)) {
                readMap.get(name).setPairedEnd(BitMagic.pack(readBases), refStart, refStart + readBases.length);
                Globals.getINSTANCE().incMERGED();
            } else {
                readMap.put(name, new Read(BitMagic.pack(readBases), refStart, refStart + readBases.length, 1));
            }
        }
        Map<Integer, Read> hashed = new HashMap<>();
        for (Read r1 : readMap.values()) {
            int hash = r1.hashCode();
            if (hashed.containsKey(hash)) {
                hashed.get(hash).incCount();
            } else {
                hashed.put(hash, r1);
            }
        }
        return hashed.values().toArray(new Read[hashed.size()]);
    }

    private static byte[] convertRead(Byte[] readSplit) {
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

//    public static Read[] parseInput(String path) {
//        byte[][] haplotypesArray = splitReadsIntoByteArrays(parseFarFile(path));
//        List<Read> hashing = new ArrayList<>();
//        for (byte[] b : haplotypesArray) {
//            boolean missing = true;
//            for (Read r : hashing) {
//                if (Arrays.equals(r.getSequence(), b)) {
//                    r.incCount();
//                    missing = false;
//                    break;
//                }
//            }
//            if (missing) {
//                hashing.add(new Read(b, 0, b.length, 1));
//            }
//        }
//        return hashing.toArray(new Read[hashing.size()]);
//    }
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
                    hashing.add(new Read(seq, begin, end, 1));
                }
            }
        } else {
            Map<Integer, Read> hashMap = new HashMap<>();
            String[] parseFarFile = parseFarFile(path);
            for (String s : parseFarFile) {
                byte[] packed = BitMagic.splitReadIntoBytes(s);
                Read r = new Read(packed, 0, s.length(), 1);
                if (hashMap.containsKey(r.hashCode())) {
                    hashMap.get(r.hashCode()).incCount();
                } else {
                    hashMap.put(r.hashCode(), new Read(packed, 0, s.length(), 1));
                }
            }
            hashing.addAll(hashMap.values());
//            byte[][] haplotypesArray = splitReadsIntoByteArrays(parseFarFile(path));
//            for (byte[] b : haplotypesArray) {
//                boolean missing = true;
//                for (Read r : hashing) {
//                    if (Arrays.equals(r.getSequence(), b)) {
//                        r.incCount();
//                        missing = false;
//                        break;
//                    }
//                }
//                if (missing) {
//                    hashing.add(new Read(b, 0, b.length, 1));
//                }
//            }
        }
        return hashing.toArray(new Read[hashing.size()]);
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

//    public static Map<String, Integer> reverse(Map<byte[], Integer> src) {
//        Map<String, Integer> dest = new HashMap<>();
//        for (Map.Entry<byte[], Integer> bb : src.entrySet()) {
//            StringBuilder sb = new StringBuilder(bb.getKey().length);
//            
//            for (byte b : bb.getKey()) {
//                sb.append(reverse(b));
//            }
//            dest.put(sb.toString(), bb.getValue());
//        }
//        return dest;
//    }
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
