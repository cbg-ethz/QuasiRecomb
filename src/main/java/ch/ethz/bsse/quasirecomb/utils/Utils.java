package ch.ethz.bsse.quasirecomb.utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class Utils extends FastaParser {

    public static String SAVEPATH = "";

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
            Rs[x] = splitReadIntoByteArray(reads[x]);
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

    public static Map<byte[], Integer> clusterReads(byte[][] reads) {
        Map<byte[], Integer> readsMap = new HashMap<>();
        for (byte[] readm : reads) {
            boolean hit = false;
            for (byte[] read : readsMap.keySet()) {
                if (Arrays.equals(read, readm)) {
                    readsMap.put(read, readsMap.get(read) + 1);
                    hit = true;
                    break;
                }
            }
            if (!hit) {
                readsMap.put(readm, 1);
            }
        }
        return readsMap;
    }

    public static void saveClusteredReads(Map<byte[], Integer> clusterReads) {
        StringBuilder sb = new StringBuilder();
        for (byte[] read : clusterReads.keySet()) {
            sb.append(clusterReads.get(read)).append("\t");
            for (byte b : read) {
                sb.append(Utils.reverse((int) b));
            }
            sb.append("\n");
        }
        Utils.saveFile(SAVEPATH + "readsClustered.txt", sb.toString());
    }

    public static byte[][] parseInput(String path) {
        byte[][] haplotypesArray = splitReadsIntoByteArrays(parseFarFile(path));
        return haplotypesArray;
    }

    public static Map<String, Integer> reverse(Map<byte[], Integer> src) {
        Map<String, Integer> dest = new HashMap<>();
        for (byte[] bb : src.keySet()) {
            StringBuilder sb = new StringBuilder(bb.length);
            for (byte b : bb) {
                sb.append(reverse(b));
            }
            dest.put(sb.toString(), src.get(bb));
        }
        return dest;
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

    public static void save(Map<String, Integer> map, String path) {
        StringBuilder sb = new StringBuilder();
        for (String read : map.keySet()) {
            for (int i = 0; i < map.get(read); i++) {
                sb.append(">").append(i).append("\n").append(read).append("\n");
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
