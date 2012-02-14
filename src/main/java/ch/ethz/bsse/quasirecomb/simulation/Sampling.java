package ch.ethz.bsse.quasirecomb.simulation;

import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class Sampling {
    
    public static String newline = System.getProperty("line.separator");
    public static Random random = new Random();

    public static Map<String, Integer> fromHaplotypes(String path, int N, int L, double epsilon, double[] hapProb, int n, String savePath) {
        return fromHaplotypes(Utils.parseFarFile(path), N, L, epsilon, hapProb, n, savePath);
    }

    public static Map<String, Integer> fromHaplotypesGlobal(String path, int N, int L, double epsilon, double[] hapProb, int n, String savePath) {
        return fromHaplotypesGlobal(Utils.parseFarFile(path), N, L, epsilon, hapProb, n, savePath);
    }

    public static Map<String, Integer> fromHaplotypesGlobal(String[] haplotypes, int N, int L, double epsilon, double[] hapProb, int n, String savePath) {
        Map<Integer, Double> freqMap = new HashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            freqMap.put(i, hapProb[i]);
        }
        Frequency<Integer> frequency = new Frequency<>(freqMap);
        Map<Integer, Double> hapFreq = new HashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            hapFreq.put(i, 0d);
        }
        Map<String, Integer> map = new HashMap<>();
        String read;
        for (int i = 0; i < N; i++) {
            int hap = frequency.roll();
            hapFreq.put(hap, hapFreq.get(hap) + 1);
            char[] readArray = new char[L];
            int start = random.nextInt((int)(L*0.9));
            int length = 0;
            do {
                length = (int) (Math.random() * 200);
                if (start + length > L) {
                    if (L - start < 30) {
                        continue;
                    }
                    length = L - start;
                    break;
                }
            } while (length < 100);
//            for (int j = 0; j < start; j++) {
//                readArray[j] = '-';
//            }
            for (int j = start; j < start + length; j++) {
                //error
                Map<Character, Double> baseMap = new HashMap<>();
                for (int v = 0; v < n; v++) {
                    char x = reverse(v);
                    if (haplotypes[hap].charAt(j) == x) {
                        baseMap.put(x, 1.0 - (n - 1.0) * epsilon);
                    } else {
                        baseMap.put(x, epsilon);
                    }
                }
                Frequency<Character> errorF = new Frequency<>(baseMap);
                readArray[j] = errorF.roll();
            }
//            if (start+length < L) {
//                for (int j = start+length; j < L; j++) {
//                    readArray[j] = '-';
//                }
//            }
            StringBuilder s = new StringBuilder(length);
            for (int j =start; j < start+length; j++) {
                s.append(readArray[j]);
            }
//            read = String.valueOf(readArray);
            read = s.toString();
            if (!map.containsKey(read)) {
                map.put(read, 0);
            }
            map.put(read, map.get(read) + 1);
        }
        StringBuilder sb = new StringBuilder();
        for (String readX : map.keySet()) {
            sb.append(map.get(readX)).append("\t").append(readX).append("\n");
        }
        Utils.saveFile(savePath + "sampledReadDistribution.txt", sb.toString());
        int z = 0;
        sb.setLength(0);
        for (String readX : map.keySet()) {
            for (int i = 0; i < map.get(readX); i++) {
                sb.append(">SAMPLED-").append(z++).append(newline).append(readX).append("\n");
                
            }
        }
        Utils.saveFile(savePath + "reads.fasta", sb.toString());
        return map;
    }

    public static Map<String, Integer> fromHaplotypes(String[] haplotypes, int N, int L, double epsilon, double[] hapProb, int n, String savePath) {
        Map<Integer, Double> freqMap = new HashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            freqMap.put(i, hapProb[i]);
        }
        Frequency<Integer> frequency = new Frequency<>(freqMap);
        Map<Integer, Double> hapFreq = new HashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            hapFreq.put(i, 0d);
        }
        Map<String, Integer> map = new HashMap<>();
        String read;
        for (int i = 0; i < N; i++) {
            int hap = frequency.roll();
            hapFreq.put(hap, hapFreq.get(hap) + 1);
            char[] readArray = new char[L];
            for (int j = 0; j < L; j++) {
                //error
                Map<Character, Double> baseMap = new HashMap<>();
                for (int v = 0; v < n; v++) {
                    char x = reverse(v);
                    if (haplotypes[hap].charAt(j) == x) {
                        baseMap.put(x, 1.0 - (n - 1.0) * epsilon);
                    } else {
                        baseMap.put(x, epsilon);
                    }
                }
                Frequency<Character> errorF = new Frequency<>(baseMap);
                readArray[j] = errorF.roll();
            }
            read = String.valueOf(readArray);
            if (!map.containsKey(read)) {
                map.put(read, 0);
            }
            map.put(read, map.get(read) + 1);
        }
        StringBuilder sb = new StringBuilder();
        for (String readX : map.keySet()) {
            sb.append(map.get(readX)).append("\t").append(readX).append("\n");
        }
        Utils.saveFile(savePath + "sampledReadDistribution.txt", sb.toString());
        return map;
    }

    private static char reverse(int v) {
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
}
