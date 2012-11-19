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
package ch.ethz.bsse.quasirecomb.simulation;

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.utils.BitMagic;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Sampling {

    public final static String newline = System.getProperty("line.separator");

    public static void fromHaplotypesGlobalPaired(String[] haplotypes, int N, int L, double epsilon, double[] hapProb, String savePath) {
        Map<Byte, Boolean> map = new HashMap<>();
        for (String h : haplotypes) {
            for (int i = 0; i < h.length(); i++) {
                map.put((byte)h.charAt(i), Boolean.TRUE);
            }
        }
        int n = map.keySet().size();
        int insertSize = 50;
        int readLength = 250;
        int fragmentSize = insertSize + 2 * readLength;
        L = haplotypes[0].length();
        Map<Integer, Double> freqMap = new ConcurrentHashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            freqMap.put(i, hapProb[i]);
        }
        Frequency<Integer> frequency = new Frequency<>(freqMap);

        Read[] reads1 = new Read[N];
        Read[] reads2 = new Read[N];
        Random rand = new Random();
        for (int i = 0; i < N; i++) {
            int hap = frequency.roll();

            int start = 0;
            int length = 0;
            length = readLength;
            if (i > ((readLength * N) / L)) {
                start = rand.nextInt(L-fragmentSize);
            }
            if (i > N - ((250d * N) / L)) {
                start = L - fragmentSize;
            }
            char[] readArray = new char[length];


            for (int j = 0; j < length; j++) {
                //error
                if (epsilon > 0d) {
                    Map<Character, Double> baseMap = new ConcurrentHashMap<>();
                    for (int v = 0; v < n; v++) {
                        char x = reverse(v);
                        if (haplotypes[hap].charAt(j + start) == x) {
                            baseMap.put(x, 1.0 - (n - 1.0) * epsilon);
                        } else {
                            baseMap.put(x, epsilon);
                        }
                    }
                    Frequency<Character> errorF = new Frequency<>(baseMap);
                    readArray[j] = errorF.roll();
                } else {
                    readArray[j] = haplotypes[hap].charAt(j + start);
                }
            }
            StringBuilder sb = new StringBuilder(length);
            for (int j = 0; j < length; j++) {
                sb.append(readArray[j]);
            }
            reads1[i] = new Read(BitMagic.splitReadIntoBytes(sb.toString()), start, start + length);


            start += readLength + insertSize;

            for (int j = 0; j < length; j++) {
                //error
                if (epsilon > 0d) {
                    Map<Character, Double> baseMap = new ConcurrentHashMap<>();
                    for (int v = 0; v < n; v++) {
                        char x = reverse(v);
                        if (haplotypes[hap].charAt(j + start) == x) {
                            baseMap.put(x, 1.0 - (n - 1.0) * epsilon);
                        } else {
                            baseMap.put(x, epsilon);
                        }
                    }
                    Frequency<Character> errorF = new Frequency<>(baseMap);
                    readArray[j] = errorF.roll();
                } else {
                    readArray[j] = haplotypes[hap].charAt(j + start);
                }
            }
            StringBuilder sb2 = new StringBuilder(length);
            for (int j = 0; j < length; j++) {
                sb2.append(readArray[j]);
            }
            reads2[i] = new Read(BitMagic.splitReadIntoBytes(sb2.toString()), start, start + length);
        }
        int z = 0;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < N; i++) {
            Read r = reads1[i];
            sb.append(">SAMPLED").append(z).append("_").append(r.getBegin()).append("-").append(r.getEnd()).append("|").append(z).append("/1").append("\n");
            sb.append(Utils.reverse(r.getSequence(), r.getLength())).append("\n");
            r = reads2[i];
            sb.append(">SAMPLED").append(z).append("_").append(r.getBegin()).append("-").append(r.getEnd()).append("|").append(z).append("/2").append("\n");
            sb.append(Utils.reverse(r.getSequence(), r.getLength())).append("\n");
            z++;
        }

        Utils.saveFile(savePath, sb.toString());
    }

    public static void fromHaplotypesGlobal(String[] haplotypes, int N, int L, double epsilon, double[] hapProb, int n, String savePath) {
        Map<Integer, Double> freqMap = new ConcurrentHashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            freqMap.put(i, hapProb[i]);
        }
        Frequency<Integer> frequency = new Frequency<>(freqMap);

        List<Read> reads = new ArrayList<>();
        String read;
        for (int i = 0; i < N; i++) {
            int hap = frequency.roll();

            int start = 0;
            int length = 0;
            for (;;) {
//                length = (int)(Math.random()*300);
                length = 400;
                if (Math.random() > .5) {
                    length -= (int) (Math.random() * 100);
                } else {
                    length += (int) (Math.random() * 100);
                }
                start = (int) (Math.random() * (L + 300 + 300));
                start -= 300;
                if (start >= L) {
                    continue;
                }
                if (start + length <= L) {
                    if (start < 0) {
                        start = 0;
                    }
                    break;
                }
                if (start + length > L) {
                    length = L - start;
                    break;
                }
            }
            System.out.println(start);
            char[] readArray = new char[length];


            for (int j = 0; j < length; j++) {
                //error
                Map<Character, Double> baseMap = new ConcurrentHashMap<>();
                for (int v = 0; v < n; v++) {
                    char x = reverse(v);
                    if (haplotypes[hap].charAt(j + start) == x) {
                        baseMap.put(x, 1.0 - (n - 1.0) * epsilon);
                    } else {
                        baseMap.put(x, epsilon);
                    }
                }
                Frequency<Character> errorF = new Frequency<>(baseMap);
                readArray[j] = errorF.roll();
            }
            StringBuilder sb = new StringBuilder(length);
            for (int j = 0; j < length; j++) {
                sb.append(readArray[j]);
            }
            read = sb.toString();
            reads.add(new Read(BitMagic.splitReadIntoBytes(read), start, start + length));
        }
        int z = 0;
        StringBuilder sb = new StringBuilder();
        for (Read r : reads) {
            sb.append(">SAMPLED").append(z++).append("_").append(r.getBegin()).append("-").append(r.getEnd()).append("\n");
            sb.append(Utils.reverse(r.getSequence(), r.getLength())).append("\n");
        }
        Utils.saveFile(savePath, sb.toString());
    }

    public static String[] fromHaplotypesCross(String[] haplotypes, int N, int L, double epsilon, double[] hapProb, int n, String savePath) {
        Map<Integer, Double> freqMap = new ConcurrentHashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            freqMap.put(i, hapProb[i]);
        }
        Frequency<Integer> frequency = new Frequency<>(freqMap);
        Map<Integer, Double> hapFreq = new ConcurrentHashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            hapFreq.put(i, 0d);
        }
        String[] reads = new String[N];
        String read;
        for (int i = 0; i < N; i++) {
            int hap = frequency.roll();
            hapFreq.put(hap, hapFreq.get(hap) + 1);
            char[] readArray = new char[L];
            for (int j = 0; j < L; j++) {
                //error
                Map<Character, Double> baseMap = new ConcurrentHashMap<>();
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
            reads[i] = read;
        }
        Globals.getINSTANCE().setHAPLOTYPE_ARRAY_EMPIRICAL(new String[N / 10]);
        Map<String, Integer> map = new ConcurrentHashMap<>();
        for (int i = 0; i < N; i++) {
            if (i < N / 10) {
                Globals.getINSTANCE().getHAPLOTYPE_ARRAY_EMPIRICAL()[i] = reads[i];
            } else {
                if (!map.containsKey(reads[i])) {
                    map.put(reads[i], 0);
                }
                map.put(reads[i], map.get(reads[i]) + 1);
            }
        }
        StringBuilder sb = new StringBuilder();
        StringBuilder sb2 = new StringBuilder();
        int z = 0;
        for (Map.Entry<String, Integer> readX : map.entrySet()) {
            sb.append(readX.getValue()).append("\t").append(readX.getKey()).append("\n");
            for (int i = 0; i < readX.getValue(); i++) {
                sb2.append(">SAMPLED-").append(z++).append(newline).append(readX).append("\n");
            }
        }
        Utils.saveFile(savePath + "sampledReadDistribution.txt", sb.toString());
        Utils.saveFile(savePath + "reads.fasta", sb2.toString());
        return reads;
    }

    public static Map<String, Integer> fromHaplotypes(String[] haplotypes, int N, int L, double epsilon, double[] hapProb, int n, String savePath) {
        Map<Integer, Double> freqMap = new ConcurrentHashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            freqMap.put(i, hapProb[i]);
        }
        Frequency<Integer> frequency = new Frequency<>(freqMap);
        Map<Integer, Double> hapFreq = new ConcurrentHashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            hapFreq.put(i, 0d);
        }
        Map<String, Integer> map = new ConcurrentHashMap<>();
        L = haplotypes[0].length();
        String read;
        for (int i = 0; i < N; i++) {
            int hap = frequency.roll();
//            hapFreq.put(hap, hapFreq.get(hap) + 1);
            char[] readArray = new char[L];
            for (int j = 0; j < L; j++) {
                //error
                Map<Character, Double> baseMap = new ConcurrentHashMap<>();
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
        StringBuilder sb2 = new StringBuilder();
        int z = 0;
        for (Map.Entry<String, Integer> readX : map.entrySet()) {
            sb.append(readX.getValue()).append("\t").append(readX.getKey()).append("\n");
            for (int i = 0; i < readX.getValue(); i++) {
                sb2.append(">SAMPLED-").append(z++).append(newline).append(readX.getKey()).append("\n");
            }
        }
        Utils.saveFile(savePath + "_dist", sb.toString());
        Utils.saveFile(savePath +".fasta", sb2.toString());
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
