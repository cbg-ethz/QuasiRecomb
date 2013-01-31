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
package ch.ethz.bsse.quasirecomb.simulation;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.utils.BitMagic;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.FutureTask;
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Simulator {

    public final static String newline = System.getProperty("line.separator");

    public static void fromHaplotypesGlobalPaired(String[] haplotypes, int N, int L, double epsilon, double[] hapProb, String savePath) {
        Map<Byte, Boolean> map = new HashMap<>();
        for (String h : haplotypes) {
            for (int i = 0; i < h.length(); i++) {
                map.put((byte) h.charAt(i), Boolean.TRUE);
            }
        }
        int n = map.keySet().size();
        int insertSize = 200;
        int readLength = 250;
        int fragmentSize = insertSize + 2 * readLength;
        L = haplotypes[0].length();
        Map<Integer, Double> freqMap = new ConcurrentHashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            freqMap.put(i, hapProb[i]);
        }
        Frequency<Integer> frequency = new Frequency<>(freqMap);

        Random rand = new Random();

        List<FutureTask<Pair<Read, Read>>> taskList = new ArrayList<>();
        int coverage[] = new int[L];
        for (int i = 0; i < N; i++) {
            int start = 0;
            int length = 0;
            length = readLength;
            if (i > ((readLength * N) / L)) {
                start = rand.nextInt(L - fragmentSize);
            }
            if (i > N - ((250d * N) / L)) {
                start = L - fragmentSize;
            }
            for (int j = 0; j < readLength; j++) {
                coverage[start + j]++;
            }
            int start2 = start;
            if (Globals.getINSTANCE().isOVERLAP()) {
                start2 += readLength + rand.nextInt(2 * insertSize) - insertSize;
                for (int j = 0; j < readLength; j++) {
                    coverage[start2 + j]++;
                }
            } else {
                start2 += readLength + insertSize;
            }
            int hap = frequency.roll();
            FutureTask<Pair<Read, Read>> futureTask_1 = new FutureTask<>(new CallableSimulator(length, epsilon, n, haplotypes, hap, start, start2));
            taskList.add(futureTask_1);
            Globals.getINSTANCE().getExecutor().execute(futureTask_1);
            Globals.getINSTANCE().print("Preparation\t" + Math.round((100d * i) / N) + "%");
        }
        System.out.println("");


        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < taskList.size(); j++) {
            FutureTask<Pair<Read, Read>> futureTask = taskList.get(j);
            try {
                Pair<Read, Read> pair = futureTask.get();
                Read r = pair.getValue0();
                sb.append("@Read").append(j).append("\n");
                String s = Utils.reverse(r.getSequence(), r.getLength()).replaceAll("-", "");
                sb.append(s).append("\n");
                sb.append("+\n");
                for (int i = 0; i < s.length(); i++) {
                    sb.append("I");
                }
                sb.append("\n");
                r = pair.getValue1();
                sb.append("@Read").append(j).append("\n");
                s = Utils.reverse(r.getSequence(), r.getLength()).replaceAll("-", "");
                sb.append(s).append("\n");
                sb.append("+\n");
                for (int i = 0; i < s.length(); i++) {
                    sb.append("I");
                }
                sb.append("\n");
                Globals.getINSTANCE().print("Simulation\t" + Math.round((100d * j) / taskList.size()) + "%");
            } catch (InterruptedException | ExecutionException ex) {
                System.err.println("Problem");
                System.err.println(ex);
            }
        }
        System.out.println("");
        if (savePath.endsWith(".fastq")) {
            Utils.saveFile(savePath, sb.toString());
        } else {
            Utils.saveFile(savePath + ".fastq", sb.toString());
        }
        sb.setLength(0);
        for (int i = 0; i < L; i++) {
            sb.append(coverage[i]).append("\n");
        }
        Utils.saveFile(savePath + "_coverage.txt", sb.toString());

    }

    public static void fromHaplotypesGlobalAmplicon(String[] haplotypes, int N, int L, double epsilon, double[] hapProb, String savePath, String amplicons, String ampliconDistribution, int length) {
        int ampliconPos[];
        String[] split = amplicons.split(",");
        ampliconPos = new int[split.length];
        int x = 0;
        for (String s : split) {
            ampliconPos[x++] = Integer.parseInt(s);
        }

        double ampliconDist[];
        split = ampliconDistribution.split(",");
        ampliconDist = new double[split.length];
        x = 0;
        double sum = 0d;
        for (String s : split) {
            ampliconDist[x++] = Double.parseDouble(s);
            sum += ampliconDist[x - 1];
        }
        if (sum != 1d && Math.abs(sum - 1d) > 1e-6) {
            throw new RuntimeException("Amplicon distribution do not add up to 1, instead to " + sum);
        }

        Map<Integer, Double> ampliconMap = new ConcurrentHashMap<>();
        for (int i = 0; i < ampliconDist.length; i++) {
            ampliconMap.put(i, ampliconDist[i]);
        }
        Frequency<Integer> ampliconFreq = new Frequency<>(ampliconMap);


        Map<Byte, Boolean> map = new HashMap<>();
        for (String h : haplotypes) {
            for (int i = 0; i < h.length(); i++) {
                map.put((byte) h.charAt(i), Boolean.TRUE);
            }
        }
        int n = map.keySet().size();
        L = haplotypes[0].length();
        Map<Integer, Double> freqMap = new ConcurrentHashMap<>();
        for (int i = 0; i < hapProb.length; i++) {
            freqMap.put(i, hapProb[i]);
        }
        Frequency<Integer> frequency = new Frequency<>(freqMap);

        List<FutureTask<Read>> taskList = new ArrayList<>();
        int coverage[] = new int[L];
        for (int i = 0; i < N; i++) {
            int l = length;
            int start = ampliconPos[ampliconFreq.roll()];
            if (length + start > haplotypes[0].length()) {
                l = haplotypes[0].length() - start;
            }

            for (int j = 0; j < l; j++) {
                coverage[start + j]++;
            }
            int hap = frequency.roll();
            FutureTask<Read> futureTask_1 = new FutureTask<>(new CallableSimulatorSingle(l, epsilon, n, haplotypes, hap, start));
            taskList.add(futureTask_1);
            Globals.getINSTANCE().getExecutor().execute(futureTask_1);
            Globals.getINSTANCE().print("Preparation\t" + Math.round((100d * i) / N) + "%");
        }
        System.out.println("");


        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < taskList.size(); j++) {
            FutureTask<Read> futureTask = taskList.get(j);
            try {
                Read r = futureTask.get();
                sb.append(">SAMPLED").append(j).append("_").append(r.getBegin()).append("-").append(r.getEnd()).append("|").append(j).append("/1").append("\n");
                sb.append(Utils.reverse(r.getSequence(), r.getLength())).append("\n");
                Globals.getINSTANCE().print("Simulation\t" + Math.round((100d * j) / taskList.size()) + "%");
            } catch (InterruptedException | ExecutionException ex) {
                System.err.println("Problem");
                System.err.println(ex);
            }
        }
        System.out.println("");
        if (savePath.endsWith(".fasta")) {
            Utils.saveFile(savePath, sb.toString());
        } else {
            Utils.saveFile(savePath + ".fasta", sb.toString());
        }
        sb.setLength(0);
        for (int i = 0; i < L; i++) {
            sb.append(coverage[i]).append("\n");
        }
        Utils.saveFile(savePath + "_coverage.txt", sb.toString());

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
                    char x = Utils.reverseChar(v);
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
                    char x = Utils.reverseChar(v);
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
                    char x = Utils.reverseChar(v);
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
        Utils.saveFile(savePath + ".fasta", sb2.toString());
        return map;
    }
}
