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
package ch.ethz.bsse.quasirecomb.modelsampling;

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.*;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public final class ModelSampling extends Utils {

    private String savePath = "";
    private int L;
    private int n;
    private final int K;
    private double[][][] rho;
    private double[] pi;
    private double[][][] mu;
    private double[][] tauOmega;
    private Frequency<Integer>[][] rhoArray;
    private Frequency<Byte>[][] muArray;
    private Map<byte[], Integer> reads;
    private int amount = Globals.getINSTANCE().getSAMPLING_NUMBER();
    private int[] recombPerObservation;
    private Map<String, Double> hexMap = new HashMap<>();
    private StringBuilder sb = new StringBuilder();
    private StringBuilder startStopSB = new StringBuilder();
    private int[] coverage;

    public ModelSampling(OptimalResult or, String savePath) {
        this.K = or.getK();
        this.L = or.getL();
        this.n = or.getn();
        this.rho = or.getRho();
        this.pi = or.getPi();
        this.mu = or.getMu();
        this.rhoArray = new Frequency[L - 1][K];
        this.muArray = new Frequency[L][K];
        this.recombPerObservation = new int[amount];
        this.tauOmega = or.getTauOmega();
        this.coverage = new int[L];
        this.savePath = savePath;
        this.start();
    }

    public ModelSampling(String string, String path) {
        OptimalResult or = null;

        try {
            FileInputStream fis = new FileInputStream(string);
            try (ObjectInputStream in = new ObjectInputStream(fis)) {
                or = (OptimalResult) in.readObject();
            }
        } catch (IOException | ClassNotFoundException ex) {
            System.err.println(ex);
        }
        this.savePath = path;
        this.K = or.getK();
        this.L = or.getL();
        this.n = or.getn();
        this.rho = or.getRho();
        this.pi = or.getPi();
        this.mu = or.getMu();
        this.tauOmega = or.getTauOmega();
        this.rhoArray = new Frequency[L - 1][K];
        this.muArray = new Frequency[L][K];
        this.recombPerObservation = new int[amount];
        this.coverage = new int[L];
        this.start();
    }

    public void start() {
        if (!new File(savePath).exists()) {
            if (!new File(savePath).mkdirs()) {
                throw new RuntimeException("Cannot create directory: " + savePath);
            }
        }
        startStopSB.append("start\tstop\n");
        this.reads = new HashMap<>();
        byte[][] readArray = new byte[amount][L];
        for (int i = 0; i < amount; i++) {
            readArray[i] = single(i);
        }

        for (int i = 0; i < amount; i++) {
            byte[] read = readArray[i];
            boolean hit = false;
            for (byte[] readm : reads.keySet()) {
                if (Arrays.equals(read, readm)) {
                    reads.put(readm, reads.get(readm) + 1);
                    hit = true;
                    break;
                }
            }
            if (!hit) {
                reads.put(read, 1);
            }
        }

        int i = 0;
        for (Object o : sortMapByValue(reads).keySet()) {
            byte[] read = (byte[]) o;
            sb.append(">read").append(i++).append("_").append(((double) reads.get(read)) / amount).append("\n");
            for (int r : read) {
                sb.append(reverse(r));
            }
            sb.append("\n");
        }
    }
    
    public void saveQuasispeciesOnly(String path) {
        Utils.saveFile(path, sb.toString());
    }

    public void save() {
        Utils.saveFile(savePath + "quasispecies.fasta", sb.toString());
        Utils.saveFile(savePath + "support" + File.separator + "startStop.txt", startStopSB.toString());
        StringBuilder coverageSB = new StringBuilder();
        coverageSB.append("x\ty\n");
        for (int i = 0; i < L; i++) {
            coverageSB.append(i).append("\t").append(coverage[i]).append("\n");
        }
        Utils.saveFile(savePath+"support"+File.separator + "simCov.txt", coverageSB.toString());
    }

    public byte[] single(int currentI) {
        byte[] read = new byte[L];
        Map<Integer, Double> piMap = new HashMap<>();
        for (int k = 0; k < K; k++) {
            piMap.put(k, pi[k]);
        }
        Frequency<Integer> piF = new Frequency<>(piMap);
        int k = piF.roll();
        int oldk = k;

        Map<Integer, Double> startMap = new HashMap<>();
        for (int j = 0; j < L + 1; j++) {
            startMap.put(j, tauOmega[0][j]);
        }
        Frequency<Integer> startF = new Frequency<>(startMap);
        int start = startF.roll();
        startStopSB.append(start).append("\t");

        Map<Integer, Double> stopMap = new HashMap<>();
        for (int j = start + 1; j < L + 1; j++) {
            stopMap.put(j, tauOmega[1][j]);
        }
        Frequency<Integer> stopF = new Frequency<>(stopMap);
        int stop = stopF.roll();
        startStopSB.append(stop).append("\n");

        for (int j = start; j < stop; j++) {
            this.coverage[j]++;
        }

        for (int j = 0; j < L; j++) {
            if (j > 0) {
                Map<Integer, Double> rhoMap = new HashMap<>();
                for (int l = 0; l < K; l++) {
                    rhoMap.put(l, rho[j - 1][k][l]);
                }
                Frequency<Integer> rhoF = new Frequency<>(rhoMap);
                rhoArray[j - 1][k] = rhoF;
                k = rhoArray[j - 1][k].roll();
                if (oldk != k) {
                    this.recombPerObservation[currentI] += 1;
                }
                oldk = k;
            }
            if (muArray[j][k] == null) {
                Map<Byte, Double> muMap = new HashMap<>();
                for (byte v = 0; v < n; v++) {
                    muMap.put(v, mu[j][k][v]);
                }
                Frequency<Byte> muF = new Frequency<>(muMap);
                muArray[j][k] = muF;
            }
            read[j] = muArray[j][k].roll();
        }
        return read;
    }

    public Map<String, Integer> getReadsReversed() {
        Map<String, Integer> m = new HashMap<>();
        for (byte[] b : this.reads.keySet()) {
            m.put(reverse(b), this.reads.get(b));
        }
        return m;
    }

    public Map<String, Double> getHexMap() {
        return hexMap;
    }

    public void printQuasispecies() {
        System.out.println(sb.toString());
    }

    public Map<byte[], Integer> getReads() {
        return reads;
    }

    public static Map sortMapByValue(Map map) {
        List listForSort;
        Map sortedList = new LinkedHashMap();
        listForSort = new LinkedList(map.entrySet());
        Collections.sort(listForSort, new Comparator() {
            @Override
            public int compare(Object value1, Object value2) {
                return ((Comparable) ((Map.Entry) (value2)).getValue()).compareTo(((Map.Entry) (value1)).getValue());
            }
        });
        Iterator itret = listForSort.iterator();
        while (itret.hasNext()) {
            Map.Entry entry = (Map.Entry) itret.next();
            sortedList.put(entry.getKey(), entry.getValue());
        }
        return sortedList;
    }

    public int getK() {
        return K;
    }
}

class ValueComparator implements Comparator, Serializable {

    Map base;
    private static final long serialVersionUID = 1L;

    public ValueComparator(Map base) {
        this.base = base;
    }

    @Override
    public int compare(Object a, Object b) {
        if (((Integer) base.get(a)).intValue() < ((Integer) base.get(b)).intValue()) {
            return 1;
        } else if (((Integer) base.get(a)).intValue() == ((Integer) base.get(b)).intValue()) {
            return 0;
        } else {
            return -1;
        }
    }
}