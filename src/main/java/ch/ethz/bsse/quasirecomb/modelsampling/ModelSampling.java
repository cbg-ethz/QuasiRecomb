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
package ch.ethz.bsse.quasirecomb.modelsampling;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.TauOmega;
import ch.ethz.bsse.quasirecomb.informationholder.Threading;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.StatusUpdate;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import static ch.ethz.bsse.quasirecomb.utils.Utils.reverse;
import com.google.common.collect.Lists;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;

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
    private TauOmega tauOmega;
    private Frequency<Integer>[][] rhoArray;
    private Frequency<Byte>[][] muArray;
    private Map<byte[], Integer> reads;
    private int amount = 10000;
    private int[] recombPerObservation;
    private Map<String, Double> hexMap = new HashMap<>();
    private StringBuilder sb = new StringBuilder();
    private Map<String, Double> map = new LinkedHashMap<>();
    private final OptimalResult or;
    private int[] coverage;

    public ModelSampling(OptimalResult or, String savePath) {
//        this.amount = or.getN();
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
        this.or = or;
        this.start();
    }

    public ModelSampling(String string, String path) {
        try {
            FileInputStream fis = new FileInputStream(string);
            try (ObjectInputStream in = new ObjectInputStream(fis)) {
                or = (OptimalResult) in.readObject();
            }
        } catch (IOException | ClassNotFoundException ex) {
            System.err.println(ex);
            Utils.error();
            throw new IllegalStateException("Optimal result could not be parsed");
        }
        if (or == null) {
            throw new IllegalStateException("Optimal result could not be parsed");
        }
//        this.amount = or.getN();
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
        StringBuilder readBuilder = new StringBuilder();
        StringBuilder startStopSB = new StringBuilder();
        startStopSB.append("start\tstop\tstart2\tstop2\n");
        this.reads = new HashMap<>();
        byte[][] readArray = new byte[amount][L];

        boolean paired = tauOmega.isPaired();
        List<Future<List<SampledRead>>> readFutures = Lists.newArrayListWithExpectedSize(amount);

        Map<Integer, Double> startMap = new HashMap<>();
        for (int j = 0; j < L; j++) {
            if (this.tauOmega.getTauWatsonMap().get(j) == null || this.tauOmega.getTauWatsonMap().get(j).isEmpty()) {
                continue;
            }
            double sum = 0d;
            for (Map.Entry<Integer, Double> e : this.tauOmega.getTauWatsonMap().get(j).entrySet()) {
                sum += e.getValue();
            }
            startMap.put(j, sum);
        }
        Frequency<Integer> startF = new Frequency<>(startMap);

        if (Globals.getINSTANCE().isSAMPLE_READS()) {
            System.out.println("");
            double counter = 0;
            for (int i = 0; i < amount; i += Globals.getINSTANCE().getSTEPS()) {
                int b = i + Globals.getINSTANCE().getSTEPS();
                if (b >= amount) {
                    b = amount;
                }
                readFutures.add(Threading.getINSTANCE().getExecutor().submit(new SingleModelSampling(or, tauOmega, paired, i, b, startF)));
                counter += Globals.getINSTANCE().getSTEPS();
                StatusUpdate.getINSTANCE().print("Sampling Reads\t\t" + (Math.round((counter / amount) * 100)) + "%");
            }
            int x = 0;
            for (Future<List<SampledRead>> f : readFutures) {
                try {
                    List<SampledRead> srList = f.get();
                    for (SampledRead sr : srList) {
                        startStopSB.append(sr.getWatsonStart()).append("\t");
                        startStopSB.append(sr.getWatsonEnd());

                        readBuilder.append(">read").append(x).append("\n");
                        for (byte r : sr.getWatsonReadBases()) {
                            readBuilder.append(reverse(r));
                        }
                        readBuilder.append("\n");

                        for (int j = sr.getWatsonStart(); j < sr.getWatsonEnd(); j++) {
                            coverage[j]++;
                        }

                        if (paired && sr.getCrickStart() > 0) {
                            startStopSB.append("\t").append(sr.getCrickStart()).append("\t");
                            startStopSB.append(sr.getCrickEnd()).append("\n");

                            readBuilder.append(">read").append(x).append("\n");
                            for (byte r : sr.getCrickReadBases()) {
                                readBuilder.append(reverse(r));
                            }
                            readBuilder.append("\n");

                            for (int j = sr.getCrickStart(); j < sr.getCrickEnd(); j++) {
                                coverage[j]++;
                            }
                        } else {
                            startStopSB.append("\n");
                        }
                    }
                } catch (InterruptedException | ExecutionException ex) {
                    Logger.getLogger(ModelSampling.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            Utils.saveFile(savePath + "sampledReads.fasta", readBuilder.toString());
            Utils.saveFile(savePath + "support" + File.separator + "startStop.txt", startStopSB.toString());
            StringBuilder coverageSB = new StringBuilder();
            coverageSB.append("x\ty\n");
            for (int i = 0; i < L; i++) {
                coverageSB.append(i).append("\t").append(coverage[i]).append("\n");
            }
            Utils.saveFile(savePath + "support" + File.separator + "simCov.txt", coverageSB.toString());
        }
        System.out.println("");

        List<Future<byte[]>> readFuturesFullLength = Lists.newArrayListWithExpectedSize(amount);
        double counterFullLength = 0;
        for (int i = 0; i < amount; i++) {
            readFuturesFullLength.add(Threading.getINSTANCE().getExecutor().submit(new HaplotypeSampling(or)));
            StatusUpdate.getINSTANCE().print("Sampling Haplotypes\t" + (Math.round((counterFullLength++ / amount) * 100)) + "%");
        }
        int y = 0;
        for (Future<byte[]> f : readFuturesFullLength) {
            try {
                readArray[y++] = f.get();
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(ModelSampling.class.getName()).log(Level.SEVERE, null, ex);
            }
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
            double f = ((double) reads.get(read)) / amount;
            if (f > Globals.getINSTANCE().getCUTOFF()) {
                sb.append(">read").append(i++).append("_").append(f).append("\n");
                StringBuilder sb2 = new StringBuilder();
                for (int r : read) {
                    sb2.append(reverse(r));
                }
                map.put(sb2.toString(), f);
                sb.append(sb2);
                sb.append("\n");
            }
        }
        if (Globals.getINSTANCE().isSAMPLE_PROTEINS()) {
            saveProteinHaplotypes(0);
            saveProteinHaplotypes(1);
            saveProteinHaplotypes(2);
        }
    }

    private void saveProteinHaplotypes(int frame) {
        Map<String, Double> proteins = new LinkedHashMap<>();
        for (Map.Entry<String, Double> e : this.map.entrySet()) {
            String p = dna2protein(e.getKey(), frame);
            if (proteins.containsKey(p)) {
                proteins.put(p, proteins.get(p) + e.getValue());
            } else {
                proteins.put(p, e.getValue());
            }
        }
        int i = 0;
        StringBuilder sbp = new StringBuilder();
        for (Object o : sortMapByValue(proteins).keySet()) {
            String read = (String) o;
            if (read.contains("?")) {
                continue;
            }
            double f = ((double) proteins.get(read));
            if (f > Globals.getINSTANCE().getCUTOFF()) {
                sbp.append(">read").append(i++).append("_").append(f).append("\n");
                sbp.append(read);
                sbp.append("\n");
            }
        }
        Utils.saveFile(savePath + "quasispecies_protein_" + frame + ".fasta", sbp.toString());
    }

    public Map<String, Double> getMap() {
        return map;
    }

    public void saveQuasispeciesOnly(String path) {
        Utils.saveFile(path, sb.toString());
    }

    public void save() {
        Utils.saveFile(savePath + "quasispecies.fasta", sb.toString());
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

    public int getK() {
        return K;
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

    public static String dna2protein(String dna, int frame) {
        StringBuilder protein = new StringBuilder();
        for (int i = frame; i < dna.length(); i += 3) {
            if (i + 3 >= dna.length()) {
                break;
            }
            String codon = dna.substring(i, i + 3);
            switch (codon) {
                case "GCT":
                case "GCC":
                case "GCA":
                case "GCG":
                    protein.append("A");
                    break;
                case "TTA":
                case "TTG":
                case "CTT":
                case "CTC":
                case "CTA":
                case "CTG":
                    protein.append("L");
                    break;
                case "CGT":
                case "CGC":
                case "CGA":
                case "CGG":
                case "AGA":
                case "AGG":
                    protein.append("R");
                    break;
                case "AAA":
                case "AAG":
                    protein.append("K");
                    break;
                case "AAT":
                case "AAC":
                    protein.append("N");
                    break;
                case "ATG":
                    protein.append("M");
                    break;
                case "GAT":
                case "GAC":
                    protein.append("D");
                    break;
                case "TTT":
                case "TTC":
                    protein.append("F");
                    break;
                case "TGT":
                case "TGC":
                    protein.append("C");
                    break;
                case "CCT":
                case "CCC":
                case "CCA":
                case "CCG":
                    protein.append("P");
                    break;
                case "CAA":
                case "CAG":
                    protein.append("Q");
                    break;
                case "TCT":
                case "TCC":
                case "TCA":
                case "TCG":
                case "AGT":
                case "AGC":
                    protein.append("S");
                    break;
                case "GAA":
                case "GAG":
                    protein.append("E");
                    break;
                case "ACT":
                case "ACC":
                case "ACA":
                case "ACG":
                    protein.append("T");
                    break;
                case "GGT":
                case "GGC":
                case "GGA":
                case "GGG":
                    protein.append("G");
                    break;
                case "TGG":
                    protein.append("W");
                    break;
                case "CAT":
                case "CAC":
                    protein.append("H");
                    break;
                case "TAT":
                case "TAC":
                    protein.append("Y");
                    break;
                case "ATT":
                case "ATC":
                case "ATA":
                    protein.append("I");
                    break;
                case "GTT":
                case "GTC":
                case "GTA":
                case "GTG":
                    protein.append("V");
                    break;
                case "TAA":
                case "TGA":
                case "TAG":
                    protein.append("-");
                    break;
                default:
                    protein.append("?");
                    break;
            }
        }
        return protein.toString();
    }
}

class ValueComparator implements Comparator, Serializable {

    private static final long serialVersionUID = 1L;
    Map base;

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
