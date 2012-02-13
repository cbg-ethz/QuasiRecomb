package ch.ethz.bsse.quasirecomb.modelsampling;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public final class ModelSampling extends Utils {

    private String savePath = "";
    private int L;
    private int n;
    private final int K;
    private double[][][] rho;
    private double[] pi;
    private double[][][] H;
    private Frequency<Integer>[][] rhoArray;
    private Frequency<Byte>[][] muArray;
    private Map<byte[], Integer> reads;
    private int amount = 10000;
    private int[] recombPerObservation;
    private Map<String, Double> hexMap = new HashMap<>();

    public ModelSampling(int L, int n, int K, double[][][] rho, double[] pi, double[][][] mu, String savePath) {
        this.L = L;
        this.n = n;
        this.K = K;
        this.rho = rho;
        this.pi = pi;
        this.H = mu;
        this.rhoArray = new Frequency[L - 1][K];
        this.muArray = new Frequency[L][K];
        this.savePath = savePath;
        this.recombPerObservation = new int[amount];
        this.start();
    }

    public ModelSampling(String string, String path, int amount) {
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
        this.reads = or.getReads();
        this.pi = or.getPi();
        this.H = or.getMu();
        this.rhoArray = new Frequency[L - 1][K];
        this.muArray = new Frequency[L][K];
        this.amount = amount;
        this.recombPerObservation = new int[amount];
        this.start();
    }

    public void start() {
        new File(savePath).mkdirs();
//        saveFile(savePath + "simu-K-" + K + ".txt", "K:" + K);
        if (new File(savePath + "K" + K + "-result.txt").exists()) {
            new File(savePath + "K" + K + "-result.txt").renameTo(new File(savePath + "K" + K + "-opt.txt"));
        }
        this.reads = new HashMap<>();
        ValueComparator bvc = new ValueComparator(reads);
        TreeMap<byte[], Integer> sorted_map = new TreeMap(bvc);
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

        sorted_map.putAll(reads);
        StringBuilder sb = new StringBuilder();
        int sum = 0;
        for (byte[] read : reads.keySet()) {
            sb.append(((double) reads.get(read))).append("\t");
            sum += ((double) reads.get(read));
            for (int r : read) {
                sb.append(reverse(r));
            }
            sb.append("\n");
        }
        Utils.saveFile(savePath + "haplotypeEnumNoSort-K" + K + ".txt", sb.toString());
        Utils.saveFile(savePath + "haplotypeAmount-K" + K + ".txt", "" + sum);

        sb.setLength(0);
        for (byte[] read : sorted_map.keySet()) {
            sb.append(((double) reads.get(read)) / amount).append("\t");
            for (int r : read) {
                sb.append(reverse(r));
            }
            sb.append("\n");
        }
        Utils.saveFile(savePath + "haplotypeDist-K" + K + ".txt", sb.toString());
        sb.setLength(0);


        for (byte[] read : reads.keySet()) {
            sb.append(">").append(((double) reads.get(read)) / amount).append("\n");
            for (int r : read) {
                sb.append(reverse(r));
            }
            sb.append("\n");
        }
        Utils.saveFile(savePath + "haplotypeDistNoSort-K" + K + ".fasta", sb.toString());
        sb.setLength(0);

        for (byte[] read : sorted_map.keySet()) {
            sb.append(">").append(((double) reads.get(read)) / amount).append("\n");
            for (int r : read) {
                sb.append(reverse(r));
            }
            sb.append("\n");
        }
        Utils.saveFile(savePath + "haplotypeDist-K" + K + ".fasta", sb.toString());
        sb.setLength(0);
        for (byte[] read : sorted_map.keySet()) {
            sb.append(((double) reads.get(read)) / amount).append("\t");
            try {
                MessageDigest digest = java.security.MessageDigest.getInstance("MD5");
                digest.update(read);
                String result = "";
                byte[] b = digest.digest();
                for (int i = 0; i < b.length; i++) {
                    result +=
                            Integer.toString((b[i] & 0xff) + 0x100, 16).substring(1);
                }
                sb.append(result);
                hexMap.put(result, ((double) reads.get(read)) / amount);
            } catch (NoSuchAlgorithmException ex) {
                Logger.getLogger(ModelSampling.class.getName()).log(Level.SEVERE, null, ex);
            }
            sb.append("\n");
        }
        Utils.saveFile(savePath + "haplotypeDist-K" + K + ".hex", sb.toString());
//        sb.setLength(0);
//        int sum = 0;
//        int hits = 0;
//        for (int i = 0; i < amount; i++) {
//            sb.append("#").append(i).append(":").append(this.recombPerObservation[i]).append("\n");
//            if (this.recombPerObservation[i] > 0) {
//                hits++;
//            }
//            sum += this.recombPerObservation[i];
//        }
//        sb.append("#sum:").append(sum).append("\n");
//        sb.append("#p:").append((double) sum / amount / L).append("\n");
//        sb.append("#recombinated observations:").append(hits).append("\n");
//        Utils.saveFile(savePath + "simuRecombRate.txt", sb.toString());
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
//                if (Math.random() > (1.0 - (K - 1.0) * rho[j - 1])) {
                //TODO: Test if that works with new rho approach
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
//                }
            }
            if (muArray[j][k] == null) {
                Map<Byte, Double> muMap = new HashMap<>();
                for (byte v = 0; v < n; v++) {
                    muMap.put(v, H[j][k][v]);
                }
                Frequency<Byte> muF = new Frequency<>(muMap);
                muArray[j][k] = muF;
            }
            read[j] = muArray[j][k].roll();
        }
        return read;
    }

    public Map<byte[], Integer> getReads() {
        return reads;
    }

    public Map<String, Double> getHexMap() {
        return hexMap;
    }
}

class ValueComparator implements Comparator {

    Map base;

    public ValueComparator(Map base) {
        this.base = base;
    }

    @Override
    public int compare(Object a, Object b) {

        if ((Integer) base.get(a) < (Integer) base.get(b)) {
            return 1;
        } else if ((Integer) base.get(a) == (Integer) base.get(b)) {
            return 0;
        } else {
            return -1;
        }
    }
}