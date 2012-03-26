package ch.ethz.bsse.quasirecomb.modelsampling;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;


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
    private StringBuilder sb = new StringBuilder();
    private TreeMap<byte[], Integer> sorted_map;

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
        this.reads = new HashMap<>();
        ValueComparator bvc = new ValueComparator(reads);
        this.sorted_map = new TreeMap(bvc);
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

//        sorted_map.putAll(reads);
        int i = 0;
//        for (byte[] read : sorted_map.keySet()) {
        for (Object o : sortMapByValue(reads).keySet()) {
            byte[] read = (byte[]) o;
            sb.append(">read").append(i++).append("_").append(((double) reads.get(read)) / amount).append("\n");
            for (int r : read) {
                sb.append(reverse(r));
            }
            sb.append("\n");
        }
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
                    muMap.put(v, H[j][k][v]);
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
        int i = 0;
        for (byte[] read : sorted_map.keySet()) {
            System.out.println(reads.get(read) + "\t" + read);
        }
    }

    public Map<byte[], Integer> getReads() {
        return reads;
    }

    public static Map sortMapByValue(Map map) {



        List listForSort = null;

        Map sortedList = new LinkedHashMap();

        listForSort = new LinkedList(map.entrySet());

        Collections.sort(listForSort, new Comparator() {
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