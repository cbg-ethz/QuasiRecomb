package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.RestartWorker;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Responsible for the start of the repeats within the EM algorithm.
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class EM extends Utils {

    private StringBuilder sb = new StringBuilder();
    private double BIC_opt = Double.NEGATIVE_INFINITY;
    private int iterations = 0;
    private int N;
    private int K;
    private int L;
    private int n;
    private Map<byte[], Integer> reads;
    private byte[][] haplotypesArray;
    private OptimalResult or;

    protected EM(int N, int L, int K, int n, String[] reads, byte[][] haplotypesArray) {
        this(K, N, L, n, Utils.clusterReads(Utils.splitReadsIntoByteArrays(reads)), haplotypesArray);
    }

    protected EM(int N, int L, int K, int n, Map<byte[], Integer> reads, byte[][] haplotypesArray) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.reads = reads;
        this.haplotypesArray = haplotypesArray;
        this.blackbox();
    }

    protected EM(int N, int L, int K, int n, byte[][] reads, byte[][] haplotypesArray) {
        this(K, N, L, n, Utils.clusterReads(reads), haplotypesArray);
    }

    public void shorten(double value) {
        if (value < 1e-5) {
            sb.append("0      ");
        } else if (value == 1.0) {
            sb.append("1      ");
        } else {
            String t = "" + value;
            String r;
            if (t.length() > 5) {
                r = t.substring(0, 7);
                if (t.contains("E")) {
                    r = r.substring(0, 4);
                    r += "E" + t.split("E")[1];
                }
                sb.append(r);
            } else {
                sb.append(value);
            }
        }
    }

    private void blackbox() {
        List<OptimalResult> ors;
        double maxLLH = Double.NEGATIVE_INFINITY;
//        boolean givenMinLLH = !Double.isInfinite(Globals.MIN_LLH);
//        if (!givenMinLLH) {
//            if (Globals.PARALLEL_RESTARTS) {
//                ors = Globals.fjPool.invoke(new RestartWorker(N, K, L, n, reads, haplotypesArray, Globals.FILTER_LLH, 0, 1));
//            } else {
//                ors = new LinkedList<>();
//                for (int i = 0; i < 1; i++) {
//                    ors.add(new SingleEM(N, K, L, n, reads, haplotypesArray, Globals.FILTER_LLH).getOptimalResult());
//                }
//            }
//
//            for (OptimalResult tmp : ors) {
//                maxLLH = Math.max(maxLLH, tmp.getLlh());
//            }
//            Globals.maxMAX_LLH(maxLLH);
//        } else {
//            Globals.maxMAX_LLH(Globals.MIN_LLH);
//        }
        System.out.println("--------------------");
        ors = Globals.fjPool.invoke(new RestartWorker(N, K, L, n, reads, haplotypesArray, Globals.DELTA_LLH, 0, Globals.REPEATS));
        Globals.PARALLEL_JHMM = true;
        double maxBIC = Double.NEGATIVE_INFINITY;
        for (OptimalResult tmp : ors) {
            if (tmp.getBIC() > maxBIC) {
                maxBIC = tmp.getBIC();
                or = tmp;
            }
        }
        if (!Globals.NO_REFINE) {
            SingleEM bestEM = new SingleEM(N, K, L, n, reads, haplotypesArray, 1e-10, or);
            this.or = bestEM.getOptimalResult();
        }
        this.saveEM();
    }

    private void saveEM() {
        this.saveBestEM(false);
    }

    public void saveBestEM(boolean bestK) {
        if (Globals.DEBUG) {
            System.out.println("BIC:" + or.getBIC());
        }
        sb.setLength(0);
        sb.append("#loglikelihood:").append(or.getLlh()).append("\n");
        sb.append("#iterations:").append(iterations).append("\n");
        sb.append("#BIC:").append(or.getBIC()).append("\n");
        sb.append("#EPS:").append(or.getEps()).append("\n");
        sb.append("#PI:\n");
        for (int k = 0; k < or.getK(); k++) {
            sb.append("##").append(or.getPi()[k]).append("\n");
        }
        sb.append("#RHO PRIOR:\n");
        for (int j = 0; j < or.getL() - 1; j++) {
            sb.append(j + 1).append("P\t");
            for (int k = 0; k < K; k++) {
                sb.append(Arrays.toString(or.getPriorRho()[j][k]));
                sb.append("\t");
            }
            sb.append("\n");
        }
        sb.append("#RHO:\n");
        for (int j = 0; j < or.getL() - 1; j++) {
            sb.append(j + 1).append("\t");
            for (int k = 0; k < K; k++) {
                sb.append("[");
                for (int l = 0; l < K; l++) {
                    shorten(or.getRho()[j][k][l]);
                    if (l + 1 < K) {
                        sb.append(", ");
                    }

                }
                sb.append("]\t");
            }
            sb.append("\n");
        }
        sb.append("#MU:\n");
        for (int j = 0; j < or.getL(); j++) {
            sb.append("##j:").append(j).append("\t");
            for (byte[] b : haplotypesArray) {
                sb.append(reverse((int) b[j]));
            }
            sb.append("|");
            for (int k = 0; k < or.getK(); k++) {
                double max = Double.MIN_VALUE;
                Map<Double, Integer> m = new HashMap<>();
                for (int v = 0; v < or.getMu()[0][0].length; v++) {
                    max = Math.max(max, or.getMu()[j][k][v]);
                    m.put(or.getMu()[j][k][v], v);
                }
                sb.append(reverse(m.get(max))).append("-");
            }
            for (int k = 0; k < or.getK(); k++) {
                sb.append("[");
                for (int v = 0; v < n; v++) {
                    shorten(or.getMu()[j][k][v]);
                    if (v + 1 < n) {
                        sb.append(", ");
                    }

                }
                sb.append("]\t");
            }
            sb.append("\n");
        }
        try {
            String s = Globals.savePath + "optimumJavaK" + (bestK ? "" : K);
            FileOutputStream fos = new FileOutputStream(s);
            try (ObjectOutputStream out = new ObjectOutputStream(fos)) {
                out.writeObject(or);
            }
        } catch (IOException ex) {
            System.out.println("Optimum Java saving\n" + ex.getMessage());
        }
    }

    /**
     * Returns the StringBuilder with the information of the best result.
     *
     * @return StringBuilder with result
     */
    public StringBuilder getSb() {
        return sb;
    }

    /**
     * Returns the best BIC.
     *
     * @return best BIC
     */
    public double getBIC_opt() {
        return BIC_opt;
    }

    /**
     * The double 3D array of mu for the best result.
     *
     * @return mu
     */
    public double[][][] getMu_opt() {
        return or.getMu();
    }

    /**
     * The double array of pi for the best result.
     *
     * @return pi
     */
    public double[] getPi_opt() {
        return or.getPi();
    }

    /**
     * The double 3D array of rho for the best result.
     *
     * @return rho
     */
    public double[][][] getRho_opt() {
        return or.getRho();
    }
}
