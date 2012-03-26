package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.Map;

/**
 * Selects the best model among the specified range of generators.
 * 
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class ModelSelection {

    private int Kmin;
    private int Kmax;
    private int N;
    private int L;
    private int n;
    private int bestK;
    private Map<byte[], Integer> clusterReads;
    private byte[][] haplotypesArray;
    private double[][][] mu = null;
    private double[][][] rho = null;
    private double[] pi = null;

    public ModelSelection(Map<byte[], Integer> clusterReads, int Kmin, int Kmax, int N, int L, int n, byte[][] haplotypesArray) {
        this.Kmax = Kmax;
        this.Kmin = Kmin;
        this.N = N;
        this.L = L;
        this.n = n;
        this.clusterReads = clusterReads;
        this.haplotypesArray = haplotypesArray;
        this.start();
    }

    private void start() {
        double optBIC = 0;
        for (int k = Kmin; k <= Kmax; k++) {
            if (!Globals.rho0force || k == 1) {
                checkRho0(k);
            }
            EM em = new EM(this.N, this.L, k, this.n, this.clusterReads, this.haplotypesArray);
            Utils.saveFile(Globals.savePath + "K" + k + "-result.txt", em.getSb().toString());
            if (em.getBIC_opt() > optBIC || optBIC == 0) {
                optBIC = em.getBIC_opt();
                this.bestK = k;
                em.saveBestEM(true);
                this.mu = em.getMu_opt();
                this.pi = em.getPi_opt();
                this.rho = em.getRho_opt();
            }
            ModelSampling ms = new ModelSampling(L, n, k, em.getRho_opt(), em.getPi_opt(), em.getMu_opt(), Globals.savePath);
            Map<byte[], Integer> reads = ms.getReads();
            Double result = 0.0;
//                for (byte base : map.keySet()) {
//                    Double frequency = (double) map.get(base) / l;
//                    result -= frequency * Math.log(frequency);
//                }
        }
    }

    private static void checkRho0(int K) {
        if (K == 1) {
            Globals.rho0 = true;
        } else {
            Globals.rho0 = false;
        }
    }

    public int getBestK() {
        return bestK;
    }

    public double[][][] getMu() {
        return mu;
    }

    public double[] getPi() {
        return pi;
    }

    public double[][][] getRho() {
        return rho;
    }
}
