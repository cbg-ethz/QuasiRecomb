package ch.ethz.bsse.quasirecomb.informatioholder;

import java.io.Serializable;
import java.util.Map;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class OptimalResult implements Serializable {

    private int N;
    private int K;
    private int L;
    private int n;
    private Map<byte[], Integer> reads;
    private byte[][] haplotypesArray;
    private double[][][] rho;
    private double[][][] priorRho;
    private double[] pi;
    private double[][][] mu;
    private double llh;
    private double[] eps;
    private double BIC;
    private boolean [][] nneqPosCount;

    public OptimalResult(int N, int K, int L, int n, Map<byte[], Integer> reads, byte[][] haplotypesArray, double[][][] rho, double[] pi, double[][][] mu, double llh, double BIC, double[][][] priorRho, double[] eps, boolean[][] nneqPosCount) {
        this.N = N;
        this.K = K;
        this.L = L;
        this.n = n;
        this.reads = reads;
        this.haplotypesArray = haplotypesArray;
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.llh = llh;
        this.BIC = BIC;
        this.priorRho = priorRho;
        this.eps = eps;
        this.nneqPosCount = nneqPosCount;
    }

    public int getK() {
        return K;
    }

    public int getL() {
        return L;
    }

    public int getN() {
        return N;
    }

    public byte[][] getHaplotypesArray() {
        return haplotypesArray;
    }

    public double getLlh() {
        return llh;
    }

    public double[][][] getMu() {
        return mu;
    }

    public int getn() {
        return n;
    }

    public double[] getPi() {
        return pi;
    }

    public Map<byte[], Integer> getReads() {
        return reads;
    }

    public double[][][] getRho() {
        return rho;
    }

    public double getBIC() {
        return BIC;
    }

    public double[][][] getPriorRho() {
        return priorRho;
    }

    public double[] getEps() {
        return eps;
    }

    public boolean[][] getNneqPosCount() {
        return nneqPosCount;
    }
}
