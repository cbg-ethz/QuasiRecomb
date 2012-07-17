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
package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorker;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorkerRecalc;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Offers a HMM with the E and M step.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class JHMM {

    private final int N;
    private final int L;
    private final int K;
    private final int n;
    //rho[j][l][k] := transition prob. at position j, for k given l
    private double[][][] rho;
    private double[][][] priorRho;
    private double[] pi;
    private double[][][] mu;
    private double[] eps;
    private double[] antieps;
    private double[][] nJK;
    private double[][][] nJKL;
    private double[][][] nJKV;
    private double[][] nVB;
    private double[] neqPos;
    private double[] nneqPos;
    private double[] nJeq;
    private double[] nJneq;
    private double neq;
    private double nneq;
    private double loglikelihood;
    private double likelihood;
    private double Q;
    private Map<ReadHMM, Integer> readHMMMap;
    private int mapSize;
    private Map<byte[], Integer> clusterReads;
    private int restart = 0;

    public JHMM(String[] reads, int K, int n) {
        this(Utils.splitReadsIntoByteArrays(reads), K, n, Globals.ESTIMATION_EPSILON);
    }

    public JHMM(String[] reads, int K, int n, double epsilon) {
        this(Utils.splitReadsIntoByteArrays(reads), K, n, epsilon);
    }

    public JHMM(byte[][] reads, int K, int n) {
        this(reads, K, n, Globals.ESTIMATION_EPSILON);
    }

    public JHMM(byte[][] reads, int K, int n, double epsilon) {
        this(reads, reads.length, reads[0].length, K, n, epsilon);
    }

    public JHMM(OptimalResult or) {
        throw new IllegalAccessError("JHMM");
//        this.N = or.getN();
//        this.L = or.getL();
//        this.K = or.getK();
//        this.n = or.getn();
//        this.clusterReads = or.getReads();
//        this.rho = or.getRho();
//        this.mu = or.getMu();
//        this.eps = Globals.ESTIMATION_EPSILON;
//        this.antieps = (1 - (n - 1) * Globals.ESTIMATION_EPSILON);
//        this.pi = or.getPi();
//        this.priorRho = or.getPriorRho();
//        this.start();
//        this.calculate();
    }

    public JHMM(Map<byte[], Integer> reads, int N, int L, int K, int n, double epsilon) {
        this(reads, N, L, K, n, epsilon,
                (1 - (n - 1) * epsilon),
                Random.generateInitRho(L - 1, K),
                Random.generateInitPi(K),
                Random.generateMuInit(L, K, n));
    }

    public JHMM(byte[][] reads, int N, int L, int K, int n, double epsilon) {
        this(reads, N, L, K, n, epsilon,
                (1 - (n - 1) * epsilon),
                Random.generateInitRho(L - 1, K),
                Random.generateInitPi(K),
                Random.generateMuInit(L, K, n));
    }

    public JHMM(byte[][] reads, int N, int L, int K, int n, double eps, double antieps, double[][][] rho, double[] pi, double[][][] mu) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.clusterReads = Utils.clusterReads(reads);
        this.rho = rho;
        this.mu = mu;
        this.uniformEpsilon(eps, antieps);
        this.pi = pi;
        this.priorRho = rho;
        this.start();
        this.calculate();
    }

    private JHMM(Map<byte[], Integer> reads, int N, int L, int K, int n, double eps, double antieps, double[][][] rho, double[] pi, double[][][] mu) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.clusterReads = reads;
        this.rho = rho;
        this.mu = mu;
        this.uniformEpsilon(eps, antieps);
        this.pi = pi;
        this.priorRho = rho;
        this.start();
        this.calculate();
    }

    private void uniformEpsilon(double eps, double antieps) {
        this.eps = new double[L];
        this.antieps = new double[L];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {

                this.eps[j] = eps;
                this.antieps[j] = antieps;
            }
        }
    }

    private void calculateLoglikelihood() {
        this.loglikelihood = 0d;
        for (ReadHMM r : this.readHMMMap.keySet()) {
            for (int j = 0; j < L; j++) {
                this.loglikelihood += Math.log(r.getC(j)) * this.readHMMMap.get(r);
            }
        }
    }

    private void start() {
//        long time = System.currentTimeMillis();
        byte[][] uniqueReads = new byte[clusterReads.keySet().size()][L];
        for (byte[] read : this.clusterReads.keySet()) {
            uniqueReads[mapSize++] = read;
        }
        this.readHMMMap = new HashMap<>();
        ReadHMM[] rArray = new ReadHMM[uniqueReads.length];
        if (Globals.PARALLEL_JHMM) {
            rArray = Globals.fjPool.invoke(new ReadHMMWorker(this, uniqueReads, rho, pi, mu, eps, antieps, K, L, n, 0, mapSize)).toArray(new ReadHMM[mapSize]);
        } else {
            for (int i = 0; i < uniqueReads.length; i++) {
                rArray[i] = new ReadHMM(L, K, n, uniqueReads[i], rho, pi, mu, eps, antieps);
            }
        }
        for (ReadHMM r : rArray) {
            this.readHMMMap.put(r, clusterReads.get(r.getRead()));
        }
//        System.out.println("R\t: " + (System.currentTimeMillis() - time));
    }

    public void restart() {
        this.restart++;
//        long time = System.currentTimeMillis();
        if (Globals.PARALLEL_JHMM) {
            Globals.fjPool.invoke(new ReadHMMWorkerRecalc(this.readHMMMap.keySet().toArray(new ReadHMM[this.readHMMMap.keySet().size()]), rho, pi, mu, eps, 0, mapSize));
        } else {
            for (ReadHMM r : this.readHMMMap.keySet()) {
                r.recalc(rho, pi, mu, eps);
            }
        }
//        System.out.println("R\t: " + (System.currentTimeMillis() - time));
        this.calculate();
    }

    private void calculate() {
        this.eStep();
        this.calculateLoglikelihood();
        this.mStep();
    }

    private void eStep() {
//        long time = System.currentTimeMillis();
        this.nJK = new double[L][K];
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nVB = new double[n][n];
        this.nJeq = new double[L];
        this.nJneq = new double[L];
        this.neqPos = new double[L];
        this.nneqPos = new double[L];
        for (Iterator<ReadHMM> it = this.readHMMMap.keySet().iterator(); it.hasNext();) {
            ReadHMM r = it.next();
            int times = this.readHMMMap.get(r);
            for (int j = 0; j < L; j++) {
                for (int k = 0; k < K; k++) {
                    this.nJK[j][k] += r.gamma(j, k) * times;
                    if (j > 0) {
                        for (int l = 0; l < K; l++) {
                            this.nJKL[j][k][l] += r.xi(j, k, l) * times;
                            if (k == l) {
                                this.nJeq[j] += r.xi(j, k, l) * times;
                            } else {
                                this.nJneq[j] += r.xi(j, k, l) * times;
                            }
                        }
                    }
                    for (int v = 0; v < n; v++) {
                        this.nJKV[j][k][v] += r.gamma(j, k, v) * times;
                    }
                }
                for (int v = 0; v < n; v++) {
                    byte b = r.getRead()[j];
                    for (int k = 0; k < K; k++) {
                        if (v != b) {
                            this.nneqPos[j] += r.gamma(j, k, v) * times;
                        } else {
                            this.neqPos[j] += r.gamma(j, k, v) * times;
                        }

                    }
                }
//        System.out.println("E\t: " + (System.currentTimeMillis() - time));
            }
        }
    }

    private double[][][] calcMu() {
        double[][][] mu_ = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double sumV = 0d;
                double sumP = 0d;
                for (int v = 0; v < n; v++) {
                    mu_[j][k][v] = this.f(this.getnJKV(j, k, v) + Globals.ALPHA_H[j][v]);
                    sumV += this.getnJKV(j, k, v);
                    sumP += Globals.ALPHA_H[j][v];
                }
                sumV = this.f(sumV + sumP);
                double divisor = 0d;
                if (sumV != 0) {
                    for (int v = 0; v < n; v++) {
                        mu_[j][k][v] /= sumV;
                        divisor += mu_[j][k][v];
                    }
                    for (int v = 0; v < n; v++) {
                        mu_[j][k][v] /= divisor;
                    }
                } else {
                    System.out.println("mu");
                }
            }
        }
        return mu_;
    }

    private double[][][] calcRho() {
        double[][][] rho_ = new double[L - 1][K][K];
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double divisor = 0d;
                double sum = 0d;
                for (int l = 0; l < K; l++) {
                    if (Double.isNaN(this.getnJKL(j, k, l))) {
                        System.out.println("rho");
                    }
                    rho_[j - 1][k][l] = this.f(this.getnJKL(j, k, l) + Globals.ALPHA_Z);
                    sum += this.getnJKL(j, k, l);
                }
                sum = this.f(sum + K * Globals.ALPHA_Z);

                for (int l = 0; l < K; l++) {
                    rho_[j - 1][k][l] /= sum;
                    divisor += rho_[j - 1][k][l];
                }
                for (int l = 0; l < K; l++) {
                    rho_[j - 1][k][l] /= divisor;
                }
            }
        }
        return rho_;
    }

    private double f(double upsilon) {
        if (upsilon == 0d) {
            return 0d;
        }
        double result = Math.exp(rho_phi(upsilon));
        if (Double.isNaN(result)) {
            Utils.error();
            System.out.println("XXX:" + upsilon);
            System.exit(9);
        }
        return result;
    }

    private double rho_phi(double upsilon) {
        double x = -1d;
        try {
            x = (upsilon > 7) ? rho_g(upsilon - .5) : (rho_phi(upsilon + 1) - 1 / upsilon);

        } catch (StackOverflowError e) {
            System.err.println(upsilon);
            System.err.println(e);
        }
        return x;
    }

    private double rho_g(double x) {
        return Math.log(x) + .04167 * Math.pow(x, -2) - .00729 * Math.pow(x, -4) + .00384 * Math.pow(x, -6) - .00413 * Math.pow(x, -8);
    }

    private double[] calcPi() {
        double[] pi_ = new double[K];
        double sumK = 0d;
        for (int k = 0; k < K; k++) {
            pi_[k] = this.getnJK(0, k);
            sumK += pi_[k];
        }
        for (int k = 0; k < K; k++) {
            pi_[k] /= sumK;
        }
        return pi_;
    }

    private void mStep() {
        if (!Globals.rho0) {
            this.rho = this.calcRho();
        }
        this.pi = this.calcPi();
        this.mu_old = this.mu;
        this.mu = this.calcMu();
        
        
        
        double ew = 0d;
        int nonzero = 0;
        for (int j = 0; j < L; j++) {
            if (this.nneqPos[j] != 0d) {
                ew += this.nneqPos[j] / N;
                nonzero++;
            }
        }
        ew /= nonzero;
        double a = 20;
        double b = (-a * ew + a + 2 * ew - 1) / ew;
//        if (Globals.FIX_EPSILON) {
            for (int j = 0; j < L; j++) {
                this.eps[j] = f(ew + a) / f((N * (n - 1)) + a + b);
            }
//        } else {
//            for (int j = 0; j < L; j++) {
//                this.eps[j] = this.nneqPos[j] / (N * (n - 1));
//            }
//        }
        
//        if (!Globals.FIX_EPSILON) {
//            for (int j = 0; j < L; j++) {
//                this.eps[j] = this.nneqPos[j] / (N * (n - 1));
//            }
//        }
    }
    public double[][][] mu_old;

    public double getnJK(int j, int k) {
        return nJK[j][k];
    }

    public double getnJKL(int j, int k, int l) {
        return nJKL[j][k][l];
    }

    public double getnJKV(int j, int k, int v) {
        return nJKV[j][k][v];
    }

    public double getnJeq(int j) {
        return nJeq[j];
    }

    public double getnJneq(int j) {
        return nJneq[j];
    }

    public double[][] getnVB() {
        return nVB;
    }

    public double getNeq() {
        return neq;
    }

    public double getNneq() {
        return nneq;
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

    public double[] getAntieps() {
        return antieps;
    }

    public double[] getEps() {
        return eps;
    }

    public int getn() {
        return n;
    }

    public double getLoglikelihood() {
        return loglikelihood;
    }

    public double getQ() {
        return Q;
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

    public Map<ReadHMM, Integer> getReadHMMMap() {
        return readHMMMap;
    }

    public double getLikelihood() {
        return likelihood;
    }

    public double[][][] getPrior_rho() {
        return priorRho;
    }

    public int getRestart() {
        return restart;
    }
}
