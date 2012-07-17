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

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorker;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorkerRecalc;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;

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
    //rho[j][k][l] := transition prob. at position j, for l given k
    private double[][][] rho;
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
    private Read[] reads;
    private ReadHMM[] readHMMArray;
    private int restart = 0;
    private int readCount;

    public JHMM(Read[] reads, int N, int L, int K, int n, double epsilon) {
        this(reads, N, L, K, n, epsilon,
                (1 - (n - 1) * epsilon),
                Random.generateInitRho(L - 1, K),
                Random.generateInitPi(K),
                Random.generateMuInit(L, K, n));
    }

    private JHMM(Read[] reads, int N, int L, int K, int n, double eps, double antieps, double[][][] rho, double[] pi, double[][][] mu) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.reads = reads;
        this.rho = rho;
        this.mu = mu;
        this.pi = pi;

        this.eps = new double[L];
        this.antieps = new double[L];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                this.eps[j] = eps;
                this.antieps[j] = antieps;
            }
        }

        this.start();
        this.calculate();
    }

    private void calculateLoglikelihood() {
        this.loglikelihood = 0d;
        for (ReadHMM r : this.readHMMArray) {
            for (int j = 0; j < L; j++) {
                this.loglikelihood += Math.log(r.getC(j)) * r.getCount();
            }
        }
    }

    private void start() {
        long time = System.currentTimeMillis();
        readHMMArray = new ReadHMM[reads.length];
        this.readCount = this.reads.length;
        if (Globals.PARALLEL_JHMM) {
            readHMMArray = Globals.fjPool.invoke(new ReadHMMWorker(this, reads, rho, pi, mu, eps, antieps, K, L, n, 0, this.readCount)).toArray(new ReadHMM[this.readCount]);
        } else {
            for (int i = 0; i < reads.length; i++) {
                readHMMArray[i] = new ReadHMM(L, K, n, reads[i], rho, pi, mu, eps, antieps);
            }
        }
        System.out.println("R\t: " + (System.currentTimeMillis() - time));
    }

    public void restart() {
        this.restart++;
        long time = System.currentTimeMillis();
        if (Globals.PARALLEL_JHMM) {
            Globals.fjPool.invoke(new ReadHMMWorkerRecalc(this.readHMMArray, rho, pi, mu, eps, antieps, 0, this.readCount));
        } else {
            for (ReadHMM r : this.readHMMArray) {
                r.recalc(rho, pi, mu, eps, antieps);
            }
        }
        System.out.println("R\t: " + (System.currentTimeMillis() - time));
        this.calculate();
    }

    private void calculate() {
        this.eStep();
        this.calculateLoglikelihood();
        this.mStep();
    }

    private void eStep() {
        long time = System.currentTimeMillis();
        this.nJK = new double[L][K];
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nVB = new double[n][n];
        this.nJeq = new double[L];
        this.nJneq = new double[L];
        this.neqPos = new double[L];
        this.nneqPos = new double[L];
        for (ReadHMM r : this.readHMMArray) {
            int times = r.getCount();
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
                    byte b = r.getSequence()[j];
                    for (int k = 0; k < K; k++) {
                        if (v != b) {
                            this.nneqPos[j] += r.gamma(j, k, v) * times;
                        } else {
                            this.neqPos[j] += r.gamma(j, k, v) * times;
                        }

                    }
                }
            }
        }
        System.out.println("E\t: " + (System.currentTimeMillis() - time));
    }

    private void calcMu() {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double AH = Globals.ALPHA_H;
                double divisor;
                double sum;
                do {
                    sum = 0d;
                    divisor = 0d;
                    for (int v = 0; v < n; v++) {
                        mu[j][k][v] = this.f(this.getnJKV(j, k, v) + AH);
                        sum += this.getnJKV(j, k, v);
                    }
                    sum = this.f(sum + n * AH);
                    if (sum > 0) {
                        for (int v = 0; v < n; v++) {
                            mu[j][k][v] /= sum;
                            divisor += mu[j][k][v];
                        }
                        if (divisor > 0) {
                            for (int v = 0; v < n; v++) {
                                mu[j][k][v] /= divisor;
                            }
                        } else {
                            AH *= 10;
                        }
                    } else {
                        AH *= 10;
                    }
                } while (sum == 0 || divisor == 0);
            }
        }
    }

    private void calcRho() {
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double AZ = Globals.ALPHA_Z;
                double divisor;
                double sum;
                do {
                    sum = 0d;
                    divisor = 0d;
                    for (int l = 0; l < K; l++) {
                        rho[j - 1][k][l] = this.f(this.getnJKL(j, k, l) + AZ);
                        sum += this.getnJKL(j, k, l);
                    }
                    sum = this.f(sum + K * AZ);
                    if (sum > 0) {
                        for (int l = 0; l < K; l++) {
                            rho[j - 1][k][l] /= sum;
                            divisor += rho[j - 1][k][l];
                        }
                        if (divisor > 0) {
                            for (int l = 0; l < K; l++) {
                                rho[j - 1][k][l] /= divisor;
                            }
                        } else {
                            AZ *= 10;
                        }
                    } else {
                        AZ *= 10;
                    }
                } while (sum == 0 || divisor == 0);
            }
        }
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
        double sumK = 0d;
        for (int k = 0; k < K; k++) {
            pi[k] = this.getnJK(0, k);
            sumK += pi[k];
        }
        for (int k = 0; k < K; k++) {
            pi[k] /= sumK;
        }
        return pi;
    }

    private void mStep() {
        if (!Globals.NO_RECOMB) {
            this.calcRho();
        }
        this.calcPi();
        this.calcMu();



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

//        if (!Globals.FIX_EPSILON) {
//            for (int j = 0; j < L; j++) {
//                this.eps[j] = this.nneqPos[j] / (N * (n - 1));
//            }
//        } else {
        for (int j = 0; j < L; j++) {
            this.eps[j] = f(ew + a) / f((N * (n - 1)) + a + b);
            this.antieps[j] = 1 - (n - 1) * eps[j];
        }
//    }
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

    public ReadHMM[] getReadHMMArray() {
        return this.readHMMArray;
    }

    public double getLikelihood() {
        return likelihood;
    }

    public int getRestart() {
        return restart;
    }
}
