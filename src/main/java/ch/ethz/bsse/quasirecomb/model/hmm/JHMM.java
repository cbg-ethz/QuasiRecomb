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

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.EInfo;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorker;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorkerRecalc;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.math.BigDecimal;
import java.util.List;
import org.javatuples.Pair;

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
    private double[][] tauOmega;
    private double[][][] rho;
    private double[] pi;
    private double[][][] mu;
    private double[] eps;
    private double[] antieps;
    private double[][] nJK;
    private double[][][] nJKL;
    private double[][][] nJKV;
    private double[] neqPos;
    private double[] nneqPos;
    private double[] nJeq;
    private double[] nJneq;
    private double loglikelihood;
    private Read[] reads;
    private ReadHMM[] readHMMArray;
    private int restart = 0;
    private int readCount;
    private int coverage[];
    private int parametersChanged = 0;
    private boolean paired;

    public JHMM(Read[] reads, int N, int L, int K, int n, double epsilon) {
        this(reads, N, L, K, n, epsilon,
                (1 - (n - 1) * epsilon),
                Random.generateInitRho(L - 1, K),
                Random.generateInitPi(L, K),
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

        this.nJK = new double[L][K];
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nJeq = new double[L];
        this.nJneq = new double[L];
        this.neqPos = new double[L];
        this.nneqPos = new double[L];
        this.coverage = new int[L];
        this.eps = new double[L];
        this.antieps = new double[L];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                this.eps[j] = eps;
                this.antieps[j] = antieps;
            }
        }
        for (Read r : reads) {
            if (r.isPaired()) {
                this.paired = true;
                break;
            }
        }
        if (this.paired) {
            this.tauOmega = new double[4][L];
        } else {
            this.tauOmega = new double[2][L];

        }
        this.init();
        this.start();
    }

    private void init() {
        int[] tau1 = new int[L];
        int[] tau2 = new int[L];
        int[] omega1 = new int[L];
        int[] omega2 = new int[L];
        for (Read r : reads) {
            for (int i = r.getWatsonBegin(); i < r.getWatsonEnd(); i++) {
                this.coverage[i - Globals.getINSTANCE().getALIGNMENT_BEGIN()]++;
            }
            tau1[r.getWatsonBegin() - Globals.getINSTANCE().getALIGNMENT_BEGIN()]++;
            omega1[r.getWatsonEnd() - Globals.getINSTANCE().getALIGNMENT_BEGIN()]++;
            if (r.isPaired()) {
                for (int i = r.getCrickBegin(); i < r.getCrickEnd(); i++) {
                    this.coverage[i - Globals.getINSTANCE().getALIGNMENT_BEGIN()]++;
                }
                tau2[r.getCrickBegin() - Globals.getINSTANCE().getALIGNMENT_BEGIN()]++;
                omega2[r.getCrickEnd() - Globals.getINSTANCE().getALIGNMENT_BEGIN()]++;
            }
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < L; i++) {
            this.tauOmega[0][i] = tau1[i] / (double) N;
            this.tauOmega[1][i] = omega1[i] / (double) N;
            if (this.paired) {
                this.tauOmega[2][i] = tau2[i] / (double) N;
                this.tauOmega[3][i] = omega2[i] / (double) N;
            }
            sb.append(this.tauOmega[0][i]);
            sb.append("\t");
            sb.append(this.tauOmega[1][i]);
            sb.append("\t");
            sb.append(this.tauOmega[2][i]);
            sb.append("\t");
            sb.append(this.tauOmega[3][i]);
            sb.append("\n");
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "twtw", sb.toString());
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
        readHMMArray = new ReadHMM[reads.length];
        this.readCount = this.reads.length;
        if (Globals.getINSTANCE().isPARALLEL_JHMM()) {
            Globals.getINSTANCE().setPARALLEL_RESTARTS_UPPER_BOUND(Integer.MAX_VALUE);
        }
        Pair<List<ReadHMM>, EInfo> invoke = Globals.getINSTANCE().getFjPool().invoke(new ReadHMMWorker(this, reads, rho, pi, mu, eps, antieps, K, L, n, this.tauOmega, 0, this.readCount));
        this.readHMMArray = invoke.getValue0().toArray(new ReadHMM[this.readCount]);
        calculate(invoke.getValue1());
    }

    public void restart() {
        this.restart++;
        this.parametersChanged = 0;
        EInfo invoke = Globals.getINSTANCE().getFjPool().invoke(new ReadHMMWorkerRecalc(K, L, n, this.readHMMArray, rho, pi, mu, eps, antieps, 0, this.readCount));
        this.calculate(invoke);
    }

    private void calculate(EInfo e) {
        this.eStep(e);
        this.calculateLoglikelihood();
        this.mStep();
    }

    private void eStep(EInfo e) {
        System.arraycopy(e.nJeq, 0, this.nJeq, 0, L);
        System.arraycopy(e.nJeq, 0, this.nJeq, 0, L);
        System.arraycopy(e.nJneq, 0, this.nJneq, 0, L);
        System.arraycopy(e.neqPos, 0, this.neqPos, 0, L);
        System.arraycopy(e.nneqPos, 0, this.nneqPos, 0, L);
//        System.arraycopy(e.coverage, 0, this.coverage, 0, L);
        for (int j = 0; j < L; j++) {
            System.arraycopy(e.nJK[j], 0, this.nJK[j], 0, K);
            for (int k = 0; k < K; k++) {
                System.arraycopy(e.nJKL[j][k], 0, this.nJKL[j][k], 0, K);
                System.arraycopy(e.nJKV[j][k], 0, this.nJKV[j][k], 0, n);
            }
        }
    }

    private void changed(double a, double b) {
        if (Math.abs(a - b) > Globals.getINSTANCE().getPCHANGE()) {
            this.parametersChanged++;
        }
    }

    private void calcMu() {
        double[] muJKV = new double[n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double AH = Globals.getINSTANCE().getALPHA_H();
                double divisor;
                double sum;
                boolean equality;
                boolean repeat;
                double[] orig = this.nJKV[j][k];
                double preSum = 0d;
                for (int v = 0; v < n; v++) {
                    preSum = orig[v];
                }
                if (preSum == 0d) {
                    for (int v = 0; v < n; v++) {
                        mu[j][k][v] = 1d / n;
                    }
                    continue;
                }
                do {
                    repeat = false;
                    equality = true;
                    sum = 0d;
                    divisor = 0d;
                    for (int v = 0; v < n; v++) {
                        muJKV[v] = this.f(this.nJKV[j][k][v] + AH);
                        sum += this.nJKV[j][k][v];
                    }
                    sum = this.f(sum + n * AH);
                    if (sum > 0) {
                        for (int v = 0; v < n; v++) {
                            muJKV[v] /= sum;
                            divisor += muJKV[v];
                        }
                        if (divisor > 0) {
                            for (int v = 0; v < n; v++) {
                                muJKV[v] /= divisor;
                            }
                        } else {
                            AH *= 2;
                        }
                    } else {
                        AH *= 2;
                    }
                    if (divisor > 0) {
                        for (int v = 1; v < n; v++) {
                            if (muJKV[v - 1] != muJKV[v]) {
                                equality = false;
                                break;
                            }
                        }
                        if (equality) {
                            double min = Double.MAX_VALUE;
                            for (int v = 1; v < n; v++) {
                                if (this.nJKV[j][k][v] > 0) {
                                    min = Math.min(min, this.nJKV[j][k][v]);
                                }
                            }
                            for (int v = 0; v < n; v++) {
                                this.nJKV[j][k][v] /= min;
                            }
                            AH = Globals.getINSTANCE().getALPHA_H();
                            repeat = true;
                        }
                    }
                } while (sum == 0 || divisor == 0 || repeat);

                for (int v = 0; v < n; v++) {
                    this.changed(mu[j][k][v], muJKV[v]);
                    mu[j][k][v] = muJKV[v];
                }
            }
        }
    }

    private void calcRho() {
        double[] rhoJKL = new double[K];
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double AZ = Globals.getINSTANCE().getALPHA_Z();
                double divisor;
                double sum;
                do {
                    sum = 0d;
                    divisor = 0d;
                    for (int l = 0; l < K; l++) {
                        rhoJKL[l] = this.f(this.nJKL[j][k][l] + AZ);
                        sum += this.nJKL[j][k][l];
                    }
                    sum = this.f(sum + K * AZ);
                    if (sum > 0) {
                        for (int l = 0; l < K; l++) {
                            rhoJKL[l] /= sum;
                            divisor += rhoJKL[l];
                        }
                        if (divisor > 0) {
                            for (int l = 0; l < K; l++) {
                                rhoJKL[l] /= divisor;
                            }
                        } else {
                            AZ *= 10;
                        }
                    } else {
                        AZ *= 10;
                    }
                } while (sum == 0 || divisor == 0);

                boolean equality = true;
                for (int l = 1; l < K; l++) {
                    if (rhoJKL[l - 1] != rhoJKL[l]) {
                        equality = false;
                    }
                }
                if (equality) {
                    double max = -1;
                    int maxV = -1;
                    for (int l = 0; l < K; l++) {
                        if (this.nJKL[j][k][l] > max) {
                            max = this.nJKL[j][k][l];
                            maxV = l;
                        }
                    }
                    for (int l = 0; l < K; l++) {
                        if (l == maxV) {
                            rhoJKL[l] = 1;
                        } else {
                            rhoJKL[l] = 0;
                        }
                    }
                }

                for (int l = 0; l < K; l++) {
                    this.changed(rho[j - 1][k][l], rhoJKL[l]);
                    rho[j - 1][k][l] = rhoJKL[l];
                }
            }
        }
    }

    private double f(double upsilon) {
        if (upsilon == 0d) {
            return 0d;
        }
        return Math.exp(rho_phi(upsilon));
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

    private void calcPi() {
        double sumK = 0d;
        for (int k = 0; k < K; k++) {
            pi[k] = this.nJK[0][k];
            sumK += pi[k];
        }
        for (int k = 0; k < K; k++) {
            pi[k] /= sumK;
        }
    }

    private void mStep() {
        if (!Globals.getINSTANCE().isNO_RECOMB()) {
            this.calcRho();
        }
        this.calcPi();
        this.calcMu();

        if (!Globals.getINSTANCE().isFLAT_EPSILON_PRIOR()) {
            double ew = 0d;
            int nonzero = 0;
            for (int j = 0; j < L; j++) {
                if (this.nneqPos[j] != 0d) {
                    ew += this.nneqPos[j] / N;
                    nonzero++;
                }
            }
            double a = 0d, b = 0d;
            if (ew != 0d) {
                ew /= nonzero;
                a = 20;
                b = (-a * ew + a + 2 * ew - 1) / ew;
            }
//        if (!Globals.FIX_EPSILON) {
//            for (int j = 0; j < L; j++) {
//                this.eps[j] = this.nneqPos[j] / (N * (n - 1));
//            }
//        } else {
            for (int j = 0; j < L; j++) {
                if (ew == 0d) {
                    this.eps[j] = 0;
                    this.antieps[j] = 1;
                } else {
                    this.eps[j] = f(ew + a) / f((coverage[j] * (n - 1)) + a + b);
                    this.antieps[j] = 1 - (n - 1) * eps[j];
                }
            }
//    }
        }
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

    public double[] getEps() {
        return eps;
    }

    public double getLoglikelihood() {
        return loglikelihood;
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

    public int getRestart() {
        return restart;
    }

    public int getParametersChanged() {
        return parametersChanged;
    }
}