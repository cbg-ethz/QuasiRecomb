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

import cc.mallet.types.Dirichlet;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.EInfo;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorker;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorkerRecalc;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;
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
            this.tauOmega = new double[4][L + 1];
        } else {
            this.tauOmega = new double[2][L + 1];

        }
        this.init();
        this.start();
    }

    private void init() {
        int[] tau1 = new int[L + 1];
        int[] tau2 = new int[L + 1];
        int[] omega1 = new int[L + 1];
        int[] omega2 = new int[L + 1];
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
        for (int i = 0; i < L + 1; i++) {
            this.tauOmega[0][i] = tau1[i] / (double) N;
            this.tauOmega[1][i] = omega1[i] / (double) N;
            sb.append(this.tauOmega[0][i]);
            sb.append("\t");
            sb.append(this.tauOmega[1][i]);
            sb.append("\t");
            if (this.paired) {
                this.tauOmega[2][i] = tau2[i] / (double) N;
                this.tauOmega[3][i] = omega2[i] / (double) N;
                sb.append(this.tauOmega[2][i]);
                sb.append("\t");
                sb.append(this.tauOmega[3][i]);
            }
            sb.append("\n");
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "twtw", sb.toString());
    }

    private void calculateLoglikelihood() {
        this.loglikelihood = 0d;
        for (ReadHMM r : this.readHMMArray) {
            for (int j = 0; j < L; j++) {
                this.loglikelihood += Math.log(r.getC(j)) * r.getCount();
                if (Double.isNaN(loglikelihood)) {
                    System.out.println("LLH NaN");
                }
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
        double[] muJKV;
        for (int j = 0; j < L; j++) {
            muJKV = new double[n];
            for (int k = 0; k < K; k++) {

//                double AH = Globals.getINSTANCE().getALPHA_H();
                double AH = 1e-4;
                double divisor;
                double sum = 0d;

                double max = Double.MIN_VALUE;
                for (int v = 0; v < n; v++) {
                    max = Math.max(max, this.nJKV[j][k][v]);
                }
                if (max < 1) {
                    for (int v = 0; v < n; v++) {
                        sum += this.nJKV[j][k][v];
                    }
                    if (sum > 0) {
                        for (int v = 0; v < n; v++) {
                            muJKV[v] = this.nJKV[j][k][v] / sum;
                        }
                    } else {
//                        muJKV = new Dirichlet(n).nextDistribution();
                        for (int v = 0; v < n; v++) {
                            muJKV[v] = 1d / n;
                        }
                    }
                    for (int v = 0; v < n; v++) {
                        this.changed(mu[j][k][v], muJKV[v]);
                        mu[j][k][v] = muJKV[v];
                        if (Double.isNaN(muJKV[v])) {
                            System.out.println("R MU plain nan, j " + j + ", k " + k);
                            for (int i = 0; i < n; i++) {
                                System.out.println(this.nJKV[j][k][v] + "\t" + muJKV[i]);
                            }
                            System.out.println("sum:\t" + sum);
                        }
                    }
                } else {
                    do {
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
                                AH *= 10;
                            }
                        } else {
                            AH *= 10;
                        }
                        if (Double.isNaN(sum) || Double.isNaN(divisor)) {
                            System.out.println("a");
                        }
                    } while (sum == 0 || divisor == 0);
                    max = Double.MIN_VALUE;
                    for (int v = 0; v < n; v++) {
                        max = Math.max(max, muJKV[v]);
                    }
                    if (max == 1d) {
                        for (int v = 0; v < n; v++) {
                            if (muJKV[v] < 1) {
                                muJKV[v] = 0d;
                            }
                        }
                    }
                    for (int v = 0; v < n; v++) {
                        this.changed(mu[j][k][v], muJKV[v]);
                        mu[j][k][v] = muJKV[v];
                        if (Double.isNaN(muJKV[v])) {
                            System.out.println("R MU reg nan");
                            System.exit(0);
                        }
                    }
                }
            }
        }
    }

    private void calcRho() {
        double[] rhoJKL = new double[K];
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double AZ = 0.001;
                double divisor;
                double sum = 0d;

                double max = Double.MIN_VALUE;
                for (int l = 0; l < K; l++) {
                    max = Math.max(max, this.nJKL[j][k][l]);
                }
                if (max < 1e-10) {
                    for (int l = 0; l < K; l++) {
                        sum += this.nJKL[j][k][l];
                    }
                    if (sum > 0) {
                        for (int l = 0; l < K; l++) {
                            rhoJKL[l] = this.nJKL[j][k][l] / sum;
                        }
//                    } else {
//                        for (int l = 0; l < K; l++) {
//                            if (k == l) {
//                                rhoJKL[l] = 1d;
//                            } else {
//                                rhoJKL[l] = 0d;
//                            }
//                        }
                    }
                    for (int l = 0; l < K; l++) {
                        this.changed(rho[j - 1][k][l], rhoJKL[l]);
                        rho[j - 1][k][l] = rhoJKL[l];
                        if (Double.isNaN(rhoJKL[l])) {
                            System.out.println("R RHO plain nan, j " + j + ", k " + k);
                            for (int i = 0; i < n; i++) {
                                System.out.println(this.nJKV[j - 1][k][l] + "\t" + rhoJKL[i]);
                            }
                            System.out.println("sum:\t" + sum);
                            System.exit(0);
                        }
                    }
                } else {
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

                    for (int l = 0; l < K; l++) {
                        this.changed(rho[j - 1][k][l], rhoJKL[l]);
                        rho[j - 1][k][l] = rhoJKL[l];
                    }
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
        return Dirichlet.digamma(upsilon);
    }

    private void calcPi() {
        double sumK = 0d;
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                pi[k] += this.nJK[j][k];
                sumK += this.nJK[j][k];
            }
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
                }
            }
            double a = 0d, b = 0d;
            if (ew != 0d) {
                ew /= L;
                a = 20;
                b = (-a * ew + a + 2 * ew - 1) / ew;
            }
            for (int j = 0; j < L; j++) {
                if (ew == 0d) {
                    this.eps[j] = 0;
                    this.antieps[j] = 1;
                } else {
                    this.eps[j] = f(this.nneqPos[j] + a) / f((coverage[j] * (n - 1)) + a + b);
                    if (this.eps[j] > 1d/n) {
                        this.eps[j] = 0.05;
                    }
                    this.antieps[j] = 1 - (n - 1) * eps[j];
                }
            }
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

    public int getMuFlats() {
        int flats = 0;
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double max = 0;
                double sum = 0;
                for (int v = 0; v < n; v++) {
                    max = Math.max(this.mu[j][k][v], max);
                    sum += this.mu[j][k][v];
//                    max = Math.max(this.nJKV[j][k][v], max);
//                    sum += this.nJKV[j][k][v];
                }
                if (max != sum) {
                    flats++;
                }
            }
        }
        return flats;
    }

    public int getNjkvFlats() {
        int flats = 0;
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double max = 0;
                double sum = 0;
                for (int v = 0; v < n; v++) {
                    max = Math.max(this.nJKV[j][k][v], max);
                    sum += this.nJKV[j][k][v];
                }
                if (max != sum) {
                    flats++;
                }
            }
        }
        return flats;
    }

    public int getRhoFlats() {
        int flats = 0;
        for (int j = 0; j < L - 1; j++) {
            for (int k = 0; k < K; k++) {
                double max = 0;
                double sum = 0;
                for (int l = 0; l < K; l++) {
                    max = Math.max(this.rho[j][k][l], max);
                    sum += this.rho[j][k][l];
                }
                if (max != sum) {
                    flats++;
                }
            }
        }
        return flats;
    }

    public int getNjklFlats() {
        int flats = 0;
        for (int j = 0; j < L - 1; j++) {
            for (int k = 0; k < K; k++) {
                double max = 0;
                double sum = 0;
                for (int l = 0; l < K; l++) {
                    max = Math.max(this.nJKL[j][k][l], max);
                    sum += this.nJKL[j][k][l];
                }
                if (max != sum) {
                    flats++;
                }
            }
        }
        return flats;
    }

    public double[][] getTauOmega() {
        return tauOmega;
    }
}
