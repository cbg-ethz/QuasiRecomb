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
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorker;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.ReadHMMWorkerRecalc;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import com.google.common.util.concurrent.AtomicDouble;
import com.google.common.util.concurrent.AtomicDoubleArray;
import java.util.List;

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
    private AtomicDoubleArray[][] nJKL;
    private AtomicDoubleArray[][] nJKV;
    private AtomicDoubleArray neqPos;
    private AtomicDoubleArray nneqPos;
    private AtomicDoubleArray nJeq;
    private AtomicDoubleArray nJneq;
    private AtomicDouble loglikelihood;
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

        this.nJKL = new AtomicDoubleArray[L][K];
        this.nJKV = new AtomicDoubleArray[L][K];
        this.nJeq = new AtomicDoubleArray(L);
        this.nJneq = new AtomicDoubleArray(L);
        this.neqPos = new AtomicDoubleArray(L);
        this.nneqPos = new AtomicDoubleArray(L);
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                this.nJKL[j][k] = new AtomicDoubleArray(K);
                this.nJKV[j][k] = new AtomicDoubleArray(n);
            }
        }
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
        this.loglikelihood = new AtomicDouble(0d);
        for (ReadHMM r : this.readHMMArray) {
            for (int j = 0; j < L; j++) {
                double a = r.getC(j);
                double x = r.getC(j);
                double y = r.getC(j);
                double b = Math.log(a);
                double bc = Math.log(r.getC(j));
                double c = b * r.getCount();
                double cc = Math.log(r.getC(j)) * r.getCount();
                if (Double.isInfinite(c)) {
                    System.err.println(r.getC(j));
                    System.err.println(Math.log(r.getC(j)));
                    double aa = Math.log(r.getC(j));
                    System.err.println(aa);
                    Object xx = Math.log(r.getC(j));
                    System.err.println(xx);
                }
                double update = (Math.log(r.getC(j)) * r.getCount());
                if (Double.isInfinite(update)) {
                    System.out.println(update);
                    System.out.println("LLH oo");
                }
                this.loglikelihood.addAndGet(Math.log(r.getC(j)) * r.getCount());

                if (Double.isInfinite(loglikelihood.doubleValue())) {
                    System.out.println("LLH oo");
                }
                if (Double.isNaN(loglikelihood.doubleValue())) {
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
        List<ReadHMM> invoke = Globals.getINSTANCE().getFjPool().invoke(new ReadHMMWorker(this, reads, 0, this.readCount));
        this.readHMMArray = invoke.toArray(new ReadHMM[this.readCount]);
        calculate();
    }

    public void restart() {
        this.restart++;
        this.parametersChanged = 0;
//        Pair<List<ReadHMM>, EInfo> invoke = Globals.getINSTANCE().getFjPool().invoke(new ReadHMMWorker(this, reads, 0, this.readCount));
//        this.readHMMArray = invoke.getValue0().toArray(new ReadHMM[this.readCount]);
//        this.calculate(invoke.getValue1());
        Globals.getINSTANCE().getFjPool().invoke(new ReadHMMWorkerRecalc(this, this.readHMMArray, 0, this.readCount));
        this.calculate();
    }

    private void calculate() {
        this.calculateLoglikelihood();
        this.mStep();

        this.nJKL = new AtomicDoubleArray[L][K];
        this.nJKV = new AtomicDoubleArray[L][K];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                this.nJKL[j][k] = new AtomicDoubleArray(K);
                this.nJKV[j][k] = new AtomicDoubleArray(n);
            }
        }
        this.nJeq = new AtomicDoubleArray(L);
        this.nJneq = new AtomicDoubleArray(L);
        this.neqPos = new AtomicDoubleArray(L);
        this.nneqPos = new AtomicDoubleArray(L);
//        for (int j = 0; j < L; j++) {
//            this.nneqPos[j] = 0d;
//            this.neqPos[j] = 0d;
//        }
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

//                double AH = Globals.getINSTANCE().getALPHA_H();
                double AH = 1e-4;
                double divisor;
                double sum = 0d;

                double max = Double.MIN_VALUE;
                for (int v = 0; v < n; v++) {
                    max = Math.max(max, this.nJKV[j][k].get(v));
                }
                if (max < 1) {
                    for (int v = 0; v < n; v++) {
                        sum += this.nJKV[j][k].get(v);
                    }
                    if (sum > 0) {
                        for (int v = 0; v < n; v++) {
                            muJKV[v] = this.nJKV[j][k].get(v) / sum;
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
                                System.out.println(this.nJKV[j][k].get(v) + "\t" + muJKV[i]);
                            }
                            System.out.println("sum:\t" + sum);
                        }
                    }
                } else {
                    do {
                        sum = 0d;
                        divisor = 0d;
                        for (int v = 0; v < n; v++) {
                            muJKV[v] = this.f(this.nJKV[j][k].get(v) + AH);
                            sum += this.nJKV[j][k].get(v);
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
                    max = Math.max(max, this.nJKL[j][k].get(l));
                }
                if (max < 1e-10) {
                    for (int l = 0; l < K; l++) {
                        sum += this.nJKL[j][k].get(l);
                    }
                    if (sum > 0) {
                        for (int l = 0; l < K; l++) {
                            rhoJKL[l] = this.nJKL[j][k].get(l) / sum;
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
                                System.out.println(this.nJKV[j - 1][k].get(l) + "\t" + rhoJKL[i]);
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
                            rhoJKL[l] = this.f(this.nJKL[j][k].get(l) + AZ);
                            sum += this.nJKL[j][k].get(l);
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
                for (int v = 0; v < n; v++) {
                    pi[k] += this.nJKV[j][k].get(v);
                    sumK += this.nJKV[j][k].get(v);
                }
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
//            int nonzero = 0;
            for (int j = 0; j < L; j++) {
                if (this.nneqPos.get(j) != 0d) {
                    ew += this.nneqPos.get(j) / N;
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
                    this.eps[j] = f(this.nneqPos.get(j) + a) / f((coverage[j] * (n - 1)) + a + b);
                    if (this.eps[j] > 1d / n) {
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

    public int getn() {
        return n;
    }

    public double[] getEps() {
        return eps;
    }

    public double[] getAntieps() {
        return antieps;
    }

    public double getLoglikelihood() {
        return loglikelihood.get();
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
                if (Math.abs(max - sum) < 1e-8) {
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
                    max = Math.max(this.nJKV[j][k].get(v), max);
                    sum += this.nJKV[j][k].get(v);
                }
                if (Math.abs(max - sum) < 1e-8) {
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
                if (Math.abs(max - sum) < 1e-8) {
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
                    max = Math.max(this.nJKL[j][k].get(l), max);
                    sum += this.nJKL[j][k].get(l);
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

    public void addnJKV(int j, int k, int v, double value) {
//        synchronized (this.nJKV[j][k]) {
        this.nJKV[j][k].addAndGet(v, value);
//        }
    }

    public void addnJKL(int j, int k, int l, double value) {
//        synchronized (this.nJKL[j][k]) {
        this.nJKL[j][k].addAndGet(l, value);
//        }
    }

    public void addnJeq(int l, double value) {
//        synchronized (this.nJeq) {
        this.nJeq.addAndGet(l, value);
//        }
    }

    public void addnJneq(int l, double value) {
//        synchronized (this.nJneq) {
        this.nJneq.addAndGet(l, value);
//        }
    }

    public void addneqPos(int l, double value) {
//        synchronized (this.neqPos[l]) {
        this.neqPos.getAndAdd(l, value);
//        }
    }

    public void addnneqPos(int l, double value) {
//        synchronized (this.nneqPos[l]) {
        this.nneqPos.addAndGet(l, value);
//        }
    }
}
