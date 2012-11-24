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
import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.CallableReadHMM;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.FutureTask;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Offers a HMM with the E and M step.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class JHMM {

    private int N;
    private int L;
    private int K;
    private int n;
    //rho[j][k][l] := transition prob. at position j, for l given k
    private double[][] tauOmega;
    private double[][][] rho;
    private double[] pi;
    private double[][][] mu;
    private double[] eps;
    private double[] antieps;
    private double[][][] nJKL;
    private double[][][] nJKV;
    private double[] nneqPos;
    private double loglikelihood;
    private double[] muPrior;
    private Read[] reads;
    private int restart = 0;
    private int coverage[];
    private int parametersChanged = 0;
    private boolean paired;
    private Map<Integer, TempJHMMStorage> garage = new ConcurrentHashMap<>();
    private final List<Integer> available = new ArrayList<>();

    public JHMM(Read[] reads, int N, int L, int K, int n, double epsilon) {
        this(reads, N, L, K, n, epsilon,
                Random.generateInitRho(L - 1, K),
                Random.generateInitPi(L, K),
                Random.generateMuInit(L, K, n));
    }

    public JHMM(Read[] reads, int N, int L, int K, int n, double eps, double[][][] rho, double[] pi, double[][][] mu) {
        this.eps = new double[L];
        this.antieps = new double[L];
        for (int j = 0; j < L; j++) {
            this.eps[j] = eps;
            this.antieps[j] = (1 - (n - 1) * eps);
        }
        this.prepare(reads, N, L, K, n, rho, pi, mu);
    }

    public JHMM(Read[] reads, int N, int L, int K, int n, double[] eps, double[][][] rho, double[] pi, double[][][] mu) {
        this.eps = eps;
        this.antieps = new double[L];
        for (int j = 0; j < L; j++) {
            this.antieps[j] = (1 - (n - 1) * eps[j]);
        }
        this.prepare(reads, N, L, K, n, rho, pi, mu);
    }

    private void prepare(Read[] reads, int N, int L, int K, int n, double[][][] rho, double[] pi, double[][][] mu) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.reads = reads;
        this.rho = rho;
        this.mu = mu;
        this.pi = pi;

        this.muPrior = new double[n];
        for (int i = 0; i < n; i++) {
            this.muPrior[i] = 0.001;
        }

        this.coverage = new int[L];
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
                this.coverage[i - Globals.getINSTANCE().getALIGNMENT_BEGIN()] += r.getCount();
            }
            tau1[r.getWatsonBegin() - Globals.getINSTANCE().getALIGNMENT_BEGIN()] += r.getCount();
            omega1[r.getWatsonEnd() - Globals.getINSTANCE().getALIGNMENT_BEGIN()] += r.getCount();
            if (r.isPaired()) {
                for (int i = r.getCrickBegin(); i < r.getCrickEnd(); i++) {
                    this.coverage[i - Globals.getINSTANCE().getALIGNMENT_BEGIN()] += r.getCount();
                }
                tau2[r.getCrickBegin() - Globals.getINSTANCE().getALIGNMENT_BEGIN()] += r.getCount();
                omega2[r.getCrickEnd() - Globals.getINSTANCE().getALIGNMENT_BEGIN()] += r.getCount();
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
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "twtw", sb.toString());
    }

    private void start() {
        calculate();
    }

    public void restart() {
        this.restart++;
        this.parametersChanged = 0;
        this.calculate();
    }

    private void compute() {
        if (Globals.getINSTANCE().isSTORAGE()) {
            available.clear();
            garage.clear();
            for (int i = 0; i < Globals.getINSTANCE().getCpus(); i++) {
                available.add(i);
                garage.put(i, new TempJHMMStorage(L, K, n, i));
            }
        }
        this.loglikelihood = 0d;
        List<FutureTask<Double>> taskList = new ArrayList<>();

        for (int i = 0; i < reads.length; i++) {
            FutureTask<Double> futureTask_1 = new FutureTask<>(new CallableReadHMM(this, reads[i]));
            taskList.add(futureTask_1);
            Globals.getINSTANCE().getExecutor().execute(futureTask_1);
        }

        for (int j = 0; j < taskList.size(); j++) {
            FutureTask<Double> futureTask = taskList.get(j);
            try {
                loglikelihood += futureTask.get();
                Globals.getINSTANCE().printPercentage(K, (double) j / reads.length);
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(JHMM.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        if (Globals.getINSTANCE().isSTORAGE()) {
            if (available.size() != Globals.getINSTANCE().getCpus()) {
                throw new IllegalStateException("Not all storages have been returned");
            }
            Iterator<TempJHMMStorage> iterator = this.garage.values().iterator();
            TempJHMMStorage store = iterator.next();
            while (iterator.hasNext()) {
                store.merge(iterator.next());
            }
            for (int j = 0; j < L; j++) {
                for (int k = 0; k < K; k++) {
                    for (int l = 0; l < K; l++) {
                        this.nJKL[j][k][l] += store.getnJKL()[j][k][l];
                    }
                    for (int v = 0; v < n; v++) {
                        this.nJKV[j][k][v] += store.getnJKV()[j][k][v];
                    }
                }
                this.nneqPos[j] += store.getNneqPos()[j];
            }
        }
    }

    public TempJHMMStorage getStorage() {
        synchronized (this.available) {
            if (available.iterator().hasNext()) {
                Integer i = available.iterator().next();
                available.remove(i);
                return garage.get(i);
            } else {
                throw new IllegalStateException("No storages left");
            }
        }
    }

    public void free(int id) {
        synchronized (this.available) {
            this.available.add(id);
        }
    }

    private void calculate() {
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nneqPos = new double[L];
        Globals.getINSTANCE().resetTimer();
        this.compute();
        this.mStep();
    }

    private void changed(double a, double b) {
        if (Math.abs(a - b) > Globals.getINSTANCE().getPCHANGE()) {
            this.parametersChanged++;
        }
    }

    private void calcMu() {
        double[] muJKV;
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                muJKV = Regularizations.regularizeOnce(this.nJKV[j][k], restart, muPrior, 100);
                for (int v = 0; v < n; v++) {
                    this.changed(mu[j][k][v], muJKV[v]);
                    mu[j][k][v] = muJKV[v];
                    if (Double.isNaN(muJKV[v])) {
                        System.out.println("R nan, j " + j + ", k " + k);
                        for (int i = 0; i < n; i++) {
                            System.out.println(this.nJKV[j][k][v] + "\t" + muJKV[i]);
                        }
                    }
                }
            }
        }
    }

    private void calcRho() {
        double[] rhoJKL;
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {

                double max = 0d;
                int lPrime = -1;
                for (int l = 0; l < K; l++) {
                    if (this.nJKL[j][k][l] > max) {
                        max = this.nJKL[j][k][l];
                        lPrime = l;
                    }
                }

                double[] rhoPrior = new double[K];
                for (int l = 0; l < K; l++) {
                    if (k == l) {
                        if (Math.abs(max - 1d) < 1e-8 && lPrime != k) {
                            rhoPrior[l] = 1;
                        } else {
                            rhoPrior[l] = 0.001;
                        }
                    } else {
                        rhoPrior[l] = 0.0001;
                    }
                }
                rhoJKL = Regularizations.regularizeOnce(this.nJKL[j][k], restart, rhoPrior, 100);

                for (int l = 0; l < K; l++) {
                    this.changed(rho[j - 1][k][l], rhoJKL[l]);
                    rho[j - 1][k][l] = rhoJKL[l];
                }
            }
        }
    }

    private void calcPi() {
        double sumK = 0d;
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
                    pi[k] += this.nJKV[j][k][v];
                    sumK += this.nJKV[j][k][v];
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
            if (Globals.getINSTANCE().isUNINFORMATIVE_EPSILON_PRIOR()) {
                for (int j = 0; j < L; j++) {
                    if (this.nneqPos[j] != 0d) {
                        ew += this.nneqPos[j] / N;
                    }
                }
                if (ew != 0d) {
                    ew /= L;
                }
            } else {
                ew = .008;
            }
            double a = 20;
            double b = (-a * ew + a + 2 * ew - 1) / ew;
            for (int j = 0; j < L; j++) {
                if (ew == 0d) {
                    this.eps[j] = 0;
                    this.antieps[j] = 1;
                } else {
                    this.eps[j] = Regularizations.f(this.nneqPos[j] + a) / Regularizations.f((coverage[j] * (n - 1)) + a + b);
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
                if (max < sum) {
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
                if (max < sum) {
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
                if (max < sum) {
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
                if (max < sum) {
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
