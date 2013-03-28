/**
 * Copyright (c) 2011-2013 Armin Töpfer
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

import ch.ethz.bsse.quasirecomb.distance.KullbackLeibler;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.CallableReadHMMList;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.javatuples.Pair;
import org.javatuples.Triplet;

/**
 * Offers a HMM with the E and M step.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class JHMM extends Garage {

    protected int Kmin;
    protected int N;
    protected int L;
    protected int K;
    protected int n;
    protected double[][] snv;
    //rho[j][k][l] := transition prob. at position j, for l given k
    protected double[][][] rho;
    protected double[] pi;
    protected double[][][] mu;
    protected double[] eps;
    protected double[] antieps;
    protected double loglikelihood;
    protected double[][][] nJKL;
    protected double[][][] nJKV;
    protected double[] nneqPos;
    protected double[] muPrior;
    protected Read[] allReads;
    protected int restart = 0;
    protected int coverage[];
    protected int muChanged = 0;
    protected int rhoChanged = 0;
    protected boolean paired;
    protected double[][][] rho_old;
    protected double[][][] mu_old;
    private int oldFlatMu = -1;
    private boolean biasMu = false;
    private int biasCounter = 0;
    private int unBiasCounter = 0;
    private int s = 0;

    public JHMM(Read[] reads, int N, int L, int K, int n, double epsilon, int Kmin) {
        this(reads, N, L, K, n, epsilon,
                Random.generateInitRho(L - 1, K),
                Random.generateInitPi(L, K),
                Random.generateMuInit(L, K, n), Kmin);
    }

    public JHMM(Read[] reads, int N, int L, int K, int n, double eps, double[][][] rho, double[] pi, double[][][] mu, int Kmin) {
        this.eps = new double[L];
        this.antieps = new double[L];
        for (int j = 0; j < L; j++) {
            this.eps[j] = eps;
            this.antieps[j] = (1 - (n - 1) * eps);
        }
        this.Kmin = Kmin;
        this.prepare(reads, N, L, K, n, rho, pi, mu);
        this.compute();
    }

    public JHMM(Read[] reads, int N, int L, int K, int n, double[] eps, double[][][] rho, double[] pi, double[][][] mu, int Kmin) {
        this.eps = eps;
        this.antieps = new double[L];
        for (int j = 0; j < L; j++) {
            this.antieps[j] = (1 - (n - 1) * eps[j]);
        }
        this.Kmin = Kmin;
        this.prepare(reads, N, L, K, n, rho, pi, mu);
        this.compute();
    }

    private void compute() {
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nneqPos = new double[L];
        this.eStep();
        this.mStep();
        s++;
    }
    List<Callable<Double>> callables = new LinkedList<>();

    private void eStep() {
        clearGarage(L, K, n);
        this.loglikelihood = 0d;
        if (callables.isEmpty()) {
            final int readAmount = allReads.length;
            int max;
            if (Globals.getINSTANCE().getSTEPS() == 2) {
                max = (int) (readAmount / (Runtime.getRuntime().availableProcessors() - 1));
            } else {
                max = Globals.getINSTANCE().getSTEPS();
            }
            for (int i = 0; i < readAmount; i += max) {
                int b = i + max;
                if (b >= readAmount) {
                    b = readAmount;
                }
                callables.add(new CallableReadHMMList(this, Arrays.copyOfRange(allReads, i, b)));
                Globals.getINSTANCE().printPercentage(K, (double) i / readAmount, Kmin);
            }
        }
        List<Future<Double>> results = null;
        try {
            results = Globals.getINSTANCE().getExecutor().invokeAll(callables);
        } catch (InterruptedException ex) {
            Logger.getLogger(JHMM.class.getName()).log(Level.SEVERE, null, ex);
        }

        for (int i = 0; i < results.size(); i++) {
            try {
                Double llh = results.get(i).get();
                loglikelihood += llh;
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(JHMM.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        updateExpectedCounts();
    }

    private void updateExpectedCounts() {
        if (Globals.getINSTANCE().isSTORAGE()) {
            TempJHMMStorage store = mergeGarage();
            this.nJKL = new double[L][K][K];
            this.nJKV = new double[L][K][n];
            this.nneqPos = new double[L];
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

    public void computeSNVPosterior() {
        double[][][] nJKV_norm = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double min = Double.MAX_VALUE;
                for (int v = 0; v < n; v++) {
                    min = Math.min(this.nJKV[j][k][v], min);
                }
                for (int v = 0; v < n; v++) {
                    nJKV_norm[j][k][v] = this.nJKV[j][k][v] - min;
                }
            }
        }
        for (int j = 0; j < L; j++) {
            for (int v = 0; v < n; v++) {
                this.snv[j][v] = 0;
                for (int k = 0; k < K; k++) {
                    this.snv[j][v] += nJKV_norm[j][k][v] / this.coverage[j];
                }
            }
        }
    }

    private void maximizeMu() {
        double[] muJKV;
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double eta = 0;
                if (Globals.getINSTANCE().getINTERPOLATE_MU() > 0) {
                    eta = Math.pow(Math.pow(s, 2) + 2, -Globals.getINSTANCE().getINTERPOLATE_MU());
                    muJKV = Regularizations.step(this.nJKV[j][k], this.mu[j][k], eta, true);
                } else {
                    muJKV = Regularizations.ml(this.nJKV[j][k]);
                }
                double mult = Globals.getINSTANCE().getMULT_MU();
                double[] muPriorLocal = new double[n];
                System.arraycopy(this.muPrior, 0, muPriorLocal, 0, n);
                boolean repeat = true;
                do {
                    muJKV = Regularizations.regularizeOnce(muJKV, restart, muPriorLocal, mult);
                    double prev = muJKV[0];
                    for (int v = 1; v < n; v++) {
                        if (prev != muJKV[v] && prev > 0d) {
                            repeat = false;
                        }
                    }
                    if (repeat) {
                        mult *= 2;
                        for (int v = 0; v < n; v++) {
                            muPriorLocal[v] *= 10;
                        }
                    }
                    if (mult > 1000) {
                        repeat = false;
                    }
                } while (repeat);
                if (Globals.getINSTANCE().isMAX()) {
                    double max = 0;
                    int maxV = 0;
                    for (int v = 0; v < n; v++) {
                        if (muJKV[v] > max) {
                            max = muJKV[v];
                            maxV = v;
                        }
                    }
                    for (int v = 0; v < n; v++) {
                        muJKV[v] = (v == maxV) ? 1 : 0;
                    }
                }
                for (int v = 0; v < n; v++) {
                    this.changedMu(mu[j][k][v], muJKV[v]);
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

    private boolean maximizeRho() {
        boolean forceRho = false;
        double[] rhoJKL = null;
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {
                double eta = 0;
                if (Globals.getINSTANCE().getINTERPOLATE_RHO() > 0) {
                    eta = Math.pow(Math.pow(s, 2) + 2, -Globals.getINSTANCE().getINTERPOLATE_RHO());
                    rhoJKL = Regularizations.step(this.nJKL[j][k], this.rho[j - 1][k], eta, false);
                } else {
                    rhoJKL = Regularizations.ml(this.nJKL[j][k]);
                }
                double max = 0d;
                int lPrime = -1;
                double sum = 0d;
                double[] intermediate = new double[K];
                for (int l = 0; l < K; l++) {
                    intermediate[l] = this.nJKL[j][k][l];
                    sum += intermediate[l];
                }
                for (int l = 0; l < K; l++) {
                    intermediate[l] /= sum;
                    if (intermediate[l] > max) {
                        max = intermediate[l];
                        lPrime = l;
                    }
                }

                double mult = Globals.getINSTANCE().getMULT_RHO();
                double[] rhoPrior = new double[K];
                boolean fix = false;
                for (int l = 0; l < K; l++) {
                    if (k == l) {
                        if (max > 0.5 && lPrime != k && Globals.getINSTANCE().isSPIKERHO()) {
                            rhoPrior[l] = 100;
                            fix = true;
                            forceRho = true;
                            break;
                        } else {
                            rhoPrior[l] = Globals.getINSTANCE().getALPHA_Z() * 10;
                        }
                    } else {
                        rhoPrior[l] = Globals.getINSTANCE().getALPHA_Z();
                    }
                }
                if (!fix) {
                    rhoJKL = Regularizations.regularizeOnceRho(k, rhoJKL, restart, rhoPrior, mult);
                    lPrime = -1;
                    max = -1;
                    for (int l = 0; l < K; l++) {
                        if (rhoJKL[l] > max) {
                            max = rhoJKL[l];
                            lPrime = l;
                        }
                    }
                    if (max > 0.5 && lPrime != k && Globals.getINSTANCE().isSPIKERHO()) {
                        fix = true;
                        forceRho = true;
                    }
                }
                if (fix) {
                    rhoJKL = new double[K];
                    for (int l = 0; l < K; l++) {
                        rhoJKL[l] = l == k ? 1 : 0;
                    }
                }
                for (int l = 0; l < K; l++) {
                    this.changedRho(rho[j - 1][k][l], rhoJKL[l]);
                    rho[j - 1][k][l] = rhoJKL[l];
                }
            }
        }
        return forceRho;
    }

    private void maximizePi() {
        double[] piTmp = new double[K];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
                    piTmp[k] += this.nJKV[j][k][v];
                }
            }
        }

        double eta = Math.pow(s + 2, -1);
        pi = Regularizations.step(piTmp, pi, eta, false);
    }

    private void mStep() {
        this.computeSNVPosterior();
        boolean forceRho = false;
        if (!Globals.getINSTANCE().isNO_RECOMB()) {
            forceRho = this.maximizeRho();
        }
        this.maximizePi();
        this.maximizeMu();
        if (forceRho) {
            for (int j = 0; j < L; j++) {
                for (int k = 0; k < K; k++) {
                    double sum = 0;
                    for (int v = 0; v < n; v++) {
                        this.mu[j][k][v] += Math.random() / 100d;
                        sum += this.mu[j][k][v];
                    }
                    for (int v = 0; v < n; v++) {
                        this.mu[j][k][v] /= sum;
                    }
                }
            }
        }
        if (Globals.getINSTANCE().isBIAS_MU()) {
            int currentFlatMu = getMuFlats();
            if (biasCounter < 1) {
                this.biasMu = true;
            } else {
                this.unBiasCounter++;
                this.biasMu = false;
                if (unBiasCounter > 200) {
                    this.biasCounter = 0;
                    this.unBiasCounter = 0;
                }
            }
            if (oldFlatMu == currentFlatMu) {
                if (biasMu) {
                    this.biasMu();
                    System.out.print("BIAS\t");
                    biasCounter++;
                } else {
                    System.out.print("UNB\t");
                }
            } else {
                System.out.print("DIFF\t");
            }
            oldFlatMu = currentFlatMu;
        }
        if (!Globals.getINSTANCE().isFLAT_EPSILON_PRIOR()) {
            this.maximizeEps();
        }
    }

    private void biasMu() {
        for (int j = 0; j < L; j++) {
            boolean flat = false;
            for (int k = 0; k < K; k++) {
                double max = 0;
                double sum = 0;
                for (int v = 0; v < n; v++) {
                    max = Math.max(this.mu[j][k][v], max);
                    sum += this.mu[j][k][v];
                }
                if (max < sum) {
                    flat = true;
                    break;
                }
            }
            if (flat) {
                for (int k = 0; k < K; k++) {
                    double sum = 0;
                    for (int v = 0; v < n; v++) {
                        this.mu[j][k][v] += Math.random() / 10d;
                        sum += this.mu[j][k][v];
                    }
                    for (int v = 0; v < n; v++) {
                        this.mu[j][k][v] /= sum;
                    }
                }
            }
        }
    }

    private void maximizeEps() {
        double a = 20;
        double b = 2357;//(-a * ew + a + 2 * ew - 1) / ew;//double ew = .008;
        for (int j = 0; j < L; j++) {
            this.eps[j] = Regularizations.f(this.nneqPos[j] + a) / Regularizations.f((coverage[j] * (n - 1)) + a + b);
            if (this.eps[j] > 1d / n) {
                this.eps[j] = 0.05;
            }
            this.antieps[j] = 1 - (n - 1) * eps[j];
        }
    }

    public void restart() {
        this.restart++;
        this.muChanged = 0;
        this.rhoChanged = 0;
        this.compute();
    }

    public Triplet<Integer, Integer, Double> minKL() {
        Set<Pair<Integer, Integer>> a = new HashSet<>();
        double min = Double.MAX_VALUE;
        Pair<Integer, Integer> argMin = null;
        for (int k = 0; k < K; k++) {
            for (int l = 0; l < K; l++) {
                if (k != l) {
                    if (!a.contains(Pair.with(k, l)) && !a.contains(Pair.with(l, k))) {
                        a.add(Pair.with(k, l));
                        double d = KullbackLeibler.symmetric(mu, k, l);
                        if (d < min) {
                            min = d;
                            argMin = Pair.with(k, l);
                        }
                    }
                }
            }
        }
        return Triplet.with(argMin.getValue0(), argMin.getValue1(), min);
    }

    protected void changedMu(double a, double b) {
        if (Math.abs(a - b) > Globals.getINSTANCE().getPCHANGE()) {
            this.muChanged++;
        }
    }

    protected void changedRho(double a, double b) {
        if (Math.abs(a - b) > Globals.getINSTANCE().getPCHANGE()) {
            this.rhoChanged++;
        }
    }

    protected final void prepare(Read[] reads, int N, int L, int K, int n, double[][][] rho, double[] pi, double[][][] mu) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.allReads = reads;
        this.rho = rho;
        this.rho_old = rho;
        this.mu = mu;
        this.mu_old = mu;
        this.pi = pi;
        this.snv = new double[L][n];
        this.muPrior = new double[n];
        for (int i = 0; i < n; i++) {
            this.muPrior[i] = Globals.getINSTANCE().getALPHA_H();
        }

        this.coverage = Globals.getINSTANCE().getTAU_OMEGA().getCoverage();
        this.paired = Globals.getINSTANCE().isPAIRED();
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

    public int getMuChanged() {
        return muChanged;
    }

    public int getRhoChanged() {
        return rhoChanged;
    }

    public double[][] getSnv() {
        return this.snv;
    }
}
