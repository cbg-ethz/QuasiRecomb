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

import ch.ethz.bsse.quasirecomb.distance.KullbackLeibler;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.CallableReadHMM;
import ch.ethz.bsse.quasirecomb.utils.Random;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
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
public class JHMM extends JHMMBasics {

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

    public void restart() {
        this.restart++;
        this.muChanged = 0;
        this.rhoChanged = 0;
        this.compute();
    }
    int currentIndex = 0;

//    private void prepareReads() {
//        final int readCount = this.allReads.length;
//        final int stepSize = Globals.getINSTANCE().getSTEPS();
//
//        if (currentIndex >= 0 && currentIndex < readCount) {
//            if (currentIndex + stepSize >= readCount) {
//                currentReads = new Read[readCount - currentIndex];
//                for (int j = 0; j < readCount - currentIndex; j++) {
//                    currentReads[j] = allReads[currentIndex + j];
//                }
//            } else {
//                currentReads = new Read[stepSize];
//                for (int j = 0; j < stepSize; j++) {
//                    currentReads[j] = allReads[currentIndex + j];
//                }
//            }
//            currentIndex += stepSize;
//        } else {
//            Collections.shuffle(Arrays.asList(this.allReads));
//            currentIndex = 0;
//        }
//    }
    private void compute() {
//        prepareReads();
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nneqPos = new double[L];
        Globals.getINSTANCE().resetTimer();
        this.eStep();
        this.mStep();
        s++;
    }
    int s = 0;

    private void clearGarage() {
        if (Globals.getINSTANCE().isSTORAGE()) {
            available.clear();
            garage.clear();
            for (int i = 0; i < Globals.getINSTANCE().getCpus(); i++) {
                available.add(i);
                garage.put(i, new TempJHMMStorage(L, K, n, i));
            }
        }
    }

    private void mergeGarage() {
        if (Globals.getINSTANCE().isSTORAGE()) {
            if (available.size() != Globals.getINSTANCE().getCpus()) {
                throw new IllegalStateException("Not all storages have been returned");
            }
            Iterator<TempJHMMStorage> iterator = this.garage.values().iterator();
            TempJHMMStorage store = iterator.next();
            while (iterator.hasNext()) {
                store.merge(iterator.next());
            }
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

    private void eStep() {
        clearGarage();
        this.loglikelihood = 0d;
        List<Future<Double>> results = new ArrayList<>();

        for (int i = 0; i < allReads.length; i++) {
            results.add(Globals.getINSTANCE().getExecutor().submit(new CallableReadHMM(this, allReads[i])));
            Globals.getINSTANCE().printPercentage(K, (double) i / allReads.length, Kmin);
        }
        Globals.getINSTANCE().getExecutor().shutdown();
        try {
            while (!Globals.getINSTANCE().getExecutor().awaitTermination(1, TimeUnit.SECONDS)) {
                TimeUnit.MILLISECONDS.sleep(10);
                System.out.println("sleeping");
            }
        } catch (InterruptedException e) {
            System.err.println("e");
            Utils.error();
            System.exit(0);
        }

        for (int i = 0; i < results.size(); i++) {
            try {
                Double llh = results.get(i).get();
                loglikelihood += llh;
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(JHMM.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        mergeGarage();
        Globals.renewExecutor();
    }

//    private void maximizeStep() {
//        for (int j = 0; j < L; j++) {
//            for (int k = 0; k < K; k++) {
//                double eta = Math.pow(s++ + 2, -Globals.getINSTANCE().getINTERPOLATE_MU());
//                this.mu[j][k] = Regularizations.step(this.nJKV[j][k], this.mu[j][k], eta);
//                if (j > 0) {
//                    eta = Math.pow(s++ + 2, -Globals.getINSTANCE().getINTERPOLATE_RHO());
//                    this.rho[j - 1][k] = Regularizations.step(this.nJKL[j][k], this.rho[j - 1][k], eta);
//                }
//            }
//        }
//    }
//    private boolean iterationCompleted() {
//        return currentIndex >= this.allReads.length;
//    }
    private void maximizeMu() {
        double[] muJKV;
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
//                if (Globals.getINSTANCE().isML()) {
//                    muJKV = Regularizations.ml(this.nJKV[j][k]);
                double eta = Math.pow(s + 2, -Globals.getINSTANCE().getINTERPOLATE_MU());
                muJKV = Regularizations.step(this.nJKV[j][k], this.mu[j][k], eta);
//                } else {
//                if (iterationCompleted()) {
                double mult = Globals.getINSTANCE().getMULT_MU();
                double[] muPriorLocal = new double[n];
                System.arraycopy(this.muPrior, 0, muPriorLocal, 0, n);
                boolean repeat = true;
                do {
                    muJKV = Regularizations.regularizeOnce(muJKV, restart, muPriorLocal, mult);
                    double prev = muJKV[0];
                    for (int v = 1; v < n; v++) {
                        if (prev != muJKV[v]) {
                            repeat = false;
                        }
                    }
                    if (repeat) {
                        mult *= 2;
                        for (int v = 0; v < n; v++) {
                            muPriorLocal[v] *= 10;
//                            repeats++;
                        }
                    }
                    if (mult > 1000) {
                        repeat = false;
                    }
                } while (repeat);
//                }
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
//                }
            }
        }
    }

    private boolean maximizeRho() {
//        int repeats = 0;
        boolean forceRho = false;
        double[] rhoJKL = null;
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {
//                if (Globals.getINSTANCE().isML()) {
                double eta = Math.pow(s + 2, -Globals.getINSTANCE().getINTERPOLATE_RHO());
                rhoJKL = Regularizations.step(this.nJKL[j][k], this.rho[j - 1][k], eta);
//                rhoJKL = Regularizations.ml(this.nJKL[j][k]);
//                } else {
//                if (iterationCompleted()) {
//                System.out.println("BAM");
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
//                                if (mult > 10) {
//                                    mult = 10;
//                                }
                            rhoPrior[l] = 100;
                            fix = true;
                            forceRho = true;
                            break;
//                                repeats++;
                        } else {
                            rhoPrior[l] = Globals.getINSTANCE().getALPHA_Z() * 10;
                        }
                    } else {
                        rhoPrior[l] = Globals.getINSTANCE().getALPHA_Z();
                    }
                }
//                rhoJKL = null;
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
//            }
//        }
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
        pi = Regularizations.step(piTmp, pi, eta);
    }
    private int oldFlatMu = -1;
    private boolean biasMu = false;
    private int biasCounter = 0;
    private int unBiasCounter = 0;

    private void mStep() {
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
//            }
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
//                    max = Math.max(this.nJKV[j][k][v], max);
//                    sum += this.nJKV[j][k][v];
                }
                if (max < sum) {
                    flat = true;
                    break;
                }
            }
            if (flat) {
                for (int k = 0; k < K; k++) {
//                    this.mu[j][k] = Random.muDir.nextDistribution();
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
        double ew = .008;
//        if (Globals.getINSTANCE().isUNINFORMATIVE_EPSILON_PRIOR()) {
//            for (int j = 0; j < L; j++) {
//                if (this.nneqPos[j] != 0d) {
//                    ew += this.nneqPos[j] / N;
//                }
//            }
//            if (ew != 0d) {
//                ew /= L;
//            }
//        } else {
//            ew = 
//        }
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

//    @Override
//    protected JHMM clone() {
//        double[][][] muClone = new double[L][K][n];
//        double[][][] rhoClone = new double[L][K][K];
//        double[] epsClone = Arrays.copyOf(eps, eps.length);
//        double[] piClone = Arrays.copyOf(pi, pi.length);
//        for (int j = 0; j < L; j++) {
//            for (int k = 0; k < K; k++) {
//                System.arraycopy(this.mu[j][k], 0, muClone[j][k], 0, n);
//                System.arraycopy(this.rho[j][k], 0, rhoClone[j][k], 0, L);
//
//            }
//        }
//        return new JHMM(reads, N, L, K, n, epsClone, rhoClone, piClone, muClone, Kmin);
//    }
    public void merge() {
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
}
