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
import java.util.Arrays;
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
public class JHMM extends JHMMBasics {

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
        this.compute();
    }

    public JHMM(Read[] reads, int N, int L, int K, int n, double[] eps, double[][][] rho, double[] pi, double[][][] mu) {
        this.eps = eps;
        this.antieps = new double[L];
        for (int j = 0; j < L; j++) {
            this.antieps[j] = (1 - (n - 1) * eps[j]);
        }
        this.prepare(reads, N, L, K, n, rho, pi, mu);
        this.compute();
    }

    public void restart() {
        this.restart++;
        this.parametersChanged = 0;
        this.compute();
    }

    private void compute() {
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nneqPos = new double[L];
        Globals.getINSTANCE().resetTimer();
        this.eStep();
        this.mStep();
    }

    private void eStep() {
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

    private void maximizeMu() {
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

    private void maximizeRho() {
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

    private void maximizePi() {
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
            this.maximizeRho();
        }
        this.maximizePi();
        this.maximizeMu();
        if (!Globals.getINSTANCE().isFLAT_EPSILON_PRIOR()) {
            this.maximizeEps();
        }
    }

    private void maximizeEps() {
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

    @Override
    protected JHMM clone() {
        double[][][] muClone = new double[L][K][n];
        double[][][] rhoClone = new double[L][K][K];
        double[] epsClone = Arrays.copyOf(eps, eps.length);
        double[] piClone = Arrays.copyOf(pi, pi.length);
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                System.arraycopy(this.mu[j][k], 0, muClone[j][k], 0, n);
                System.arraycopy(this.rho[j][k], 0, rhoClone[j][k], 0, L);

            }
        }
        return new JHMM(reads, N, L, K, n, epsClone, rhoClone, piClone, muClone);
    }
}
