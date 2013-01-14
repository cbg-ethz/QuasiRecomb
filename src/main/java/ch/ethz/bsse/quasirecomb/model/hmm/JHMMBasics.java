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
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class JHMMBasics {

    protected int Kmin;
    protected int N;
    protected int L;
    protected int K;
    protected int n;
    //rho[j][k][l] := transition prob. at position j, for l given k
    protected double[][] tauOmega;
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
//    protected Read[] currentReads;
    protected Read[] allReads;
    protected int restart = 0;
    protected int coverage[];
    protected int muChanged = 0;
    protected int rhoChanged = 0;
    protected boolean paired;
    protected Map<Integer, TempJHMMStorage> garage = new ConcurrentHashMap<>();
    protected final List<Integer> available = new ArrayList<>();
    int s = 0;

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

    protected void prepare(Read[] reads, int N, int L, int K, int n, double[][][] rho, double[] pi, double[][][] mu) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.allReads = reads;
        this.rho = rho;
        this.mu = mu;
        this.pi = pi;

        this.muPrior = new double[n];
        for (int i = 0; i < n; i++) {
            this.muPrior[i] = Globals.getINSTANCE().getALPHA_H();
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
    }

    protected void init() {
        int[] tau1 = new int[L + 1];
        int[] tau2 = new int[L + 1];
        int[] omega1 = new int[L + 1];
        int[] omega2 = new int[L + 1];
        for (Read r : allReads) {
            for (int i = r.getWatsonBegin(); i < r.getWatsonEnd(); i++) {
                this.coverage[i] += r.getCount();
            }
            tau1[r.getWatsonBegin()] += r.getCount();
            omega1[r.getWatsonEnd()] += r.getCount();
            if (r.isPaired()) {
                for (int i = r.getCrickBegin(); i < r.getCrickEnd(); i++) {
                    this.coverage[i] += r.getCount();
                }
                tau2[r.getCrickBegin()] += r.getCount();
                omega2[r.getCrickEnd()] += r.getCount();
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
        if (Globals.getINSTANCE().isDEBUG()) {
            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "twtw", sb.toString());
        }
    }

    public TempJHMMStorage getStorage() {
        synchronized (this.available) {
            while (!available.iterator().hasNext()) {
                try {
                    notify();
                    TimeUnit.MILLISECONDS.sleep(10);
                    System.err.println("sleep");
                } catch (InterruptedException ex) {
                    Logger.getLogger(JHMMBasics.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
//            System.out.println("GET " + s++);
            Integer i = available.iterator().next();
            available.remove(i);
            return garage.get(i);
        }
    }

    public void free(int id) {
        synchronized (this.available) {
            this.available.add(id);
        }
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
}
