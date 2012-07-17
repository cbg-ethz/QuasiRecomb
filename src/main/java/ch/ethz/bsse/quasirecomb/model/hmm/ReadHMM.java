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

/**
 * Calculates the forward / backward algorithm for a single read.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadHMM {

    private final int L;
    private final int K;
    private final int n;
    private final Read read;
    private double[][][] rho;
    private double[] pi;
    private double[][][] mu;
    private double[] eps;
    private double[] antieps;
    private double[][] fJK;
    private double[][] bJK;
    private double[] c;
    private double[][] gJK;
    private double[][][] gJKV;
    private double[][][] xJKL;
    private double[][][] fJKV;
    private int begin;
    private int end;
    private int length;

    public ReadHMM(int L, int K, int n, Read read, double[][][] rho, double[] pi, double[][][] mu, double[] eps, double[] antieps) {
        this.L = L;
        this.K = K;
        this.n = n;
        this.read = read;
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.eps = eps;
        this.antieps = antieps;
        this.begin = read.getBegin();
        this.end = read.getEnd();
        this.length = end - begin;

        calculate();
    }

    private void calculate() {
        this.forward();
        this.backward();
        this.gammaXsi();
    }

    public boolean checkConsistency() {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                if (Double.isNaN(this.fJK[j][k])) {
                    return true;
                }
                for (int v = 0; v < n; v++) {
                    if (Double.isNaN(this.fJKV[j][k][v])) {
                        return true;
                    }
                }
            }
            if (Double.isNaN(this.c[j])) {
                return true;
            }
        }
        for (int j = 0; j < L - 1; j++) {
            for (int k = 0; k < K; k++) {
                if (Double.isNaN(this.bJK[j][k])) {
                    return true;
                }
            }
        }
        return true;
    }

    public void recalc(double[][][] rho, double[] pi, double[][][] mu, double[] epsilon) {
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.eps = epsilon;
        for (int j = 0; j < L; j++) {
            this.antieps[j] = 1 - (n - 1) * epsilon[j];
        }
        this.calculate();
    }

    private void forward() {
        this.c = new double[length];
        this.fJKV = new double[length][K][n];
        this.fJK = new double[length][K];
        for (int j = 0; j < length; j++) {
            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
                    fJKV[j][k][v] = prRjHv(begin + j, v, k) * mu[begin + j][k][v];
                    if (j == 0) {
                        fJKV[j][k][v] *= pi[k];
                    } else {
                        double sumL = 0d;
                        for (int l = 0; l < K; l++) {
                            sumL += fJK[j - 1][l] * rho[begin + j - 1][l][k];
                        }
                        fJKV[j][k][v] *= sumL;
                    }
                    c[j] += fJKV[j][k][v];
                    if (fJKV[j][k][v] < 0) {
                        System.out.println("fjkv");
                    }
                }
            }

            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
                    fJKV[j][k][v] /= c[j];
                    fJK[j][k] += fJKV[j][k][v];
                }
            }
        }
    }

    private void backward() {
        this.bJK = new double[length][K];
        for (int j = length - 1; j >= 0; j--) {
            for (int k = 0; k < K; k++) {
                if (j == length - 1) {
                    bJK[j][k] = 1d / c[length - 1];
                } else {
                    for (int l = 0; l < K; l++) {
                        double sumV = 0d;
                        for (int v = 0; v < n; v++) {
                            sumV += prRjHv(begin + j + 1, v, k) * mu[begin + j + 1][l][v];
                        }
                        bJK[j][k] += sumV * rho[begin + j][k][l] * bJK[j + 1][l];
                    }
                    if (c[j] != 0d) {
                        bJK[j][k] /= c[j];
                    }
                }
            }
        }
    }

    private void gammaXsi() {
        this.gJK = new double[length][K];
        this.gJKV = new double[length][K][n];
        this.xJKL = new double[length][K][K];
        for (int j = 0; j < length; j++) {
            for (int k = 0; k < K; k++) {
                this.gJK[j][k] = this.fJK[j][k] * this.bJK[j][k] * c[j];
                for (int v = 0; v < n; v++) {
                    this.gJKV[j][k][v] = this.fJKV[j][k][v] * this.bJK[j][k] * c[j];
                }
                this.gJK[j][k] = this.fJK[j][k] * c[j] * this.bJK[j][k];
                if (j > 0) {
                    for (int l = 0; l < K; l++) {
                        double marginalV = 0d;
                        for (int v = 0; v < n; v++) {
                            marginalV += prRjHv(j, v, k) * mu[begin + j][l][v];
                        }
                        this.xJKL[j][k][l] = this.fJK[j - 1][k] * marginalV * rho[begin + j - 1][k][l] * this.bJK[j][l];
                    }
                }
            }
        }
    }

    public double gamma(int j, int k) {
        if (j >= begin && j - begin < length) {
            return this.gJK[j - begin][k];
        }
        return 0;
    }

    final public double gamma(int j, int k, int v) {
        if (j >= begin && j - begin < length) {
            return this.gJKV[j - begin][k][v];
        }
        return 0;
    }

    final public double xi(int j, int k, int l) {
        if (j < 1) {
            throw new IllegalStateException("J >= 1?");
        }
        if (j >= begin && j - begin < length) {
            return this.xJKL[j - begin][k][l];
        }
        return 0;
    }

    private double prRjHv(int j, int v, int k) {
        return (this.getSequence()[j] == v) ? antieps[j] : eps[j];
    }

    final public double getC(int j) {
        if (j >= begin && j - begin < length) {
            return c[j - begin];
        }
        return 1;
    }

    final public byte[] getSequence() {
        return read.getSequence();
    }

    public int getCount() {
        return this.read.getCount();
    }
}
