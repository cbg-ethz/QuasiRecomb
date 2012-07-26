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
package ch.ethz.bsse.quasirecomb.utils;

import cc.mallet.types.Dirichlet;
import ch.ethz.bsse.quasirecomb.model.Globals;

/**
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class Random {

    public static double[][][] generateInitRho(int Ldec, int K) {
        double[][][] rho = new double[Ldec][K][K];
        if (!Globals.NO_RECOMB) {
            for (int j = 0; j < Ldec; j++) {
                for (int k = 0; k < K; k++) {
                    rho[j][k] = new Dirichlet(K, Globals.BETA_Z).nextDistribution();
                    int maxIndex = 0;
                    double max = 0;
                    for (int l = 0; l < K; l++) {
                        if (rho[j][k][l] > max) {
                            max = rho[j][k][l];
                            maxIndex = l;
                        }
                    }
                    rho[j][k][maxIndex] = rho[j][k][k];
                    rho[j][k][k] = max;
                }
            }
        } else {
            for (int j = 0; j < Ldec; j++) {
                for (int k = 0; k < K; k++) {
                    for (int l = 0; l < K; l++) {
                        rho[j][k][l] = (l == k) ? 1d : 0d;
                    }
                }
            }
        }
        return rho;
    }

    public static double[] generateInitPi(int L, int K) {
        if (Globals.NO_RECOMB) {
            return new Dirichlet(K, 2).nextDistribution();
        } else {
            double[] pi = new double[K];
            for (int k = 1; k <= K; k++) {
                pi[k - 1] = 1d / K;
            }
            return pi;
        }
    }

    public static double[][][] generateMuInit(int L, int K, int n) {
        double[][][] mu = new double[L][K][n];

        for (int j = L - 1; j >= 0; j--) {
            for (int k = K - 1; k >= 0; k--) {
                if (Globals.NO_RECOMB) {
                    mu[j][k] = new Dirichlet(n, 2).nextDistribution();
                } else {
                    for (int i = n - 1; i >= 0; i--) {
                        mu[j][k][i] = 1d / n;
                    }
                }
            }
        }
        return mu;
    }
}
