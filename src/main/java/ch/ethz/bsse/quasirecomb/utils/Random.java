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
package ch.ethz.bsse.quasirecomb.utils;

import cc.mallet.types.Dirichlet;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;

/**
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class Random {

    private static Dirichlet[] rhoDir;
    public static Dirichlet muDir;
    private static Dirichlet[] sigmaDir;

    public static double[][][][] generateInitSigma(int Ldec, int K) {
        if (sigmaDir == null) {
            sigmaDir[0] = new Dirichlet(new double[]{1, 0.01});
            sigmaDir[1] = new Dirichlet(new double[]{0.01, 1});
        }
        double sigma[][][][] = new double[Ldec][K][2][2];
        for (int j = 0; j < Ldec; j++) {
            for (int k = 0; k < K; k++) {
                for (int a = 0; a < 2; a++) {
                    sigma[j][k][a] = sigmaDir[a].nextDistribution();
                }
            }
        }
        return sigma;
    }

    public static double[][][] generateInitRho(int Ldec, int K) {
//        if (rhoDir == null || rhoDir.length != K) {
        rhoDir = new Dirichlet[K];
        for (int k = 0; k < K; k++) {
            double[] d = new double[K];
            for (int l = 0; l < K; l++) {
                d[l] = k == l ? 1 : 0.01;
            }
            rhoDir[k] = new Dirichlet(d);
        }
//        }
        double[][][] rho = new double[Ldec][K][K];
        if (!Globals.getINSTANCE().isNO_RECOMB()) {
            for (int j = 0; j < Ldec; j++) {
                for (int k = 0; k < K; k++) {
                    rho[j][k] = rhoDir[k].nextDistribution();
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

    public static double[][] generateInitPi(int L, int K) {
        if (Globals.getINSTANCE().isNO_RECOMB()) {
            double[][] pis = new double[L][K];
            for (int j = 0; j < L; j++) {
                pis[j] = new Dirichlet(K, 2).nextDistribution();
            }
            return pis;
        } else {
            double[][] pi = new double[L][K];
            for (int j = 0; j < L; j++) {
                for (int k = 1; k <= K; k++) {
                    pi[j][k - 1] = 1d / K;
                }
            }
            return pi;
        }
    }

    public static double[][][] generateMuInit(int L, int K, int n) {
//        if (muDir == null) {
        muDir = new Dirichlet(n, 100);
//        }
        double[][][] mu = new double[L][K][n];

        for (int j = L - 1; j >= 0; j--) {
            for (int k = K - 1; k >= 0; k--) {

//                if (Globals.getINSTANCE().isPRIORMU()) {
                double[] d = new double[n];
                for (int v = 0; v < n; v++) {
//                        d[v] = (double) Globals.getINSTANCE().getMU_PRIOR()[j][v] + 1;
                    d[v] = 1d / n;
                }
                if (Globals.getINSTANCE().isNO_RECOMB()) {
                    mu[j][k] = muDir.nextDistribution();
                } else {
                    for (int v = 0; v < n; v++) {
                        mu[j][k][v] = 1d / n;
                    }
                }
//                } else {
//                    mu[j][k] = muDir.nextDistribution();
//                }
            }
        }
        return mu;
    }
}
