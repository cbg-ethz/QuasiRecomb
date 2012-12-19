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
import ch.ethz.bsse.quasirecomb.informationholder.Globals;

/**
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class Random {

    public static Dirichlet muDir;
    private static Dirichlet[] rhoDir;

    public static double[][][] generateInitRho(int Ldec, int K) {
        if (rhoDir == null || rhoDir.length != K) {
            rhoDir = new Dirichlet[K];
            for (int k = 0; k < K; k++) {
                double[] d = new double[K];
                for (int l = 0; l < K; l++) {
                    d[l] = k == l ? 1 : 0.001;
                }
                rhoDir[k] = new Dirichlet(d);
            }
        }
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

    public static double[] generateInitPi(int L, int K) {
//        if (Globals.getINSTANCE().isNO_RECOMB()) {
        return new Dirichlet(K, 2).nextDistribution();
//        } else {
//            double[] pi = new double[K];
//            for (int k = 1; k <= K; k++) {
//                pi[k - 1] = 1d / K;
//            }
//            return pi;
//        }
    }

    public static double[][][] generateMuInit(int L, int K, int n) {
        if (muDir == null) {
            muDir = new Dirichlet(n, 2);
        }
        double[][][] mu = new double[L][K][n];

        for (int j = L - 1; j >= 0; j--) {
            for (int k = K - 1; k >= 0; k--) {
                if (Globals.getINSTANCE().isPRIORMU()) {
                    double[] d = new double[n];
                    for (int v = 0; v < n; v++) {
//                        d[v] = (double) Globals.getINSTANCE().getMU_PRIOR()[j][v] + 1;
                        d[v] = 1d/n;
                    }
                    mu[j][k] = new Dirichlet(d).nextDistribution();
                } else {
                    mu[j][k] = muDir.nextDistribution();
                }
            }
        }
        return mu;
    }
}
