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

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.informationholder.ParallelJHMMStorage;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadHMMStatic_NR {

    public static double computeFB(JHMM jhmm, Read read, ParallelJHMMStorage storage) {
        try {
            int begin = read.getBegin();

            int length = read.getLength();
            int K = jhmm.getK();
            int n = jhmm.getn();
            double[][][] fJKV = new double[length][K][n];
            double[][] fJK = new double[length][K];
            double[][] bJK = new double[length][K];
            double[] c = new double[length];
            double likelihood = 0;

            /*Forward*/
            for (int j = 0; j < length; j++) {
                final boolean hit = read.isHit(j);
                byte b = -1;
                double q = 1;
                if (hit) {
                    b = read.getBase(j);
                    q = read.getQuality(j);
                }
                int jGlobal = j + begin;

//                if (jGlobal < Globals.getINSTANCE().getTRAIN_WINDOW_BEGIN() || jGlobal >= Globals.getINSTANCE().getTRAIN_WINDOW_END()) {
//                    free(jhmm, storage);
//                    return 0;
//                }
                for (int k = 0; k < K; k++) {

                    for (int v = 0; v < n; v++) {
                        if (j == 0) {
                            fJKV[j][k][v] = jhmm.getPi()[jGlobal][k];
                        } else {
                            fJKV[j][k][v] = fJK[j - 1][k];
                        }

                        if (hit) {
                            if (q < 1) {
//                                fJKV[j][k][v] *= (b == v ? jhmm.getAntieps()[jGlobal] * (1 - (jhmm.getn() - 1) * (1 - q)) : jhmm.getEps()[jGlobal] * (1 - q));
                                fJKV[j][k][v] *= (b == v ? jhmm.getAntieps()[jGlobal] * q : jhmm.getEps()[jGlobal] * (1 - q) / (double) n);
                            } else {
                                fJKV[j][k][v] *= (b == v ? jhmm.getAntieps()[jGlobal] : jhmm.getEps()[jGlobal]);
                            }
                            fJKV[j][k][v] *= jhmm.getMu()[jGlobal][k][v];
//                            fJKV[j][k][v] *= q;
                        }
                        c[j] += fJKV[j][k][v];
                    }
                }
                c[j] = 1d / c[j];
                for (int k = 0; k < K; k++) {
                    for (int v = 0; v < n; v++) {
                        fJKV[j][k][v] *= c[j];
                        fJK[j][k] += fJKV[j][k][v];
                    }
                }
            }
            /*Backward*/
            for (int j = length - 1; j >= 0; j--) {
                final boolean hit = read.isHit(j);
                final boolean hit1 = j < (length - 1) && read.isHit(j + 1);
                if (hit) {
                    likelihood += Math.log(1d / c[j]);
                }
                int jGlobal = j + begin;
                double gammaSum = 0d;
                byte b = -1;
                double q = 1;
                if (hit1) {
                    b = read.getBase(j + 1);
                    q = read.getQuality(j + 1);
                }
                for (int k = 0; k < K; k++) {
                    if (j == length - 1) {
                        bJK[j][k] = c[j];
                    } else {
                        bJK[j][k] = 0;
                        if (hit1) {
                            double sumV = 0d;
                            for (int v = 0; v < n; v++) {
                                if (q < 1) {
//                                        sumV += (b == v ? jhmm.getAntieps()[jGlobal + 1] * (1 - (jhmm.getn() - 1) * (1 - q)) : jhmm.getEps()[jGlobal + 1] * (1 - q)) * jhmm.getMu()[jGlobal + 1][l][v];
                                    sumV += (b == v ? jhmm.getAntieps()[jGlobal + 1] * q : jhmm.getEps()[jGlobal + 1] * (1 - q) / (double) n) * jhmm.getMu()[jGlobal + 1][k][v];
                                } else {
                                    sumV += (b == v ? jhmm.getAntieps()[jGlobal + 1] : jhmm.getEps()[jGlobal + 1]) * jhmm.getMu()[jGlobal + 1][k][v];
                                }
                            }
//                                sumV *= q;
                            bJK[j][k] += sumV * bJK[j + 1][k];
                        } else {
                            bJK[j][k] += bJK[j + 1][k];
                        }
                        bJK[j][k] *= c[j];
                    }
                    if (Double.isInfinite(bJK[j][k])) {
                        //this is infinite, because the char has not been observed and there no probability to emit it
                        //thus we divide 0 by a very small number, i.e. 1e-300.
                        bJK[j][k] = 0d;
                    }
                    gammaSum += fJK[j][k] * bJK[j][k];
                }
                if (hit) {
                    b = read.getBase(j);
                    double xiSum = 0d;
                    for (int k = 0; k < K; k++) {
                        if (gammaSum == 0) {
                            for (int v = 0; v < n; v++) {
                                storage.addnJKV(jGlobal, k, v, ((double) read.getCount()) / n);
                            }
                        } else {
                            for (int v = 0; v < n; v++) {
                                double gamma = read.getCount() * fJKV[j][k][v] * bJK[j][k] / gammaSum;
                                storage.addnJKV(jGlobal, k, v, gamma);
                                if (b != v) {
                                    storage.addnneqPos(j, gamma);
                                }
                            }
                        }
                        if (j > 0) {
                            for (int l = 0; l < K; l++) {
                                double marginalV = 0d;
                                for (int v = 0; v < n; v++) {
                                    marginalV += (b == v ? jhmm.getAntieps()[jGlobal] : jhmm.getEps()[jGlobal]) * jhmm.getMu()[jGlobal][l][v];
                                }
                                double xi = fJK[j - 1][k] * jhmm.getRho()[jGlobal - 1][k][l] * marginalV * bJK[j][l];
                                xiSum += xi;
                            }
                        }
                    }
                    for (int k = 0; k < K; k++) {
                        if (xiSum != 0 && j > 0) {
                            for (int l = 0; l < K; l++) {
                                double marginalV = 0d;
                                for (int v = 0; v < n; v++) {
                                    marginalV += (b == v ? jhmm.getAntieps()[jGlobal] : jhmm.getEps()[jGlobal]) * jhmm.getMu()[jGlobal][l][v];
                                }
                                double xi = read.getCount() * fJK[j - 1][k] * jhmm.getRho()[jGlobal - 1][k][l] * marginalV * bJK[j][l] / xiSum;
                                storage.addnJKL(jGlobal, k, l, xi);
                            }
                        }
                    }
                }
            }
            likelihood *= read.getCount();

            return likelihood;
        } catch (Exception e) {
            System.err.println(e);
        }
        return 0;
    }
}
