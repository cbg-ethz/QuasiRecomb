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
import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadHMMStatic {

    public static double computeFB(JHMM jhmm, Read read) {
        try {
        TempJHMMStorage storage = jhmm.getStorage();
        int begin = read.getBegin();
        int length = read.getLength();
        double[][][] fJKV = new double[length][jhmm.getK()][jhmm.getn()];
        double[][] fJK = new double[length][jhmm.getK()];
        double[][] bJK = new double[length][jhmm.getK()];
        double[] c = new double[length];
        double likelihood = 0;

        for (int j = 0; j < length; j++) {
            int jGlobal = j + begin;
            for (int k = 0; k < jhmm.getK(); k++) {
                for (int v = 0; v < jhmm.getn(); v++) {
                    if (j == 0) {
                        fJKV[j][k][v] = jhmm.getPi()[k];
                    } else {
                        double sumL = 0d;
                        for (int l = 0; l < jhmm.getK(); l++) {
                            sumL += fJK[j - 1][l] * jhmm.getRho()[jGlobal - 1][l][k];
                        }
                        fJKV[j][k][v] = sumL;
                        if (Double.isNaN(fJKV[j][k][v])) {
                            System.err.println("x");
                        }
                    }

                    if (read.isHit(j)) {
                        fJKV[j][k][v] *= (read.getBase(j) == v ? jhmm.getAntieps()[jGlobal] : jhmm.getEps()[jGlobal]);
                        fJKV[j][k][v] *= jhmm.getMu()[jGlobal][k][v];
                    }
                    if (Double.isNaN(fJKV[j][k][v])) {
                        System.err.println("fJKV:\t" + fJKV[j][k][v]);
                        System.err.println("mu:\t" + jhmm.getMu()[jGlobal][k][v]);
                        System.exit(0);
                    }
                    c[j] += fJKV[j][k][v];
                }
            }
            if (c[j] <= 0) {
//                jhmm.free(storage.getId());
//                return 0;
                System.err.println("R");
            }
            c[j] = 1d / c[j];
            for (int k = 0; k < jhmm.getK(); k++) {
                for (int v = 0; v < jhmm.getn(); v++) {
                    fJKV[j][k][v] *= c[j];
                    fJK[j][k] += fJKV[j][k][v];
                }
            }
        }
        for (int j = length - 1; j >= 0; j--) {
            if (read.isHit(j)) {
                likelihood += Math.log(1d / c[j]);
            }
            int jGlobal = j + begin;
            double gammaSum = 0d;
            for (int k = 0; k < jhmm.getK(); k++) {
                if (j == length - 1) {
                    bJK[j][k] = c[j];
                } else {
                    bJK[j][k] = 0;
                    for (int l = 0; l < jhmm.getK(); l++) {
                        if (read.isHit(j + 1)) {
                            double sumV = 0d;
                            for (int v = 0; v < jhmm.getn(); v++) {
                                sumV += (read.getBase(j + 1) == v ? jhmm.getAntieps()[jGlobal + 1] : jhmm.getEps()[jGlobal + 1]) * jhmm.getMu()[jGlobal + 1][l][v];
                            }
                            bJK[j][k] += sumV * jhmm.getRho()[jGlobal][k][l] * bJK[j + 1][l];
                        } else {
                            bJK[j][k] += jhmm.getRho()[jGlobal][k][l] * bJK[j + 1][l];
                        }
                    }
                    if (c[j] == 0d) {
                        System.err.println("C == 0");
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
            if (read.isHit(j)) {
                double xiSum = 0d;
                for (int k = 0; k < jhmm.getK(); k++) {
                    if (Double.isNaN(gammaSum)) {
                        System.err.println("XXX");
                        System.exit(0);
                    } else if (gammaSum == 0) {
                        for (int v = 0; v < jhmm.getn(); v++) {
                            storage.addnJKV(jGlobal, k, v, ((double) read.getCount()) / jhmm.getn());
                        }
                    } else {
                        for (int v = 0; v < jhmm.getn(); v++) {
                            double gamma = read.getCount() * fJKV[j][k][v] * bJK[j][k] / gammaSum;
                            storage.addnJKV(jGlobal, k, v, gamma);
                            if (Double.isNaN(gamma)) {
                                System.out.println("#####");
                                System.exit(0);
                            }
                            if (read.getBase(j) != v) {
                                storage.addnneqPos(j, gamma);
                            }
                        }
                    }
                    if (j > 0) {
                        for (int l = 0; l < jhmm.getK(); l++) {
                            double marginalV = 0d;
                            for (int v = 0; v < jhmm.getn(); v++) {
                                marginalV += (read.getBase(j) == v ? jhmm.getAntieps()[jGlobal] : jhmm.getEps()[jGlobal]) * jhmm.getMu()[jGlobal][l][v];
                            }
                            double xi = fJK[j - 1][k] * jhmm.getRho()[jGlobal - 1][k][l] * marginalV * bJK[j][l];
                            xiSum += xi;
                        }
                    }
                }
                for (int k = 0; k < jhmm.getK(); k++) {
                    if (xiSum != 0 && j > 0) {
                        for (int l = 0; l < jhmm.getK(); l++) {
                            double marginalV = 0d;
                            for (int v = 0; v < jhmm.getn(); v++) {
                                marginalV += (read.getBase(j) == v ? jhmm.getAntieps()[jGlobal] : jhmm.getEps()[jGlobal]) * jhmm.getMu()[jGlobal][l][v];
                            }
                            double xi = read.getCount() * fJK[j - 1][k] * jhmm.getRho()[jGlobal - 1][k][l] * marginalV * bJK[j][l] / xiSum;
                            if (Double.isNaN(xi)) {
                                System.err.println("xi nan");
                                System.exit(0);
                            }
                            storage.addnJKL(jGlobal, k, l, xi);
                        }
                    }
                }
            }
        }

        likelihood *= read.getCount();

        free(jhmm,storage);
        return likelihood;
        } catch (Exception e) {
            System.err.println(e);
        }
        return 0;
    }

    private static void free(JHMM jhmm, TempJHMMStorage storage) {
        jhmm.free(storage.getId());
    }
}
