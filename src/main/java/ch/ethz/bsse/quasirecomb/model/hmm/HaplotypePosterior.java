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

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.informationholder.ParallelJHMMStorage;
import ch.ethz.bsse.quasirecomb.utils.Utils;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class HaplotypePosterior {

    public static double computeFB(OptimalResult or, String read) {
        try {
            int begin = 0;

            int length = read.length();
            byte[] readBases = Utils.splitReadIntoByteArray(read);
            double[][][] fJKV = new double[length][or.getK()][or.getn()];
            double[][] fJK = new double[length][or.getK()];
            double[][] bJK = new double[length][or.getK()];
            double[] c = new double[length];
            double likelihood = 0;

            /*Forward*/
            for (int j = 0; j < length; j++) {
                byte b = -1;
                b = readBases[j];
                int jGlobal = j + begin;
                for (int k = 0; k < or.getK(); k++) {
                    for (int v = 0; v < or.getn(); v++) {
                        if (j == 0) {
                            fJKV[j][k][v] = or.getPi()[jGlobal][k];
                        } else {
                            double sumL = 0d;
                            for (int l = 0; l < or.getK(); l++) {
                                sumL += fJK[j - 1][l] * or.getRho()[jGlobal - 1][l][k];
                            }
                            fJKV[j][k][v] = sumL;
                        }

                        fJKV[j][k][v] *= (b == v ? (1 - (or.getn()) * or.getEps()[jGlobal]) : or.getEps()[jGlobal]);
                        fJKV[j][k][v] *= or.getMu()[jGlobal][k][v];
//                            fJKV[j][k][v] *= q;
                        c[j] += fJKV[j][k][v];
                    }
                }
                c[j] = 1d / c[j];
                for (int k = 0; k < or.getK(); k++) {
                    for (int v = 0; v < or.getn(); v++) {
                        fJKV[j][k][v] *= c[j];
                        fJK[j][k] += fJKV[j][k][v];
                    }
                }
            }
            /*Backward*/
            for (int j = length - 1; j >= 0; j--) {
                likelihood += Math.log(1d / c[j]);
                int jGlobal = j + begin;
                byte b = -1;
                b = readBases[j + 1];
                for (int k = 0; k < or.getK(); k++) {
                    if (j == length - 1) {
                        bJK[j][k] = c[j];
                    } else {
                        bJK[j][k] = 0;
                        for (int l = 0; l < or.getK(); l++) {
                            double sumV = 0d;
                            for (int v = 0; v < or.getn(); v++) {
                                sumV += (b == v ? (1 - (or.getn()) * or.getEps()[jGlobal + 1]) : or.getEps()[jGlobal + 1]) * or.getMu()[jGlobal + 1][l][v];
                            }
                            bJK[j][k] += sumV * or.getRho()[jGlobal][k][l] * bJK[j + 1][l];
                        }
                        bJK[j][k] *= c[j];
                    }
                    if (Double.isInfinite(bJK[j][k])) {
                        //this is infinite, because the char has not been observed and there no probability to emit it
                        //thus we divide 0 by a very small number, i.e. 1e-300.
                        bJK[j][k] = 0d;
                    }
                }
            }
            return likelihood;
        } catch (Exception e) {
            System.err.println(e);
        }
        return 0;
    }
}
