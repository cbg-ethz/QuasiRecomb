/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;

/**
 *
 * @author toepfera
 */
public class ReadHMMStatic {

    public static double computeFB(JHMMI jhmm, Read read) {
        TempJHMMStorage storage = jhmm.getStorage();
        int begin = read.getBegin() - Globals.getINSTANCE().getALIGNMENT_BEGIN();
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
                likelihood += Math.log(1d/c[j]);
            }
            int jGlobal = j + begin;
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
            }
        }

        for (int j = 0; j < length; j++) {
            int jGlobal = j + begin;
            if (read.isHit(j)) {
                double gammaSum = 0d;
                for (int k = 0; k < jhmm.getK(); k++) {
                    gammaSum += fJK[j][k] * bJK[j][k];
                }
                for (int k = 0; k < jhmm.getK(); k++) {
                    double gammaKSum = 0d;
                    for (int v = 0; v < jhmm.getn(); v++) {
                        double gamma = read.getCount() * fJKV[j][k][v] * bJK[j][k] / gammaSum;
                        gammaKSum += gamma;
                        storage.addnJKV(jGlobal, k, v, gamma);
                        if (Double.isNaN(gamma)) {
                            System.out.println("#####");
                        }
                        if (read.getBase(j) != v) {
                            storage.addnneqPos(j, gamma);
                        }
                    }
                    if (gammaKSum == 0) {
                        for (int v = 0; v < jhmm.getn(); v++) {
                            storage.addnJKV(jGlobal, k, v, ((double) read.getCount()) / jhmm.getn());
                        }
                    }
                }
            }
        }

        for (int j = 1; j < length; j++) {
            int jGlobal = j + begin;
            if (read.isHit(j)) {
                double xiSum = 0d;
                for (int k = 0; k < jhmm.getK(); k++) {
                    for (int l = 0; l < jhmm.getK(); l++) {
                        double marginalV = 0d;
                        for (int v = 0; v < jhmm.getn(); v++) {
                            marginalV += (read.getBase(j) == v ? jhmm.getAntieps()[jGlobal] : jhmm.getEps()[jGlobal]) * jhmm.getMu()[jGlobal][l][v];
                        }
                        double xi = fJK[j - 1][k] * jhmm.getRho()[jGlobal - 1][k][l] * marginalV * bJK[j][l];
                        xiSum += xi;
                    }
                }
                for (int k = 0; k < jhmm.getK(); k++) {
                    for (int l = 0; l < jhmm.getK(); l++) {
                        double marginalV = 0d;
                        for (int v = 0; v < jhmm.getn(); v++) {
                            marginalV += (read.getBase(j) == v ? jhmm.getAntieps()[jGlobal] : jhmm.getEps()[jGlobal]) * jhmm.getMu()[jGlobal][l][v];
                        }
                        double xi = read.getCount() * fJK[j - 1][k] * jhmm.getRho()[jGlobal - 1][k][l] * marginalV * bJK[j][l] / xiSum;
                        if (Double.isNaN(xi)) {
                            System.err.println("xx");
//                                xi = 0;
                        }
                        storage.addnJKL(jGlobal, k, l, xi);
                    }
                }
            }
        }
//        for (int j = 0; j < length - 2; j++) {
//            int jGlobal = j + begin;
//            double denom = 0d;
//            for (int k = 0; k < jhmm.getK(); k++) {
//                for (int l = 0; l < jhmm.getK(); l++) {
//                    double sumV = 0d;
//                    for (int v = 0; v < jhmm.getn(); v++) {
//                        sumV += (read.getBase(j + 1) == v ? jhmm.getAntieps()[jGlobal + 1] : jhmm.getEps()[jGlobal + 1]) * jhmm.getMu()[jGlobal + 1][l][v];
//                    }
//                    denom += fJK[j][k] * jhmm.getRho()[jGlobal][k][l] * bJK[jGlobal + 1][l] * sumV;
//                }
//            }
//            double xi = 0d;
//            for (int k = 0; k < jhmm.getK(); k++) {
//                for (int l = 0; l < jhmm.getK(); l++) {
//                    double sumV = 0d;
//                    for (int v = 0; v < jhmm.getn(); v++) {
//                        sumV += (read.getBase(j + 1) == v ? jhmm.getAntieps()[jGlobal + 1] : jhmm.getEps()[jGlobal + 1]) * jhmm.getMu()[jGlobal + 1][l][v];
//                    }
//                    double xiSingle = fJK[j][k] * jhmm.getRho()[jGlobal][k][l] * bJK[jGlobal + 1][l] * sumV;
//                    storage.addnJKL(j, k, l, xiSingle);
//                    xi += xiSingle;
//                }
//            }
//        }
        likelihood *= read.getCount();

        jhmm.free(storage.getId());
        return likelihood;
    }
}
