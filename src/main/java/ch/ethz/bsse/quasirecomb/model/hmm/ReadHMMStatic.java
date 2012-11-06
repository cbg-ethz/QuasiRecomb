/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;

/**
 *
 * @author toepfera
 */
public class ReadHMMStatic {

    public static double computeFB(JHMM jhmm, Read read) {
        int begin = read.getBegin() - Globals.getINSTANCE().getALIGNMENT_BEGIN();
        int length = read.getLength();

        double[][] fJK = new double[length][jhmm.getK()];
        double[][] bJK = new double[length][jhmm.getK()];
        double[] c = new double[length];
        double[][][] fJKV = new double[length][jhmm.getK()][jhmm.getn()];
        double likelihood = 0;

        for (int j = 0; j < length; j++) {
            int jGlobal = j + begin;
            for (int k = 0; k < jhmm.getK(); k++) {
                for (int v = 0; v < jhmm.getn(); v++) {
                    if (j == 0) {
//                        System.out.println("w_in");
                        fJKV[j][k][v] = jhmm.getTauOmega()[0][j + begin];
                        fJKV[j][k][v] = jhmm.getPi()[k];
                    } else {
                        double sumL = 0d;
                        for (int l = 0; l < jhmm.getK(); l++) {
                            sumL += fJK[j - 1][l] * jhmm.getRho()[begin + j - 1][l][k];
                        }
                        try {
                            switch (read.getPosition(j)) {
                                case WATSON_HIT:
//                                System.out.println("w_hit");
                                    sumL *= 1 - jhmm.getTauOmega()[1][jGlobal];
                                    break;
                                case WATSON_OUT:
//                                System.out.println("w_out");
                                    sumL *= jhmm.getTauOmega()[1][jGlobal];
                                    break;
                                case INSERTION:
//                                System.out.println("insert");
                                    sumL *= 1 - jhmm.getTauOmega()[2][jGlobal];
                                    break;
                                case CRICK_IN:
//                                System.out.println("c_in");
                                    sumL *= jhmm.getTauOmega()[2][jGlobal];
                                    break;
                                case CRICK_HIT:
//                                System.out.println("c_hit");
                                    sumL *= 1 - jhmm.getTauOmega()[3][jGlobal];
                                    break;
                                case CRICK_OUT:
//                                System.out.println("c_out");
                                    break;
                                default:
                                    System.err.println("");
                                    break;
                            }
                        } catch (ArrayIndexOutOfBoundsException e) {
                            System.out.println("w00t");
                        }
                        fJKV[j][k][v] = sumL;
                        if (Double.isNaN(fJKV[j][k][v])) {
                            System.out.println("x");
                        }
                    }

                    if (read.isHit(j)) {
                        fJKV[j][k][v] *= (read.getBase(j) == v ? jhmm.getAntieps()[begin + j] : jhmm.getEps()[begin + j]);
                        fJKV[j][k][v] *= jhmm.getMu()[begin + j][k][v];
                    }
                    if (Double.isNaN(fJKV[j][k][v])) {
                        System.out.println("fJKV:\t" + fJKV[j][k][v]);
                        System.out.println("mu:\t" + jhmm.getMu()[begin + j][k][v]);
                        System.exit(0);
                    }
                    if (k == 0 && v == 0) {
                        c[j] = fJKV[j][k][v];
                    } else {
                        c[j] += fJKV[j][k][v];
                    }
                }
            }
            if (c[j] <= 0) {
                System.out.println("R");
            }
            for (int k = 0; k < jhmm.getK(); k++) {
                fJKV[j][k][0] /= c[j];
                fJK[j][k] = fJKV[j][k][0];
                for (int v = 1; v < jhmm.getn(); v++) {
                    fJKV[j][k][v] /= c[j];
                    fJK[j][k] += fJKV[j][k][v];
                }
            }
        }
        for (int j = length - 1; j >= 0; j--) {
            likelihood += Math.log(c[j]);
            int jGlobal = j + begin;
            for (int k = 0; k < jhmm.getK(); k++) {
                if (j == length - 1) {
                    bJK[j][k] = 1d / c[length - 1];
                } else {
                    bJK[j][k] = 0;
                    for (int l = 0; l < jhmm.getK(); l++) {
                        if (read.isHit(j + 1)) {
                            double sumV = 0d;
                            for (int v = 0; v < jhmm.getn(); v++) {
                                sumV += (read.getBase(j + 1) == v ? jhmm.getAntieps()[jGlobal + 1] : jhmm.getEps()[jGlobal + 1]) * jhmm.getMu()[jGlobal + 1][l][v];
                            }
                            switch (read.getPosition(j)) {
                                case WATSON_HIT:
//                                System.out.println("w_hit");
                                    sumV *= 1 - jhmm.getTauOmega()[1][jGlobal];
                                    break;
                                case WATSON_OUT:
//                                System.out.println("w_out");
                                    sumV *= jhmm.getTauOmega()[1][jGlobal];
                                    break;
                                case INSERTION:
//                                System.out.println("insert");
                                    sumV *= 1 - jhmm.getTauOmega()[2][jGlobal];
                                    break;
                                case CRICK_IN:
//                                System.out.println("c_in");
                                    sumV *= jhmm.getTauOmega()[2][jGlobal];
                                    break;
                                case CRICK_HIT:
//                                System.out.println("c_hit");
                                    sumV *= 1 - jhmm.getTauOmega()[3][jGlobal];
                                    break;
                                case CRICK_OUT:
//                                System.out.println("c_out");
                                    break;
                            }
                            bJK[j][k] += sumV * jhmm.getRho()[begin + j][k][l] * bJK[j + 1][l];
                        } else {
                            bJK[j][k] += jhmm.getRho()[begin + j][k][l] * bJK[j + 1][l];
                        }
                    }
                    if (c[j] != 0d) {
                        bJK[j][k] /= c[j];
                    }
                }
                if (Double.isInfinite(bJK[j][k])) {
                    //this is infinite, because the char has not been observed and there no probability to emit it
                    //thus we divide 0 by a very small number, i.e. 1e-300.
                    bJK[j][k] = 0d;
                }
                if (read.isHit(j)) {
                    for (int v = 0; v < jhmm.getn(); v++) {
                        double gamma = fJKV[j][k][v] * bJK[j][k] * c[j] * read.getCount();
                        jhmm.addnJKV(jGlobal, k, v, gamma);
                        if (Double.isNaN(gamma)) {
                            System.out.println("#####");
                        }
                        if (read.getBase(j) != v) {
                            jhmm.addnneqPos(k, gamma);
                        } else {
                            jhmm.addneqPos(k, gamma);
                        }
                    }
                }
            }
            if (read.isHit(j)) {
                if (j > 0) {
                    for (int k = 0; k < jhmm.getK(); k++) {
                        for (int l = 0; l < jhmm.getK(); l++) {
                            double xi = fJK[j - 1][k] * jhmm.getRho()[begin + j - 1][k][l] * bJK[j][l] * read.getCount();
                            double marginalV = 0d;
                            for (int v = 0; v < jhmm.getn(); v++) {
                                marginalV += (read.getBase(j) == v ? jhmm.getAntieps()[begin + j] : jhmm.getEps()[begin + j]) * jhmm.getMu()[begin + j][l][v];
                            }
                            xi *= marginalV;
                            jhmm.addnJKL(jGlobal, k, l, xi);
                        }
                    }
                }
            }
        }
        likelihood *= read.getCount();
        return likelihood;
    }
}
