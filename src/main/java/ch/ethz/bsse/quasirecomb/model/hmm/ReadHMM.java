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

/**
 * Calculates the forward / backward algorithm for a single read.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadHMM {

    private final Read read;
    private double[][] fJK;
    private double[][] bJK;
    private double[] c;
    private double[][][] fJKV;
    private int begin;
    private int length;
    private JHMM jhmm;

    public ReadHMM(JHMM jhmm, Read read) {
        this.read = read;
        this.begin = read.getBegin() - Globals.getINSTANCE().getALIGNMENT_BEGIN();
        this.length = this.read.getLength();
        this.jhmm = jhmm;
        this.fJKV = new double[length][jhmm.getK()][jhmm.getn()];
        this.fJK = new double[length][jhmm.getK()];
        this.bJK = new double[length][jhmm.getK()];
        this.c = new double[length];
        calculate();
    }

    private void calculate() {
        this.forward();
        this.backward();
    }

    public boolean checkConsistency() {
        for (int j = 0; j < jhmm.getL(); j++) {
            for (int k = 0; k < jhmm.getK(); k++) {
                if (Double.isNaN(this.fJK[j][k])) {
                    return true;
                }
                for (int v = 0; v < jhmm.getn(); v++) {
                    if (Double.isNaN(this.fJKV[j][k][v])) {
                        return true;
                    }
                }
            }
            if (Double.isNaN(this.c[j])) {
                return true;
            }
        }
        for (int j = 0; j < jhmm.getL() - 1; j++) {
            for (int k = 0; k < jhmm.getK(); k++) {
                if (Double.isNaN(this.bJK[j][k])) {
                    return true;
                }
            }
        }
        return true;
    }

    public void recalc() {
        this.calculate();
    }

    private void forward() {
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
                            switch (this.read.getPosition(j)) {
                                case WATSON_HIT:
//                                System.out.println("w_hit");
                                    sumL *= 1 - this.jhmm.getTauOmega()[1][jGlobal];
                                    break;
                                case WATSON_OUT:
//                                System.out.println("w_out");
                                    sumL *= this.jhmm.getTauOmega()[1][jGlobal];
                                    break;
                                case INSERTION:
//                                System.out.println("insert");
                                    sumL *= 1 - this.jhmm.getTauOmega()[2][jGlobal];
                                    break;
                                case CRICK_IN:
//                                System.out.println("c_in");
                                    sumL *= this.jhmm.getTauOmega()[2][jGlobal];
                                    break;
                                case CRICK_HIT:
//                                System.out.println("c_hit");
                                    sumL *= 1 - this.jhmm.getTauOmega()[3][jGlobal];
                                    break;
                                case CRICK_OUT:
//                                System.out.println("c_out");
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

                    if (this.read.isHit(j)) {
                        fJKV[j][k][v] *= (this.getSequence(j) == v ? jhmm.getAntieps()[begin + j] : jhmm.getEps()[begin + j]);
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
    }

    private void backward() {
        for (int j = length - 1; j >= 0; j--) {
            int jGlobal = j + begin;
            for (int k = 0; k < jhmm.getK(); k++) {
                if (j == length - 1) {
                    bJK[j][k] = 1d / c[length - 1];
                } else {
                    bJK[j][k] = 0;
                    for (int l = 0; l < jhmm.getK(); l++) {
                        if (this.read.isHit(j + 1)) {
                            double sumV = 0d;
                            for (int v = 0; v < jhmm.getn(); v++) {
                                sumV += (this.getSequence(j + 1) == v ? jhmm.getAntieps()[jGlobal + 1] : jhmm.getEps()[jGlobal + 1]) * jhmm.getMu()[jGlobal + 1][l][v];
                            }
                            switch (this.read.getPosition(j)) {
                                case WATSON_HIT:
//                                System.out.println("w_hit");
                                    sumV *= 1 - this.jhmm.getTauOmega()[1][jGlobal];
                                    break;
                                case WATSON_OUT:
//                                System.out.println("w_out");
                                    sumV *= this.jhmm.getTauOmega()[1][jGlobal];
                                    break;
                                case INSERTION:
//                                System.out.println("insert");
                                    sumV *= 1 - this.jhmm.getTauOmega()[2][jGlobal];
                                    break;
                                case CRICK_IN:
//                                System.out.println("c_in");
                                    sumV *= this.jhmm.getTauOmega()[2][jGlobal];
                                    break;
                                case CRICK_HIT:
//                                System.out.println("c_hit");
                                    sumV *= 1 - this.jhmm.getTauOmega()[3][jGlobal];
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
                if (this.read.isHit(j)) {
                    for (int v = 0; v < jhmm.getn(); v++) {
                        double gamma = this.fJKV[j][k][v] * this.bJK[j][k] * c[j] * this.getCount();
                        this.jhmm.addnJKV(jGlobal, k, v, gamma);
                        if (Double.isNaN(gamma)) {
                            System.out.println("#####");
                        }
                        if (this.read.getSequence()[j] != v) {
                            this.jhmm.addnneqPos(k, gamma);
                        } else {
                            this.jhmm.addneqPos(k, gamma);
                        }
                    }
                }
            }
            if (this.read.isHit(j)) {
                if (j > 0) {
                    for (int k = 0; k < jhmm.getK(); k++) {
                        for (int l = 0; l < jhmm.getK(); l++) {
                            double xi = this.fJK[j - 1][k] * jhmm.getRho()[begin + j - 1][k][l] * this.bJK[j][l] * this.getCount();
                            double marginalV = 0d;
                            for (int v = 0; v < jhmm.getn(); v++) {
                                marginalV += (this.getSequence(j) == v ? jhmm.getAntieps()[begin + j] : jhmm.getEps()[begin + j]) * jhmm.getMu()[begin + j][l][v];
                            }
                            xi *= marginalV;
                            this.jhmm.addnJKL(jGlobal, k, l, xi);
                            if (k == l) {
                                this.jhmm.addnJeq(jGlobal, xi);
                            } else {
                                this.jhmm.addnJneq(jGlobal, xi);
                            }
                        }
                    }
                }
            }
        }
    }

    final public double getC(int j) {
        if (j >= begin && j - begin < length) {
            return c[j - begin];
        }
        return 1;
    }

    final public byte getSequence(int j) {
        return this.read.getBase(j);
    }

    public int getCount() {
        return this.read.getCount();
    }

    public int getLength() {
        return length;
    }

    public int getBegin() {
        return begin;
    }

    public Read getRead() {
        return read;
    }
}
