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
package ch.ethz.bsse.quasirecomb.distance;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class KullbackLeibler {

    public static double single(double[] p, double[] q) {
        double d = 0d;
        for (int i = 0; i < p.length; i++) {
            if (p[i] > 1e-10 && q[i] > 1e-10) {
                d += p[i] * Math.log(p[i] / q[i]);
            }
        }
        return d;
    }

//    public static double array(double[][] p, double[][] q) {
//        double d = 0d;
//        for (int i = 0; i < p.length; i++) {
//            d += single(p[i], q[i]);
//        }
//        return d;
//    }
    public static double generators(double[][][] mu, int k, int l) {
        double d = 0d;
        for (int i = 0; i < mu.length; i++) {
            d += single(mu[i][k], mu[i][l]);
            if (Double.isInfinite(d)) {
                System.err.println("mu is infinite");
                single(mu[i][k], mu[i][l]);
            }
        }
        return d;
    }

    public static double symmetric(double[][][] mu, int k, int l) {
        double d = generators(mu, k, l) + generators(mu, l, k);
        return d;
    }
    
    public static double shannonEntropy(double[][][] mu, int k) {
        double e = 0d;
        for (int j = 0; j < mu.length; j++) {
            for (int v = 0; v < mu[j][k].length; v++) {
                if (mu[j][k][v] > 0) {
                    e += Math.log(mu[j][k][v]) + mu[j][k][v];
                }
            }
        }
        return -e;
    }
}
