/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.distance;

/**
 *
 * @author XLR
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
                System.err.println("");
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
