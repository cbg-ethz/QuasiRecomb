package ch.ethz.bsse.quasirecomb.utils;

import cc.mallet.types.Dirichlet;
import ch.ethz.bsse.quasirecomb.model.Globals;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class Random {

    private static final Dirichlet dir5 = new Dirichlet(5, Globals.BETA_H);
    private static final Dirichlet dir4 = new Dirichlet(4, Globals.BETA_H);
    private static final Dirichlet dir3 = new Dirichlet(3, Globals.BETA_H);
    private static final Dirichlet dir2 = new Dirichlet(2, Globals.BETA_H);
    private static Dirichlet dirPhi;

//    public static double[][][] generatePriorRho(int Ldec, int K) {
//        double[][][] rho = new double[Ldec][K][K];
//        if (!Globals.rho0) {
//            if (dirPhi == null) {
//                dirPhi = new Dirichlet(K, Globals.BETA_Z);
//            }
//            double[] phi = dirPhi.nextDistribution();
//            for (int j = 0; j < Ldec; j++) {
//                for (int k = 0; k < K; k++) {
//                    rho[j][k] = dirPhi.nextDistribution();
//                }
//            }
//        } else {
//            for (int j = 0; j < Ldec; j++) {
//                for (int k = 0; k < K; k++) {
//                    for (int l = 0; l < K; l++) {
//                        rho[j][k][l] = (l == k) ? 1d : 0d;
//                    }
//                }
//            }
//        }
//        return rho;
//    }

    public static double[][][] generateInitRho(int Ldec, int K) {
        double[][][] rho = new double[Ldec][K][K];
        if (!Globals.rho0) {
            for (int j = 0; j < Ldec; j++) {
                for (int k = 0; k < K; k++) {
                    rho[j][k] = new Dirichlet(K, Globals.BETA_Z).nextDistribution();
                    int maxIndex = 0;
                    double max = 0;
                    for (int l = 0; l < K; l++) {
                        if (rho[j][k][l] > max) {
                            max = rho[j][k][l];
                            maxIndex = l;
                        }
                    }
                    rho[j][k][maxIndex] = rho[j][k][k];
                    rho[j][k][k] = max;
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
    private static Map<Integer, Dirichlet> piGen = new HashMap<>();

    public static double[] generateInitPi(int K) {
        if (!piGen.containsKey(K)) {
            double[] d = new double[K];
            for (int k = 1; k <= K; k++) {
//            d[k - 1] = 1d / K;
                d[k - 1] = 2d;
            }
            piGen.put(K, new Dirichlet(d));
        }
        return piGen.get(K).nextDistribution();
//        return d;
    }

    public static double[] generateMuVector(int n) {
        double[] d = null;
        if (n == 5) {
            d = dir5.nextDistribution();
        } else if (n == 4) {
            d = dir4.nextDistribution();
        } else if (n == 3) {
            d = dir3.nextDistribution();
        } else if (n == 2) {
            d = dir2.nextDistribution();
        }
//        d = new double[n];
//        for (int i = 0; i < n; i++) {
//            d[i] = 1d / n;
//        }
        return d;
    }

    public static double[][][] generateMuInit(int L, int K, int n) {
        double[][][] eArray = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                eArray[j][k] = Random.generateMuVector(n);
            }
        }
        return eArray;
    }
}
