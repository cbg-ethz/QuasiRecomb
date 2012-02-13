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

    private static final Dirichlet dir5 = new Dirichlet(5, 2d);
    private static final Dirichlet dir4 = new Dirichlet(4, 2d);
    private static final Dirichlet dir3 = new Dirichlet(3, 2d);
    private static final Dirichlet dir2 = new Dirichlet(2, 2d);
    private static Dirichlet dirPhi;

    public static double[][][] generatePriorRho(int Ldec, int K) {
        double[][][] rho = new double[Ldec][K][K];
        if (!Globals.rho0) {
            if (dirPhi == null) {
                dirPhi = new Dirichlet(Ldec, Globals.PRIOR_ALPHA);
            }
            double[] phi = dirPhi.nextDistribution();
            for (int j = 0; j < Ldec; j++) {
                for (int k = 0; k < K; k++) {
                    double[] d = new double[K];
                    double sum = 0d;
                    for (int l = 0; l < K; l++) {
                        if (l < k) {
//                            d[l] = 0.05d;
//                            d[l] = phi[j];
                            d[l] = phi[j] / (k + 1);
                        } else if (l > k) {
//                            d[l] = 0.05d;
//                            d[l] = phi[j];
                            d[l] = phi[j] / (l + 1);
                        } else {
//                            d[l] = 1d-(K-1)*0.05d;
                            d[l] = ((1 - phi[j]));
                        }
                        if (d[l] == 0d) {
                            d[l] = 1e-200;
                        }
                        sum += d[l];
                    }
                    for (int l = 0; l < K; l++) {
                        d[l] /= sum;
                    }
                    try {
//                        rho[j][k] = new Dirichlet(d).nextDistribution();
                        rho[j][k] = d;
                    } catch (IllegalArgumentException e) {
                        System.err.println(Arrays.toString(d));
                        System.err.println(e);
                        System.exit(0);
                    }
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

    public static double[][][] generateInitRho(int Ldec, int K) {
//        double[][][] rho = new double[Ldec][K][K];
//        if (!Globals.rho0) {
//            for (int j = 0; j < Ldec; j++) {
//                for (int k = 0; k < K; k++) {
//                        rho[j][k] = new Dirichlet(K, 2d).nextDistribution();
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
        return generatePriorRho(Ldec, K);
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
//        if (n == 5) {
//            d = dir5.nextDistribution();
//        } else if (n == 4) {
//            d = dir4.nextDistribution();
//        } else if (n == 3) {
//            d = dir3.nextDistribution();
//        } else if (n == 2) {
//            d = dir2.nextDistribution();
//        }
        d = new double[n];
        for (int i = 0; i < n; i++) {
            d[i] = 1d/n;
        }
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
