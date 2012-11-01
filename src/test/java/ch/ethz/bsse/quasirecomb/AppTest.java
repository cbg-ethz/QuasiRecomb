package ch.ethz.bsse.quasirecomb;

import cc.mallet.types.Dirichlet;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apfloat.Apfloat;

/**
 * Unit test for simple App.
 */
public class AppTest
        extends TestCase {

    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest(String testName) {
        super(testName);
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite() {
        return new TestSuite(AppTest.class);
    }

    private double[][][] calcMu(double[][][] nJKV) {
        int n = 4;
        int L = 1;
        int K = 1;

        double[][][] mu = new double[L][K][n];
        double[] muJKV = new double[n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {

                double AH = 1e-10;
                double divisor;
                double sum;

                double[] orig = nJKV[j][k];
                double preSum = 0d;
                for (int v = 0; v < n; v++) {
                    preSum = orig[v];
                }
                if (preSum == 0d) {
                    for (int v = 0; v < n; v++) {
                        mu[j][k][v] = 1d / n;
                    }
                    continue;
                }
                if (preSum < 1e-5) {
                    System.out.println("scale 1");
                    double min = Double.MAX_VALUE;
                    double max = Double.MAX_VALUE;
                    for (int v = 0; v < n; v++) {
                        if (nJKV[j][k][v] > 0) {
                            min = Math.min(min, nJKV[j][k][v]);
                            max = Math.max(max, nJKV[j][k][v]);
                        }
                    }
                    for (int v = 0; v < n; v++) {
                        nJKV[j][k][v] /= min;
                    }
                }

                do {
                    sum = 0d;
                    divisor = 0d;
                    for (int v = 0; v < n; v++) {
                        muJKV[v] = f(nJKV[j][k][v] + AH);
                        sum += nJKV[j][k][v];
                    }
                    sum = f(sum + n * AH);
                    if (sum > 0) {
                        for (int v = 0; v < n; v++) {
                            muJKV[v] /= sum;
                            divisor += muJKV[v];
                        }
                        if (divisor > 0) {
                            for (int v = 0; v < n; v++) {
                                muJKV[v] /= divisor;
                            }
                        } else {
                            AH *= 10;
                            System.out.println("scale up");
                        }
                    } else {
                        System.out.println("scale up 2");
                        AH *= 10;
                    }
                } while (sum == 0 || divisor == 0);

                for (int v = 0; v < n; v++) {
                    mu[j][k][v] = muJKV[v];
                }
            }
        }
        return mu;
    }

    private double f(double upsilon) {
        if (upsilon == 0d) {
            return 0d;
        }
        return Math.exp(rho_phi(upsilon));
    }

    private double rho_phi(double upsilon) {
        return Dirichlet.digamma(upsilon);
//        double x = -1d;
//        try {
//            x = (upsilon > 7) ? rho_g(upsilon - .5) : (rho_phi(upsilon + 1) - 1 / upsilon);
//
//        } catch (StackOverflowError e) {
//            System.err.println(upsilon);
//            System.err.println(e);
//        }
//        return x;
    }

    private double rho_g(double x) {
        return Math.log(x) + .04167 * Math.pow(x, -2) - .00729 * Math.pow(x, -4) + .00384 * Math.pow(x, -6) - .00413 * Math.pow(x, -8);
    }

    /**
     * Rigourous Test :-)
     */
    public void testApp() throws MathException {
        long sum1 = 0l;
        long sum2 = 0l;
        long t1;
        long t2;

        double[][] matrixData = new double[10000][5];
        double[][] matrixData2 = new double[10000][5];
        for (int i = 0; i < 10000; i++) {
            for (int j = 0; j < 5; j++) {
                matrixData[i][j] = Math.random();
                matrixData2[i][j] = Math.random();
            }
        }
        System.out.println("a");
        RealMatrix m = new Array2DRowRealMatrix(matrixData);
        RealMatrix n = new Array2DRowRealMatrix(matrixData2);
        for (int x = 0; x < 1000; x++) {
            t1 = System.nanoTime();

            RealMatrix add = m.add(n);
            double[][] data = add.getData();

            sum1 += System.nanoTime() - t1;
            System.out.println("b");
            t2 = System.nanoTime();

            double[][] matrixData3 = new double[10000][5];
            for (int i = 0; i < 10000; i++) {
                for (int j = 0; j < 5; j++) {
                    matrixData3[i][j] = matrixData2[i][j] + matrixData3[i][j];
                }
            }

            sum2 += System.nanoTime() - t2;
            System.out.println("c");
        }
        System.out.println(sum1);
        System.out.println(sum2);
        System.out.println((double) sum1 / sum2);




        //        System.out.println(      Globals.getINSTANCE().time());
        //        for (int i = 0; i < 10; i++) {
        //            long time = System.currentTimeMillis();
        //            for (int j = 0; j < 1000000; j++) {
        //                double rho_g = rho_phi(Math.random());
        //            }
        //            time = System.currentTimeMillis() - time;
        //            long time2 = System.currentTimeMillis();
        //            for (int j = 0; j < 1000000; j++) {
        //                double digamma = Dirichlet.digamma(Math.random());
        //            }
        //            time2 = System.currentTimeMillis() - time2;
        //            System.out.println(time-time2);
        //        }
        //                GammaDistributionImpl gamma = new GammaDistributionImpl(1,2);
        //        double[] sample = gamma.sample(100);
        //        System.out.println("");
        //        double[][][] a = new double[][][]{{{1 * 1e-100, 1 * 1e-100, 1 * 1e-100, 3 * 1e-98}}};
        //        System.out.println(Arrays.toString(a[0][0]));
        //        double[][][] mu = calcMu(a);
        //        for (;;) {
        //            System.out.println(Arrays.toString(mu[0][0]));
        //            mu = calcMu(mu.clone());
        //        }
        //        Read[] reads = FastaParser.parseFastaPairedEnd("C:/Users/XLR/Dropbox/simulationStudy/a1.fasta");
        //        String genome = FastaParser.parseFarFile("C:/Users/XLR/Dropbox/simulationStudy/haplotypes/haploytpes_1_1.fasta")[0];
        //        for (Read r : reads) {
        //            if (genome.contains(Utils.reverse(r.getSequence()))) {
        //                System.out.println("a");
        //            } else {
        //                System.out.println("'''''''");
        //            }
        //        }
        //        double eArray[][][] = new double[5000][500][500];
        //        long minus = 0;
        //        long plus = 0;
        //        for (int x = 0; x < 10; x++) {
        //
        //            long time = System.currentTimeMillis();
        //            for (int j = 4999; j >= 0; j--) {
        //                for (int k = 499; k >= 0; k--) {
        //                    for (int i = 499; i >= 0; i--) {
        //                        eArray[j][k][i] = 1d / 5;
        //                    }
        //                }
        //            }
        //            minus += (System.currentTimeMillis() - time);
        //            time = System.currentTimeMillis();
        //            for (int j = 0; j < 5000; j++) {
        //                for (int k = 0; k < 500; k++) {
        //                    for (int i = 0; i < 500; i++) {
        //                        eArray[j][k][i] = 1d / 5;
        //                    }
        //                }
        //            }
        //            plus += (System.currentTimeMillis() - time);
        //        }
        //        System.out.println(minus/10d);
        //        System.out.println(plus/10d);
    }
}
