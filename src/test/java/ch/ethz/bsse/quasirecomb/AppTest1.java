package ch.ethz.bsse.quasirecomb;

import cc.mallet.types.Dirichlet;
import java.util.Arrays;

/**
 * Unit test for simple App.
 */
public class AppTest1 {

    /**
     * Rigourous Test :-)
     */
    @org.junit.Test
    public void testApp() {
        long a = System.currentTimeMillis();
        Dirichlet d = new Dirichlet(4, 100);
        double[][][] mu = new double[10000][5][4];
        for (int j = 0; j < 10000; j++) {
            for (int k = 0; k < 5; k++) {
                mu[j][k] = d.nextDistribution();
            }
        }
        System.out.println(System.currentTimeMillis()-a);

//        double[] input = new double[]{540, 530, 12, 23};
//        double[] prior = new double[]{0.001, 0.001, 0.001, 0.001};
//        for (int i = 0; i < 100; i++) {
//            input = Regularizations.regularizeOnce(input, i, prior, 1);
//            System.out.println(Arrays.toString(input));
//        }
    }
}
