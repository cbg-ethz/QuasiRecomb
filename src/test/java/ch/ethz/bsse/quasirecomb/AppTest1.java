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
        Dirichlet dirichlet = new Dirichlet(4, 100000);
        for (int i = 0; i < 100; i++) {
            System.out.println(Arrays.toString(dirichlet.nextDistribution()));
        }

//        double[] input = new double[]{540, 530, 12, 23};
//        double[] prior = new double[]{0.001, 0.001, 0.001, 0.001};
//        for (int i = 0; i < 100; i++) {
//            input = Regularizations.regularizeOnce(input, i, prior, 1);
//            System.out.println(Arrays.toString(input));
//        }
    }
}
