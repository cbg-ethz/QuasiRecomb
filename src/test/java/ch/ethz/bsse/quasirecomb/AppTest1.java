package ch.ethz.bsse.quasirecomb;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.Regularizations;
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
        double[] input = new double[]{540, 530, 12, 23};
        double[] prior = new double[]{0.001, 0.001, 0.001, 0.001};
        for (int i = 0; i < 100; i++) {
            input = Regularizations.regularizeOnce(input, i, prior, 1);
            System.out.println(Arrays.toString(input));
        }
    }
}
