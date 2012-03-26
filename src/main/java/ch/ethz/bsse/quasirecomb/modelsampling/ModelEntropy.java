package ch.ethz.bsse.quasirecomb.modelsampling;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;

/**
 *
 * @author toepfera
 */
public class ModelEntropy {

    public ModelEntropy(String string) {
        OptimalResult or = null;

        try {
            FileInputStream fis = new FileInputStream(string);
            try (ObjectInputStream in = new ObjectInputStream(fis)) {
                or = (OptimalResult) in.readObject();
            }
        } catch (IOException | ClassNotFoundException ex) {
            System.err.println(ex);
        }
        double[][][] mu = or.getMu();
        double sum = 0d;
        for (int j = 0; j < or.getL(); j++) {
            for (int k = 0; k < or.getK(); k++) {
                for (int v = 0; v < or.getn(); v++) {

                    sum -= mu[j][k][v] * Math.log(mu[j][k][v]) / Math.log(4);
                }
            }
        }
        System.out.println(sum / (or.getK()*or.getL()));
    }
}
