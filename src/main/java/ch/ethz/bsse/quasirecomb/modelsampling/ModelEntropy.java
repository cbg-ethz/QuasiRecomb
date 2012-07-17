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
package ch.ethz.bsse.quasirecomb.modelsampling;

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
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
