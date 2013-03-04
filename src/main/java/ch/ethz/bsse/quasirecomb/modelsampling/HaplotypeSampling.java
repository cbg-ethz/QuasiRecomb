/**
 * Copyright (c) 2011-2013 Armin Töpfer
 *
 * This file is part of InDelFixer.
 *
 * InDelFixer is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * InDelFixer is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * InDelFixer. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.quasirecomb.modelsampling;

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class HaplotypeSampling implements Callable<byte[]> {

    private final OptimalResult or;

    public HaplotypeSampling(final OptimalResult or) {
        this.or = or;
    }

    @Override
    public byte[] call() {
        int L = or.getL();
        int n = or.getn();
        int K = or.getK();
        double[][][] rho = or.getRho();
        double[] pi = or.getPi();
        double[][][] mu = or.getMu();
        Frequency<Integer>[][] rhoArray = new Frequency[L - 1][K];
        Frequency<Byte>[][] muArray = new Frequency[L][K];

        Map<Integer, Double> piMap = new HashMap<>();
        for (int k = 0; k < K; k++) {
            piMap.put(k, pi[k]);
        }
        Frequency<Integer> piF = new Frequency<>(piMap);
        int k = piF.roll();

        byte[] read1 = new byte[L];
        for (int j = 0; j < L; j++) {
            if (j > 0) {
                Map<Integer, Double> rhoMap = new HashMap<>();
                for (int l = 0; l < K; l++) {
                    rhoMap.put(l, rho[j - 1][k][l]);
                }
                Frequency<Integer> rhoF = new Frequency<>(rhoMap);
                rhoArray[j - 1][k] = rhoF;
                k = rhoArray[j - 1][k].roll();
            }
            if (muArray[j][k] == null) {
                Map<Byte, Double> muMap = new HashMap<>();
                for (byte v = 0; v < n; v++) {
                    muMap.put(v, mu[j][k][v]);
                }
                Frequency<Byte> muF = new Frequency<>(muMap);
                muArray[j][k] = muF;
            }
            read1[j] = muArray[j][k].roll();
        }
        return read1;
    }
}
