/**
 * Copyright (c) 2011-2013 Armin Töpfer
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
import ch.ethz.bsse.quasirecomb.informationholder.TauOmega;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import com.google.common.collect.Lists;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SingleModelSampling implements Callable<List<SampledRead>> {

    private final OptimalResult or;
    private final TauOmega tauOmega;
    private final boolean paired;
    private final int amount;
    private final Frequency<Integer> startF;

    public SingleModelSampling(final OptimalResult or, final TauOmega tauOmega, final boolean paired, int a, int b, Frequency<Integer> startF) {
        this.tauOmega = tauOmega;
        this.or = or;
        this.paired = paired;
        this.amount = b - a;
        this.startF = startF;
    }

    @Override
    public List<SampledRead> call() {
        List<SampledRead> list = Lists.newArrayListWithCapacity(this.amount);
        for (int i = 0; i < amount; i++) {
            int L = or.getL();
            int n = or.getn();
            int K = or.getK();
            double[][][] rho = or.getRho();
            double[] pi = or.getPi()[0];
            double[][][] mu = or.getMu();
            Frequency<Integer>[][] rhoArray = new Frequency[L - 1][K];
            Frequency<Byte>[][] muArray = new Frequency[L][K];
            int watsonStart = 0, watsonLength = 0, watsonEnd = 0, crickStart = 0, crickEnd = 0;

            watsonStart = startF.roll();

            Map<Integer, Double> watsonLengthMap = new HashMap<>();
            for (Map.Entry<Integer, Double> e : this.tauOmega.getTauWatsonMap().get(watsonStart).entrySet()) {
                watsonLengthMap.put(e.getKey(), e.getValue());
            }
            Frequency<Integer> watsonLengthF = new Frequency<>(watsonLengthMap);
            watsonLength = watsonLengthF.roll();
            watsonEnd = watsonStart + watsonLength;
            boolean localPaired = paired && this.tauOmega.getOmegaWatsonMap().get(watsonEnd) != null;
            if (localPaired) {
                Map<Integer, Double> watsonOmegaMap = new HashMap<>();
                for (Map.Entry<Integer, Double> e : this.tauOmega.getOmegaWatsonMap().get(watsonEnd).entrySet()) {
                    watsonOmegaMap.put(e.getKey(), e.getValue());
                }
                Frequency<Integer> inLength = new Frequency<>(watsonOmegaMap);
                int insertLength = inLength.roll();
                crickStart = watsonEnd + insertLength;


                Map<Integer, Double> crickTauMap = new HashMap<>();
                if (this.tauOmega.getTauCrickMap().get(crickStart) == null) {
                    System.err.println("");
                }
                for (Map.Entry<Integer, Double> e : this.tauOmega.getTauCrickMap().get(crickStart).entrySet()) {
                    crickTauMap.put(e.getKey(), e.getValue());
                }
                Frequency<Integer> crickL = new Frequency<>(crickTauMap);
                int crickLength = crickL.roll();
                crickEnd = crickStart + crickLength;
            }
            Map<Integer, Double> piMap = new HashMap<>();
            for (int k = 0; k < K; k++) {
                piMap.put(k, pi[k]);
            }
            Frequency<Integer> piF = new Frequency<>(piMap);
            int k = piF.roll();
            List<Byte> read1 = new LinkedList<>();
            for (int j = watsonStart; j < watsonEnd; j++) {
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
                read1.add(muArray[j][k].roll());
            }
            if (localPaired) {
                List<Byte> read2 = new LinkedList<>();
                for (int j = watsonEnd; j < crickStart; j++) {
                    Map<Integer, Double> rhoMap = new HashMap<>();
                    for (int l = 0; l < K; l++) {
                        rhoMap.put(l, rho[j - 1][k][l]);
                    }
                    Frequency<Integer> rhoF = new Frequency<>(rhoMap);
                    rhoArray[j - 1][k] = rhoF;
                    k = rhoArray[j - 1][k].roll();
                }
                for (int j = crickStart; j < crickEnd; j++) {
                    Map<Integer, Double> rhoMap = new HashMap<>();
                    for (int l = 0; l < K; l++) {
                        rhoMap.put(l, rho[j - 1][k][l]);
                    }
                    Frequency<Integer> rhoF = new Frequency<>(rhoMap);
                    rhoArray[j - 1][k] = rhoF;
                    k = rhoArray[j - 1][k].roll();
                    if (muArray[j][k] == null) {
                        Map<Byte, Double> muMap = new HashMap<>();
                        for (byte v = 0; v < n; v++) {
                            muMap.put(v, mu[j][k][v]);
                        }
                        Frequency<Byte> muF = new Frequency<>(muMap);
                        muArray[j][k] = muF;
                    }
                    read2.add(muArray[j][k].roll());
                }
                list.add(new SampledRead(read1, read2, watsonStart, watsonEnd, crickStart, crickEnd));
            } else {
                list.add(new SampledRead(read1, watsonStart, watsonEnd));
            }
        }
        return list;
    }
}
