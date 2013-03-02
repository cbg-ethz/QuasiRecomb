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
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SingleModelSampling implements Callable<SampledRead> {

    private final OptimalResult or;
    private final double[][] tauOmega;
    private final boolean paired;

    public SingleModelSampling(final OptimalResult or, final double[][] tauOmega, final boolean paired) {
        this.tauOmega = tauOmega;
        this.or = or;
        this.paired = paired;
    }

    @Override
    public SampledRead call() {
        int L = or.getL();
        int n = or.getn();
        int K = or.getK();
        double[][][] rho = or.getRho();
        double[] pi = or.getPi();
        double[][][] mu = or.getMu();
        Frequency<Integer>[][] rhoArray = new Frequency[L - 1][K];
        Frequency<Byte>[][] muArray = new Frequency[L][K];
        int start = 0, stop = 0, start2 = 0, stop2 = 0;
        while (true) {
            try {
//                System.out.print(".");
                Map<Integer, Double> startMap = new HashMap<>();
                for (int j = 0; j < L + 1; j++) {
                    startMap.put(j, tauOmega[0][j]);
                }
                Frequency<Integer> startF = new Frequency<>(startMap);
                start = startF.roll();

                Map<Integer, Double> stopMap = new HashMap<>();
                for (int j = 0; j < L + 1; j++) {
                    stopMap.put(j, tauOmega[1][j]);
                }
                Frequency<Integer> stopF = new Frequency<>(stopMap);
                stop = stopF.roll();
                if (paired) {
                    Map<Integer, Double> start2Map = new HashMap<>();
                    for (int j = 0; j < L + 1; j++) {
                        start2Map.put(j, tauOmega[2][j]);
                    }
                    Frequency<Integer> start2F = new Frequency<>(start2Map);
                    start2 = start2F.roll();

                    Map<Integer, Double> stop2Map = new HashMap<>();
                    for (int j = 0; j < L + 1; j++) {
                        stop2Map.put(j, tauOmega[3][j]);
                    }
                    Frequency<Integer> stop2F = new Frequency<>(stop2Map);
                    stop2 = stop2F.roll();
                    if (start < stop && stop < start2 && start2 < stop2) {
                        break;
                    } else {
                        continue;
                    }
                } else if (start < stop) {
                    break;
                } else {
                    continue;
                }
            } catch (Exception e) {
                continue;
            }
        }
//        System.out.println("");
//        StringBuilder startStopSB = new StringBuilder();
//        startStopSB.append(start).append("\t");
//        startStopSB.append(stop).append("\t");
//        startStopSB.append(start2).append("\t");
//        startStopSB.append(stop2).append("\n");

//        for (int j = start; j < stop; j++) {
//            this.coverage[j]++;
//        }
//        for (int j = start2; j < stop2; j++) {
//            this.coverage[j]++;
//        }

        Map<Integer, Double> piMap = new HashMap<>();
        for (int k = 0; k < K; k++) {
            piMap.put(k, pi[k]);
        }
        Frequency<Integer> piF = new Frequency<>(piMap);
        int k = piF.roll();

        List<Byte> read1 = new LinkedList<>();
        for (int j = start; j < stop; j++) {
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
        if (paired) {
            List<Byte> read2 = new LinkedList<>();
            for (int j = stop; j < start2; j++) {
                if (j > 0) {
                    Map<Integer, Double> rhoMap = new HashMap<>();
                    for (int l = 0; l < K; l++) {
                        rhoMap.put(l, rho[j - 1][k][l]);
                    }
                    Frequency<Integer> rhoF = new Frequency<>(rhoMap);
                    rhoArray[j - 1][k] = rhoF;
                    k = rhoArray[j - 1][k].roll();
                }
            }
            for (int j = start2; j < stop2; j++) {
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
                read2.add(muArray[j][k].roll());
            }
            return new SampledRead(read1, read2, start, stop, start2, stop2);
        } else {
            return new SampledRead(read1, start, stop);
        }
    }
}
