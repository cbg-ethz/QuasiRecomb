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
package ch.ethz.bsse.quasirecomb.utils;

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Summary extends Utils {

    public String print(OptimalResult or) {
        StringBuilder sb = new StringBuilder();
        sb.setLength(0);
        sb.append("#loglikelihood:").append(or.getLlh()).append("\n");
        sb.append("#BIC:").append(or.getBIC()).append("\n");
        double[][][] mu = or.getMu();
//        double mue = 0d;
//        for (int j = 0; j < or.getL(); j++) {
//            for (int k = 0; k < or.getK(); k++) {
//                for (int v = 0; v < or.getn(); v++) {
//                    mue -= mu[j][k][v] * Math.log(mu[j][k][v]) / Math.log(or.getn());
//                }
//            }
//        }
//        sb.append("#MUE:").append(mue / (or.getK() * or.getL())).append("\n");
//        double rhoe = 0d;
        double[][][] rho = or.getRho();
//        for (int j = 0; j < or.getL() - 1; j++) {
//            for (int k = 0; k < or.getK(); k++) {
//                for (int l = 0; l < or.getK(); l++) {
//                    rhoe -= rho[j][k][l] * Math.log(rho[j][k][l]) / Math.log(or.getK());
//                }
//            }
//        }
//        sb.append("#RHOE:").append(mue / (or.getK() * or.getL())).append("\n");
//        sb.append("#PE:").append((mue - rhoe) / (or.getK() * or.getL())).append("\n");
        sb.append("\n#EPS:").append("\n");
        for (int j = 0; j < or.getL(); j++) {
            sb.append("##j:").append(j).append("\t").append(shorten(or.getEps()[j])).append("\n");
        }
        sb.append("\n");

        sb.append("#PI:\n");
        for (int k = 0; k < or.getK(); k++) {
            sb.append("##").append(or.getPi()[0][k]).append("\n");
        }
        sb.append("\n");
        sb.append("#RHO:\n");
        for (int j = 0; j < or.getL() - 1; j++) {
            boolean flat = false;
            for (int k = 0; k < or.getK(); k++) {
                for (int l = 0; l < or.getK(); l++) {
                    if (or.getRho()[j][k][l] > 1e-20 && or.getRho()[j][k][l] != 1d) {
                        flat = true;
                        break;
                    }
                }
                if (flat) {
                    break;
                }
            }
            if (flat) {
                sb.append(j + 1).append("\t");
                for (int k = 0; k < or.getK(); k++) {
                    sb.append("[");
                    for (int l = 0; l < or.getK(); l++) {
                        sb.append(shorten(or.getRho()[j][k][l]));
                        if (l + 1 < or.getK()) {
                            sb.append(", ");
                        }

                    }
                    sb.append("]\t");
                }
                sb.append("\n");
            }
        }
        sb.append("\n");
        sb.append("#MU:\n");
        for (int j = 0; j < or.getL(); j++) {
            boolean flat = false;
            for (int k = 0; k < or.getK(); k++) {
                for (int v = 0; v < or.getMu()[0][0].length; v++) {
                    if (or.getMu()[j][k][v] > 1e-20 && or.getMu()[j][k][v] != 1d) {
                        flat = true;
                        break;
                    }
                }
                if (flat) {
                    break;
                }
            }
            if (flat) {
                sb.append("##j:").append(j).append("\t");
//                for (byte[] b : or.getHaplotypesArray()) {
//                    sb.append(reverse((int) b[j]));
//                }
//                sb.append("|");
                for (int k = 0; k < or.getK(); k++) {
                    double max = Double.MIN_VALUE;
                    Map<Double, Integer> m = new HashMap<>();
                    for (int v = 0; v < or.getMu()[0][0].length; v++) {
                        max = Math.max(max, or.getMu()[j][k][v]);
                        m.put(or.getMu()[j][k][v], v);
                    }
                    sb.append(reverse(m.get(max))).append("-");
                }
                    double entropy = 0d;
                for (int k = 0; k < or.getK(); k++) {
                    sb.append("[");
                    for (int v = 0; v < or.getn(); v++) {
                        if (or.getMu()[j][k][v] > 0)
                        entropy -= or.getMu()[j][k][v]*Math.log(or.getMu()[j][k][v]) ;
                        sb.append(shorten(or.getMu()[j][k][v]));
                        if (v + 1 < or.getn()) {
                            sb.append(", ");
                        }
                    }
                    sb.append("]\t");
                }
//                sb.append(entropy+"\n");
                sb.append("\n");
            }
        }
        return sb.toString();
    }

    private String shorten(double value) {
        String s;
        if (value < 1e-20) {
            s = "0      ";
        } else if (value == 1.0) {
            s = "1      ";
        } else {
            String t = "" + value;
            String r;
            if (t.length() > 7) {
                r = t.substring(0, 7);
                if (t.contains("E")) {
                    r = r.substring(0, 4);
                    r += "E" + t.split("E")[1];
                }
                s = r;
            } else {
                s = String.valueOf(value);
            }
        }
        return s;
    }
}
