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
package ch.ethz.bsse.quasirecomb.model.hmm;

import cc.mallet.types.Dirichlet;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Regularizations {

    public static double[] regularizeOnce(double[] estCounts, int restart, double[] hyperParameter, double fidelity) {
//        double hyperParameter = 0.001;
        int x = estCounts.length;
        double[] regCounts = new double[x];
        double divisor;

        double sum = 0d;
        double max = Double.MIN_VALUE;
        for (int v = 0; v < x; v++) {
            sum += estCounts[v];
            max = Math.max(estCounts[v], max);
        }
        if (sum == 0) {
            for (int v = 0; v < x; v++) {
                regCounts[v] = 1d / x;
            }
            return regCounts;
        }
        if (Math.abs(max - sum) < 1e-8) {
            for (int v = 0; v < x; v++) {
                if (estCounts[v] < max) {
                    regCounts[v] = 0d;
                } else {
                    regCounts[v] = 1;
                }
            }
            return regCounts;
        }
        for (int v = 0; v < x; v++) {
            regCounts[v] = fidelity * estCounts[v] / sum;
        }

        sum = 0d;
        double hyperSum = 0d;
        divisor = 0d;
        for (int i = 0; i < x; i++) {
            regCounts[i] = f(regCounts[i] + hyperParameter[i]);
            sum += regCounts[i];
            hyperSum = hyperParameter[i];
        }
        sum = f(sum + hyperSum);
        if (sum > 0) {
            for (int i = 0; i < x; i++) {
                regCounts[i] /= sum;
                divisor += regCounts[i];
            }
            if (divisor > 0) {
                for (int i = 0; i < x; i++) {
                    regCounts[i] /= divisor;
                }
            }
        }
        if (Double.isNaN(sum) || Double.isNaN(divisor)) {
            System.out.println("reg nan");
            System.exit(0);
        }
        max = Double.MIN_VALUE;
        for (int v = 0; v < x; v++) {
            max = Math.max(max, regCounts[v]);
        }
        if (Math.abs(max - 1d) < 1e-8) {
            for (int v = 0; v < x; v++) {
                if (regCounts[v] < max) {
                    regCounts[v] = 0d;
                } else {
                    regCounts[v] = 1;
                }
            }
        } else {
//            for (int v = 0; v < x; v++) {
//                if (regCounts[v] <= 1e-7) {
//                    regCounts[v] = 0d;
//                }
//            }
            if (restart < Globals.getINSTANCE().getPERTURB()) {
                sum = 0;
                for (int l = 0; l < x; l++) {
                    regCounts[l] += Math.random() / (10 * (restart + 1));
                    sum += regCounts[l];
                }
                for (int l = 0; l < x; l++) {
                    regCounts[l] /= sum;
                }
            }
        }
        return regCounts;
    }

    public static double f(double upsilon) {
        if (upsilon == 0d) {
            return 0d;
        }
        return Math.exp(Dirichlet.digamma(upsilon));
    }
}