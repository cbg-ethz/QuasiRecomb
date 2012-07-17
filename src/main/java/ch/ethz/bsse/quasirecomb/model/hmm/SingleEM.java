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

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SingleEM {

    private long time = -1;
    private StringBuilder sb = new StringBuilder();
    private JHMM jhmm;
    private int iterations = 0;
    private int N;
    private int K;
    private int L;
    private int n;
    private double llh_opt;
    private OptimalResult or;
    private double delta;
    private Read[] reads;

    public SingleEM(int N, int K, int L, int n, Read[] reads, double delta) {
        this.N = N;
        this.K = K;
        this.L = L;
        this.n = n;
        this.delta = delta;
        this.reads = reads;
        start();
    }

    private void start() {
        this.llh_opt = Globals.getMAX_LLH();
        time(false);
        jhmm = new JHMM(reads, N, L, K, n, Globals.ESTIMATION_EPSILON);

        double llh = Double.MIN_VALUE;
        double oldllh = Double.MIN_VALUE;
        List<Double> history = new ArrayList<>();
        boolean broken = false;
        do {
            if (iterations % 10 == 0 && iterations > 0) {
                Globals.maxMAX_LLH(llh);
                this.llh_opt = Math.max(Globals.getMAX_LLH(), this.llh_opt);
            }
            history.add(llh);
            oldllh = llh;
            llh = jhmm.getLoglikelihood();
            if (iterations > 500) {
                if (history.get(iterations - 500) - llh > -1) {
                    Globals.log("break 500;\t");
                    broken = true;
                    break;
                }
            }
            log(llh);

            if (Globals.DEBUG) {
                if ((oldllh - llh) / llh == -1) {
                    Globals.log("0\t");
                } else {
                    Globals.log((oldllh - llh) / llh + "\t");
                }
                Globals.log(llh + "\n");
            }
            jhmm.restart();
            iterations++;

        } while (Math.abs((oldllh - llh) / llh) > this.delta);
        if (!broken) {
            Globals.log("\t\t");
        }

        Globals.log((String.valueOf((oldllh - llh) / llh).contains("-") ? "dist: 1e-" + (String.valueOf((oldllh - llh) / llh).split("-")[1]) : String.valueOf((oldllh - llh) / llh)) + "(" + iterations + ")" + this.llh_opt + "\tthis: " + llh + "\topt:" + this.llh_opt + "\tmax:" + Globals.getMAX_LLH());

        this.calcBic(llh);

        if (Globals.DEBUG) {
            Globals.log("####");
        }
        Globals.printPercentage(K);
    }

    public void calcBic(double llh) {
        //overview
        double BIC_current = 0;

        // calculate loglikelihood from scaling factors
        for (ReadHMM r : jhmm.getReadHMMArray()) {
            int times = r.getCount();
            for (int j = 0; j < jhmm.getL(); j++) {
                BIC_current += Math.log(r.getC(j)) * times;
            }
        }

        // count free parameters
        int freeParameters = 0;
        double ERROR = 1e-8;

        double[][][] rho = jhmm.getRho();
        double[][][] mu = jhmm.getMu();
        double[] pi = jhmm.getPi();
        double[] eps = jhmm.getEps();
        for (int j = 0; j < rho.length; j++) {
            for (int k = 0; k < rho[j].length; k++) {
                boolean different = false;
                for (int l = 1; l < rho[j][k].length; l++) {
                    if (Math.abs(rho[j][k][l - 1] - rho[j][k][l]) > ERROR) {
                        different = true;
                        break;
                    }
                }
                if (different) {
                    for (int l = 0; l < rho[j][k].length; l++) {
                        if (rho[j][k][l] > ERROR) {
                            freeParameters++;
                        }
                    }
                }
            }
        }
        for (int j = 0; j < mu.length; j++) {
            for (int k = 0; k < mu[j].length; k++) {
                boolean different = false;
                for (int v = 1; v < mu[j][k].length; v++) {
                    if (Math.abs(mu[j][k][v - 1] - mu[j][k][v]) > ERROR) {
                        different = true;
                        break;
                    }
                }
                if (different) {
                    for (int v = 0; v < mu[j][k].length; v++) {
                        if (mu[j][k][v] > ERROR) {
                            freeParameters++;
                        }
                    }
                }
            }
        }
        for (int k = 0; k < pi.length; k++) {
            if (pi[k] > ERROR) {
                freeParameters++;
            }
        }
        for (int j = 0; j < eps.length; j++) {
            if (eps[j] > ERROR) {
                freeParameters++;
            }
        }

        BIC_current -= (freeParameters / 2d) * Math.log(N);

        if (Globals.LOG_BIC) {
            Utils.appendFile(Globals.savePath + "BIC-" + K + ".txt", BIC_current + "\t" + freeParameters + "\n");
        }

        double[][][] mu_tmp = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                System.arraycopy(jhmm.getMu()[j][k], 0, mu_tmp[j][k], 0, n);
            }
        }
        this.or = new OptimalResult(N, K, L, n,
                Arrays.copyOf(jhmm.getRho(), jhmm.getRho().length),
                Arrays.copyOf(jhmm.getPi(), jhmm.getPi().length),
                mu_tmp,
                llh,
                BIC_current, jhmm.getEps(), jhmm.getRestart());
        if (llh >= llh_opt) {
            Globals.maxMAX_LLH(llh);
        }
    }

    private long time(boolean show) {
        long t = 0;
        if (time == -1) {
            time = System.currentTimeMillis();
        } else {
            t = (System.currentTimeMillis() - time);
            if (show) {
                sb.append(t).append("\t\t");
                if (Globals.DEBUG) {
                    Globals.log(iterations + "\t" + t + "\t\t");
                }
            }
            time = System.currentTimeMillis();
        }
        return t;
    }

    private void log(double llh) {
        Globals.runtime.add((int) time(true));
        sb.append(llh).append("\t\t").append("\n");
    }

    public OptimalResult getOptimalResult() {
        return or;
    }

    public double getLlh_opt() {
        return llh_opt;
    }
}
