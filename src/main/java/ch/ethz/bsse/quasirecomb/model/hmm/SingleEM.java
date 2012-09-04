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
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
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
        this.llh_opt = Globals.getINSTANCE().getMAX_LLH();
        time(false);
        jhmm = new JHMM(reads, N, L, K, n, Globals.getINSTANCE().getESTIMATION_EPSILON());

        double llh = Double.MIN_VALUE;
        double oldllh;
        List<Double> history = new ArrayList<>();
        do {
            if (iterations % 10 == 0 && iterations > 0) {
                Globals.getINSTANCE().maxMAX_LLH(llh);
                this.llh_opt = Math.max(Globals.getINSTANCE().getMAX_LLH(), this.llh_opt);
            }
            history.add(llh);
            oldllh = llh;
            llh = jhmm.getLoglikelihood();
            if (iterations > 500) {
                if (history.get(iterations - 500) - llh > -1) {
                    Globals.getINSTANCE().log("break 500;\t");
                    break;
                }
            }
            log(llh);

            if (Globals.getINSTANCE().isDEBUG()) {
                if (Math.abs((oldllh - llh) / llh - 1d) > 1e-15) {
                    Globals.getINSTANCE().log("0\t");
                } else {
                    Globals.getINSTANCE().log((oldllh - llh) / llh + "\t" + jhmm.getParametersChanged() + "\t");
                }
                Globals.getINSTANCE().log(llh + "\n");
            }
            jhmm.restart();
            iterations++;
            Globals.getINSTANCE().printPercentage(K);
        } while (Math.abs((oldllh - llh) / llh) > this.delta && jhmm.getParametersChanged() != 0);
        Globals.getINSTANCE().log("###\t" + jhmm.getParametersChanged() + "\n");
        Utils.appendFile(Globals.getINSTANCE().getSAVEPATH() + "p.txt", jhmm.getParametersChanged() + "\n");

        Globals.getINSTANCE().incPercentage();
        this.calcBic();

        if (Globals.getINSTANCE().isDEBUG()) {
            Globals.getINSTANCE().log("####");
        }
    }

    public void calcBic() {
        //overview
        double BIC_current = this.jhmm.getLoglikelihood();

        // count free parameters
        int freeParameters = 0;
        double ERROR = 1e-8;

        double[][][] rho = jhmm.getRho();
        double[][][] mu = jhmm.getMu();
        double[] pi = jhmm.getPi();
        double[] eps = jhmm.getEps();

        for (int j = 0; j < mu.length; j++) {
            for (int k = 0; k < mu[j].length; k++) {

                //mu
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

                //rho
                if (j < L - 1) {
                    if (!Globals.getINSTANCE().isNO_RECOMB()) {
                        different = false;
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
                if (eps[j] > ERROR) {
                    freeParameters++;
                }
            }

            for (int k = 0; k < pi.length; k++) {
                if (pi[k] > ERROR) {
                    freeParameters++;
                }
            }
        }
        BIC_current -= (freeParameters / 2d) * Math.log(N);
        if (Globals.getINSTANCE().isLOG_BIC()) {
            Utils.appendFile(Globals.getINSTANCE().getSAVEPATH() + "BIC-" + K + ".txt", BIC_current + "\t" + freeParameters + "\n");
        }
        double[][][] muCopy = new double[L][K][n];
        for (int j = 0;
                j < L;
                j++) {
            for (int k = 0; k < K; k++) {
                System.arraycopy(jhmm.getMu()[j][k], 0, muCopy[j][k], 0, n);
            }
        }
        this.or = new OptimalResult(N, K, L, n,
                Arrays.copyOf(jhmm.getRho(), jhmm.getRho().length),
                Arrays.copyOf(jhmm.getPi(), jhmm.getPi().length),
                muCopy,
                this.jhmm.getLoglikelihood(),
                BIC_current, Arrays.copyOf(jhmm.getEps(), jhmm.getEps().length), jhmm.getRestart());

        if (this.jhmm.getLoglikelihood()
                >= llh_opt) {
            Globals.getINSTANCE().maxMAX_LLH(this.jhmm.getLoglikelihood());
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
                if (Globals.getINSTANCE().isDEBUG()) {
                    Globals.getINSTANCE().log(iterations + "\t" + t + "\t\t");
                }
            }
            time = System.currentTimeMillis();
        }
        return t;
    }

    private void log(double llh) {
        Globals.getINSTANCE().getRuntime().add((int) time(true));
        sb.append(llh).append("\t\t").append("\n");
    }

    public OptimalResult getOptimalResult() {
        return or;
    }

    public double getLlh_opt() {
        return llh_opt;
    }
}
