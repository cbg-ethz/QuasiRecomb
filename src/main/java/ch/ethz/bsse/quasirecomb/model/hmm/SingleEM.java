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

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
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
    private int repeat;
    private double loglikelihood;

    public SingleEM(int N, int K, int L, int n, Read[] reads, double delta, int repeat) {
        this.N = N;
        this.K = K;
        this.L = L;
        this.n = n;
        this.delta = delta;
        this.reads = reads;
        this.repeat = repeat;
        time(false);
        jhmm = new JHMM(reads, N, L, K, n, Globals.getINSTANCE().getESTIMATION_EPSILON());
        start();
    }

    public SingleEM(OptimalResult or, double delta, Read[] reads) {
        this.N = or.getN();
        this.K = or.getK();
        this.L = or.getL();
        this.n = or.getn();
        this.delta = delta;
        this.reads = reads;
        this.repeat = -99;
        time(false);
        jhmm = new JHMM(reads, N, L, K, n, or.getEps(), or.getRho(), or.getPi(), or.getMu());
        start();
    }

    public String getOptimumPath() {
        if (!Globals.getINSTANCE().isSNAPSHOTS()) {
            this.snapshot();
        }

        String save = Globals.getINSTANCE().getSnapshotDir() + (Globals.getINSTANCE().isMODELSELECTION() ? "modelselection" : "training") + File.separator + "R" + (repeat < 10 ? "00" : repeat < 100 ? "0" : "") + repeat + "_K" + K + "_" + (iterations < 10 ? "000" : iterations < 100 ? "00" : iterations < 1000 ? "0" : "") + iterations;
        OptimalResult localOr = new OptimalResult(N, K, L, n,
                jhmm.getRho(),
                jhmm.getPi(),
                jhmm.getMu(),
                this.jhmm.getLoglikelihood(),
                calcBIC(), jhmm.getEps(), jhmm.getRestart(), jhmm.getTauOmega());
        Utils.saveOptimum(save + ".optimum", localOr);

        return Globals.getINSTANCE().getSnapshotDir() + (Globals.getINSTANCE().isMODELSELECTION() ? "modelselection" : "training") + File.separator + "R" + (repeat < 10 ? "00" : repeat < 100 ? "0" : "") + repeat + "_K" + K + "_" + (iterations < 10 ? "000" : iterations < 100 ? "00" : iterations < 1000 ? "0" : "") + iterations + ".optimum";
    }

    private void snapshot() {
        String save = Globals.getINSTANCE().getSnapshotDir() + (Globals.getINSTANCE().isMODELSELECTION() ? "modelselection" : "training") + File.separator + "R" + (repeat < 10 ? "00" : repeat < 100 ? "0" : "") + repeat + "_K" + K + "_" + (iterations < 10 ? "000" : iterations < 100 ? "00" : iterations < 1000 ? "0" : "") + iterations;
        OptimalResult localOr = new OptimalResult(N, K, L, n,
                jhmm.getRho(),
                jhmm.getPi(),
                jhmm.getMu(),
                this.jhmm.getLoglikelihood(),
                calcBIC(), jhmm.getEps(), jhmm.getRestart(), jhmm.getTauOmega());
//        Utils.saveOptimum(save + ".optimum", localOr);
        Summary summary = new Summary();
        Utils.saveFile(save + ".txt", summary.print(localOr));
//        Utils.saveFile(save + ".html", summary.html(localOr));
//        System.out.println("saved: " + (System.currentTimeMillis() - time));
    }

    private void start() {
        this.loglikelihood = Double.MIN_VALUE;
        if (Globals.getINSTANCE().isSNAPSHOTS()) {
            this.snapshot();
        }

        double oldllh;
        List<Double> history = new ArrayList<>();
        do {
//            System.out.println("FLATS:"+flats);
            if (iterations % 10 == 0 && iterations > 0) {
                Globals.getINSTANCE().maxMAX_LLH(loglikelihood);
                this.llh_opt = Math.max(Globals.getINSTANCE().getMAX_LLH(), this.llh_opt);
            }
            iterations++;
            history.add(loglikelihood);
            oldllh = loglikelihood;
            loglikelihood = jhmm.getLoglikelihood();
            if (iterations > 500) {
                if (history.get(iterations - 500) - loglikelihood > -1) {
                    Globals.getINSTANCE().log("break 500;\t");
                    break;
                }
            }
            log(loglikelihood);

            if (Globals.getINSTANCE().isDEBUG()) {
                if (loglikelihood < 0 && oldllh < 0) {
                    Globals.getINSTANCE().log((oldllh - loglikelihood) / loglikelihood + "\tm(" + jhmm.getMuFlats() + "|" + jhmm.getNjkvFlats() + ")\tr(" + jhmm.getRhoFlats() + "|" + jhmm.getNjklFlats() + ")\t" + jhmm.getParametersChanged() + "\t");
                } else if (loglikelihood > 0 && oldllh > 0) {
                    Globals.getINSTANCE().log((loglikelihood - oldllh) / loglikelihood + "\tm(" + jhmm.getMuFlats() + "|" + jhmm.getNjkvFlats() + ")\tr(" + jhmm.getRhoFlats() + "|" + jhmm.getNjklFlats() + ")\t" + jhmm.getParametersChanged() + "\t");
                } else if (loglikelihood > 0 && oldllh < 0) {
                    Globals.getINSTANCE().log((loglikelihood + oldllh) / loglikelihood + "\tm(" + jhmm.getMuFlats() + "|" + jhmm.getNjkvFlats() + ")\tr(" + jhmm.getRhoFlats() + "|" + jhmm.getNjklFlats() + ")\t" + jhmm.getParametersChanged() + "\t");
                    
                } 
//                Globals.getINSTANCE().log((oldllh - loglikelihood) / loglikelihood + jhmm.getParametersChanged() + "\t");
                Globals.getINSTANCE().log(loglikelihood + "\n");
            }
            jhmm.restart();
            Globals.getINSTANCE().printPercentage(K);
            if (Globals.getINSTANCE().isSNAPSHOTS()) {
                this.snapshot();
            }
        } while ((Math.abs((oldllh - loglikelihood) / loglikelihood) > this.delta && !Globals.getINSTANCE().isPDELTA()) || (Globals.getINSTANCE().isPDELTA() && jhmm.getParametersChanged() != 0));

//        } while (Math.abs((oldllh - loglikelihood) / loglikelihood) > this.delta);
//        } while (iterations <= 500);
        Globals.getINSTANCE().log("###\t" + jhmm.getParametersChanged() + "\n");

        Globals.getINSTANCE().incPercentage();
        this.calcBic();

        if (Globals.getINSTANCE().isDEBUG()) {
            Globals.getINSTANCE().log("####");
        }
    }
    
    private double calcBIC() {
        double BIC_current = this.jhmm.getLoglikelihood();
        // count free parameters
        int freeParameters = 0;
        double ERROR = 1e-15;

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
        return BIC_current;
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
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                System.arraycopy(jhmm.getMu()[j][k], 0, muCopy[j][k], 0, n);
            }
        }
        double[][] tauOmegaCopy = new double[jhmm.getTauOmega().length][L + 1];
        for (int j = 0; j < jhmm.getTauOmega().length; j++) {
            System.arraycopy(jhmm.getTauOmega()[j], 0, tauOmegaCopy[j], 0, L + 1);
        }
        this.or = new OptimalResult(N, K, L, n,
                Arrays.copyOf(jhmm.getRho(), jhmm.getRho().length),
                Arrays.copyOf(jhmm.getPi(), jhmm.getPi().length),
                muCopy,
                this.jhmm.getLoglikelihood(),
                BIC_current, Arrays.copyOf(jhmm.getEps(), jhmm.getEps().length), jhmm.getRestart(), tauOmegaCopy);

        if (this.jhmm.getLoglikelihood()
                >= llh_opt) {
            Globals.getINSTANCE().maxMAX_LLH(this.jhmm.getLoglikelihood());
        }
    }
    private List<Long> times = new ArrayList<>();

    public void printMeanTime() {
        long sum = 0;
        for (long l : times) {
            sum += l;
        }
        System.out.println("Mean:" + ((double) sum) / times.size());
    }

    private long time(boolean show) {
        long t = 0;
        if (time == -1) {
            time = System.currentTimeMillis();
        } else {
            t = (System.currentTimeMillis() - time);
            times.add(t);
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

    public double getLoglikelihood() {
        return loglikelihood;
    }
}
