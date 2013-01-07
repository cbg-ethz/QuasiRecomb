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

import ch.ethz.bsse.quasirecomb.distance.KullbackLeibler;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import org.javatuples.Triplet;

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
    private OptimalResult or;
    private double delta;
    private Read[] reads;
    private int repeat;
    private double loglikelihood;
    private int Kmin;
    private double maxBIC;
    private List<Long> times = new ArrayList<>();

    public SingleEM(int N, int K, int L, int n, Read[] reads, double delta, int repeat) {
        this.N = N;
        this.K = K;
        this.Kmin = K;
        this.L = L;
        this.n = n;
        this.delta = delta;
        this.reads = reads;
        this.repeat = repeat;
        time(false);
        if (Globals.getINSTANCE().isPRUNE()) {
            jhmm = new JHMM(reads, N, L, K * 2, n, Globals.getINSTANCE().getESTIMATION_EPSILON(), K);
        } else {
            jhmm = new JHMM(reads, N, L, K, n, Globals.getINSTANCE().getESTIMATION_EPSILON(), K);
        }
        this.K = jhmm.getK();
        start();
    }

    public SingleEM(OptimalResult or, double delta, Read[] reads) {
        this.N = or.getN();
        this.K = or.getK();
        this.Kmin = K;
        this.L = or.getL();
        this.n = or.getn();
        this.delta = delta;
        this.reads = reads;
        this.repeat = -99;
        time(false);
        jhmm = new JHMM(reads, N, L, K, n, or.getEps(), or.getRho(), or.getPi(), or.getMu(), K);
        this.K = jhmm.getK();
        start();
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
        this.loglikelihood = Double.NEGATIVE_INFINITY;
        this.maxBIC = calcBIC(jhmm);

        if (Globals.getINSTANCE().isSNAPSHOTS()) {
            this.snapshot();
        }

        this.iterate();

        Globals.getINSTANCE().log("###c(" + jhmm.getMuChanged() + "|" + jhmm.getRhoChanged() + ")\n");

        Globals.getINSTANCE().incPercentage();

        Globals.getINSTANCE().maxMAX_LLH(loglikelihood);
        this.calcBic();

        if (Globals.getINSTANCE().isDEBUG()) {
            Globals.getINSTANCE().log("####");
            Globals.getINSTANCE().log("\n");
        }
    }

    private void iterate() {
        double oldllh;
        List<Double> history = new LinkedList<>();
        do {
            Globals.getINSTANCE().minMIN_BIC(maxBIC);
            iterations++;
            history.add(loglikelihood);
            oldllh = loglikelihood;
            loglikelihood = jhmm.getLoglikelihood();
            if (Globals.getINSTANCE().isSTOP_QUICK() && Math.abs((oldllh - loglikelihood) / loglikelihood) < 1e-2 && loglikelihood != Globals.getINSTANCE().getMAX_LLH()
                    && ((loglikelihood - Globals.getINSTANCE().getMAX_LLH()) / loglikelihood) > 0.1) {
                if (Globals.getINSTANCE().isDEBUG()) {
                    System.out.println("too small");
                }
                break;
            }
            if (iterations > 500) {
                if (history.get(iterations - 500) - loglikelihood > -1) {
                    Globals.getINSTANCE().log("break 500;\t");
                    break;
                }
            }
            log(loglikelihood);
            Globals.getINSTANCE().setCURRENT_DELTA_LLH((oldllh - loglikelihood) / loglikelihood);
            if (Globals.getINSTANCE().isDEBUG()) {
                if (loglikelihood < 0 && oldllh < 0) {
                    Globals.getINSTANCE().log((oldllh - loglikelihood) / loglikelihood + "\tm(" + jhmm.getMuFlats() + "|" + jhmm.getNjkvFlats() + ")\tr(" + jhmm.getRhoFlats() + "|" + jhmm.getNjklFlats() + ")\tc(" + jhmm.getMuChanged() + "|" + jhmm.getRhoChanged() + ")\t" + ((loglikelihood - Globals.getINSTANCE().getMAX_LLH()) / loglikelihood) + "\t");
                } else if (loglikelihood > 0 && oldllh > 0) {
                    Globals.getINSTANCE().log((loglikelihood - oldllh) / loglikelihood + "\tm(" + jhmm.getMuFlats() + "|" + jhmm.getNjkvFlats() + ")\tr(" + jhmm.getRhoFlats() + "|" + jhmm.getNjklFlats() + ")\tc(" + jhmm.getMuChanged() + "|" + jhmm.getRhoChanged() + ")\t" + ((loglikelihood - Globals.getINSTANCE().getMAX_LLH()) / loglikelihood) + "\t");
                } else if (loglikelihood > 0 && oldllh < 0) {
                    Globals.getINSTANCE().log((loglikelihood + oldllh) / loglikelihood + "\tm(" + jhmm.getMuFlats() + "|" + jhmm.getNjkvFlats() + ")\tr(" + jhmm.getRhoFlats() + "|" + jhmm.getNjklFlats() + ")\tc(" + jhmm.getMuChanged() + "|" + jhmm.getRhoChanged() + ")\t" + ((loglikelihood - Globals.getINSTANCE().getMAX_LLH()) / loglikelihood) + "\t");

                }
                Globals.getINSTANCE().log(loglikelihood + "\n");
            }

            if (Globals.getINSTANCE().isPRUNE() && K > Kmin) {
                this.jhmm = prune();
            } else {
                jhmm.restart();
            }

            if (Globals.getINSTANCE().isSNAPSHOTS()) {
                this.snapshot();
            }
        } while ((Math.abs((oldllh - loglikelihood) / loglikelihood) > this.delta && !Globals.getINSTANCE().isPDELTA()) || (Globals.getINSTANCE().isPDELTA() && (jhmm.getRhoChanged() + jhmm.getMuChanged()) != 0));
    }

    private JHMM prune() {
//        double currentBIC = calcBIC(jhmm);
        double currentBIC = Double.MAX_VALUE;
        Triplet<Integer, Integer, Double> minKL = jhmm.minKL();
        int kPrime = minKL.getValue0();
        double entropyKPrime = KullbackLeibler.shannonEntropy(jhmm.mu, kPrime);
        int lPrime = minKL.getValue1();
        double entropyLPrime = KullbackLeibler.shannonEntropy(jhmm.mu, lPrime);

        int survivingGenerator = entropyKPrime <= entropyLPrime ? kPrime : lPrime;
        int dieingGenerator = entropyKPrime > entropyLPrime ? kPrime : lPrime;

        double[][][] muDeletion = new double[L][K - 1][n];
        double[][][] rhoDeletion = new double[L][K - 1][K - 1];
        double[] piDeletion = new double[K - 1];

        double[][][] muMerging = new double[L][K - 1][n];
        double[][][] rhoMerging = new double[L][K - 1][K - 1];
        double[] piMerging = new double[K - 1];

        int kResorted = 0;
        Map<Integer, Integer> resortedDeletionMap = new HashMap<>();
        for (int k = 0; k < K; k++) {
            if (k != dieingGenerator) {
                resortedDeletionMap.put(k, kResorted);
                kResorted++;
            }
        }

        kResorted = 0;
        Map<Integer, Integer> resortedMergedMap = new HashMap<>();
        for (int k = 0; k < K; k++) {
            if (k == kPrime || k == lPrime) {
                resortedMergedMap.put(k, K - 2);
            } else {
                resortedMergedMap.put(k, kResorted);
                kResorted++;
            }
        }

        //PI DELETION
        double piSumDel = 0d;
        for (int k = 0; k < K; k++) {
            if (resortedDeletionMap.containsKey(k)) {
                int kR = resortedDeletionMap.get(k);
                piDeletion[kR] = jhmm.pi[k];
                piSumDel += jhmm.pi[k];
            }
        }
        for (int k = 0; k < K - 1; k++) {
            piDeletion[k] = piDeletion[k] / piSumDel;
        }

        //PI MERGING
        double piSumMer = 0d;
        for (int k = 0; k < K; k++) {
            if (resortedMergedMap.containsKey(k)) {
                int kR = resortedMergedMap.get(k);
                piMerging[kR] = jhmm.pi[k];
                piSumMer += jhmm.pi[k];
            }
        }
        for (int k = 0; k < K - 1; k++) {
            piMerging[k] = piMerging[k] / piSumMer;
        }


        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                //MU DELETION
                if (resortedDeletionMap.containsKey(k)) {
                    int kRD = resortedDeletionMap.get(k);
                    for (int v = 0; v < n; v++) {
                        muDeletion[j][kRD][v] = jhmm.mu[j][k][v];
                    }
                }
                //MU MERGING
                if (resortedMergedMap.containsKey(k)) {
                    int kRM = resortedMergedMap.get(k);
                    for (int v = 0; v < n; v++) {
                        muMerging[j][kRM][v] += jhmm.mu[j][k][v];
                    }
                }
            }
            for (int v = 0; v < n; v++) {
                muMerging[j][K - 2][v] /= 2;
            }
        }

        //RHO DELETION
        for (int j = 0; j < L - 1; j++) {
            for (int k = 0; k < K; k++) {
                if (resortedDeletionMap.containsKey(k)) {
                    int kR = resortedDeletionMap.get(k);
                    for (int l = 0; l < K; l++) {
                        if (resortedDeletionMap.containsKey(l)) {
                            int lR = resortedDeletionMap.get(l);
                            rhoDeletion[j][kR][lR] = jhmm.rho[j][k][l];
                        }
                    }
                }
            }
        }
        //RHO DEL NORMALIZATION
        for (int j = 0; j < L - 1; j++) {
            for (int k = 0; k < K - 1; k++) {
                double sum = 0;
                for (int l = 0; l < K - 1; l++) {
                    sum += rhoDeletion[j][k][l];
                }
                if (sum > 0) {
                    for (int l = 0; l < K - 1; l++) {
                        rhoDeletion[j][k][l] /= sum;
                    }
                }
            }
        }

        //RHO MERGING
        for (int j = 0; j < L - 1; j++) {
            for (int k = 0; k < K; k++) {
                if (resortedMergedMap.containsKey(k)) {
                    int kR = resortedMergedMap.get(k);
                    for (int l = 0; l < K; l++) {
                        if (resortedMergedMap.containsKey(l)) {
                            int lR = resortedMergedMap.get(l);
                            rhoMerging[j][kR][lR] += jhmm.rho[j][k][l];
                        }
                    }
                }
            }
        }
        //RHO MER NORMALIZATION
        for (int j = 0; j < L - 1; j++) {
            for (int k = 0; k < K - 1; k++) {
                double sum = 0;
                for (int l = 0; l < K - 1; l++) {
                    sum += rhoMerging[j][k][l];
                }
                if (sum > 0) {
                    for (int l = 0; l < K - 1; l++) {
                        rhoMerging[j][k][l] /= sum;
                    }
                }
            }
        }


        JHMM merge = new JHMM(reads, N, L, K - 1, n, Arrays.copyOf(jhmm.getEps(), jhmm.getEps().length), rhoMerging, piMerging, muMerging, Kmin);
        double mergeBIC = calcBIC(merge);
        JHMM del = new JHMM(reads, N, L, K - 1, n, Arrays.copyOf(jhmm.getEps(), jhmm.getEps().length), rhoDeletion, piDeletion, muDeletion, Kmin);
        double delBIC = calcBIC(del);

        maxBIC = currentBIC;
        JHMM argMax = jhmm;
        String s = "C";
        if (delBIC < maxBIC) {
            argMax = del;
            maxBIC = delBIC;
            s = "D";
        }
        if (mergeBIC < maxBIC) {
            argMax = merge;
            maxBIC = mergeBIC;
            s = "M";
        }
        if (s.equals("C")) {
            maxBIC = calcBIC(jhmm);
            argMax.restart();
        }
        this.K = argMax.getK();
        return argMax;
    }

    private double calcBIC(JHMM jhmm) {
        // count free parameters
        double BIC_current = jhmm.getLoglikelihood();
        BIC_current -= (freeParameters(jhmm) / 2d) * Math.log(N);
        return BIC_current;
    }

    private int freeParameters(JHMM jhmm) {
        int freeParameters = 0;
        double ERROR = 1e-15;

        double[][][] rho = jhmm.getRho();
        double[][][] mu = jhmm.getMu();
        double[] pi = jhmm.getPi();
        double[] eps = jhmm.getEps();

        for (int j = 0; j < mu.length; j++) {
            for (int k = 0; k < mu[j].length; k++) {

                //mu
                for (int v = 0; v < mu[j][k].length; v++) {
                    if (mu[j][k][v] > ERROR) {
                        freeParameters++;
                    }
                }

                //rho
                if (j < L - 1) {
                    if (!Globals.getINSTANCE().isNO_RECOMB()) {
                        for (int l = 0; l < rho[j][k].length; l++) {
                            if (rho[j][k][l] > ERROR) {
                                freeParameters++;
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
        return freeParameters;
    }

    private double calcBIC() {
        double BIC_current = this.jhmm.getLoglikelihood();
        BIC_current -= (freeParameters(this.jhmm) / 2d) * Math.log(N);
        return BIC_current;
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

    public void calcBic() {
        //overview
        double BIC_current = this.jhmm.getLoglikelihood();

        BIC_current -= (freeParameters(this.jhmm) / 2d) * Math.log(N);
        if (Globals.getINSTANCE().isLOG_BIC()) {
            Utils.appendFile(Globals.getINSTANCE().getSAVEPATH() + "BIC-" + K + ".txt", BIC_current + "\t" + freeParameters(this.jhmm) + "\n");
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
    }

    public void printMeanTime() {
        long sum = 0;
        for (long l : times) {
            sum += l;
        }
        System.out.println("Mean:" + ((double) sum) / times.size());
    }

    public OptimalResult getOptimalResult() {
        return or;
    }

    public double getLoglikelihood() {
        return loglikelihood;
    }
}
