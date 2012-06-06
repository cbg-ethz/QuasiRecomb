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

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.parallel.RestartWorker;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

/**
 * Responsible for the start of the repeats within the EM algorithm.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class EM extends Utils {

    private StringBuilder sb = new StringBuilder();
    private double BIC_opt = Double.NEGATIVE_INFINITY;
    private int iterations = 0;
    private int N;
    private int K;
    private int L;
    private int n;
    private Map<byte[], Integer> reads;
    private byte[][] haplotypesArray;
    private OptimalResult or;
    private List<OptimalResult> ors;

    protected EM(int N, int L, int K, int n, String[] reads, byte[][] haplotypesArray) {
        this(K, N, L, n, Utils.clusterReads(Utils.splitReadsIntoByteArrays(reads)), haplotypesArray);
    }

    protected EM(int N, int L, int K, int n, Map<byte[], Integer> reads, byte[][] haplotypesArray) {
        this.N = N;
        this.L = L;
        this.K = K;
        this.n = n;
        this.reads = reads;
        this.haplotypesArray = haplotypesArray;
        this.blackbox();
    }

    protected EM(int N, int L, int K, int n, byte[][] reads, byte[][] haplotypesArray) {
        this(K, N, L, n, Utils.clusterReads(reads), haplotypesArray);
    }

    public void shorten(double value) {
        if (value < 1e-5) {
            sb.append("0      ");
        } else if (value == 1.0) {
            sb.append("1      ");
        } else {
            String t = "" + value;
            String r;
            if (t.length() > 5) {
                r = t.substring(0, 7);
                if (t.contains("E")) {
                    r = r.substring(0, 4);
                    r += "E" + t.split("E")[1];
                }
                sb.append(r);
            } else {
                sb.append(value);
            }
        }
    }

    private void blackbox() {
        Globals.LOG = new StringBuilder();

        if (Globals.PARALLEL_RESTARTS) {
            ors = Globals.fjPool.invoke(new RestartWorker(N, K, L, n, reads, haplotypesArray, Globals.DELTA_LLH, 0, Globals.REPEATS));
        } else {
            ors = new ArrayList<>();
            for (int i = 0; i < Globals.REPEATS; i++) {
                SingleEM sem = new SingleEM(N, K, L, n, reads, haplotypesArray, Globals.DELTA_LLH);
                ors.add(sem.getOptimalResult());
            }
        }

        Globals.PARALLEL_JHMM = true;
        double maxLLH = Double.NEGATIVE_INFINITY;
        StringBuilder restarts = new StringBuilder();
        for (OptimalResult tmp : ors) {
            restarts.append(tmp.getRestarts()).append("\n");
            if (tmp != null && tmp.getLlh() >= maxLLH) {
                maxLLH = tmp.getLlh();
                or = tmp;
            }
        }

        System.out.println("\tBIC: " + (int) or.getBIC());
//        if (!Globals.NO_REFINE) {
//            SingleEM bestEM = new SingleEM(N, K, L, n, reads, haplotypesArray, 1e-10, or);
//            this.or = bestEM.getOptimalResult();
//        }
//        this.saveEM();
        Globals.log("\n" + new Summary().print(or));
        if (Globals.LOGGING) {
            Utils.saveFile(Globals.savePath + "log_K" + K, Globals.LOG.toString());
            Utils.saveFile(Globals.savePath + "restarts_K" + K, restarts.toString());
        }
    }

    /**
     * Returns the StringBuilder with the information of the best result.
     *
     * @return StringBuilder with result
     */
    public StringBuilder getSb() {
        return sb;
    }

    /**
     * Returns the best BIC.
     *
     * @return best BIC
     */
    public double getBIC_opt() {
        return BIC_opt;
    }

    /**
     * The double 3D array of mu for the best result.
     *
     * @return mu
     */
    public double[][][] getMu_opt() {
        return or.getMu();
    }

    /**
     * The double array of pi for the best result.
     *
     * @return pi
     */
    public double[] getPi_opt() {
        return or.getPi();
    }

    public List<OptimalResult> getOrs() {
        return ors;
    }

    /**
     * The double 3D array of rho for the best result.
     *
     * @return rho
     */
    public double[][][] getRho_opt() {
        return or.getRho();
    }

    public OptimalResult getOr() {
        return or;
    }
}
