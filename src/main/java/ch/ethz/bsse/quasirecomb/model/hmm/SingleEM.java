package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
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
    private Map<byte[], Integer> reads;
    private byte[][] haplotypesArray;
    private double llh_opt;
    private OptimalResult or;
    private double delta;

    public SingleEM(int N, int K, int L, int n, Map<byte[], Integer> reads, byte[][] haplotypesArray, double delta) {
        this.N = N;
        this.K = K;
        this.L = L;
        this.n = n;
        this.reads = reads;
        this.haplotypesArray = haplotypesArray;
        this.delta = delta;
        start(null);
    }

    public SingleEM(OptimalResult or) {
        this.N = or.getN();
        this.K = or.getK();
        this.L = or.getL();
        this.n = or.getn();
        this.reads = or.getReads();
        this.haplotypesArray = or.getHaplotypesArray();
        this.delta = Globals.DELTA_LLH_HARDER;
        start(or);
    }
    public SingleEM(int N, int K, int L, int n, Map<byte[], Integer> reads, byte[][] haplotypesArray, double delta, OptimalResult or) {
        this.N = N;
        this.K = K;
        this.L = L;
        this.n = n;
        this.reads = reads;
        this.haplotypesArray = haplotypesArray;
        this.delta = delta;
        start(or);
    }

    private void start(OptimalResult givenPrior) {
        this.llh_opt = Globals.getMAX_LLH();
        time(false);
        if (givenPrior == null) {
            jhmm = new JHMM(reads, N, L, K, n, Globals.ESTIMATION_EPSILON);
        } else {
            jhmm = new JHMM(givenPrior);
        }

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
//            if (((oldllh - llh) / llh) < 1e-5 || Globals.NO_BREAK_THRESHOLD) {
//                if (iterations % Globals.MAX_PRE_BREAK == 0 && iterations > 0) {
//                    System.out.println("BIAS CHECK! This: " + llh + "\tbias: " + (llh - (this.llh_opt * Globals.BIAS)) + "\topt:" + this.llh_opt);
//                    if ((llh - (this.llh_opt * Globals.BIAS)) < this.llh_opt) {
//                        System.out.print("pre break;\t");
//                        broken = true;
//                        break;
//                    }
//                }
//            }

            if (iterations > 500) {
//                System.out.print("h: " + (history.get(iterations - 500) - llh) + "\t");
                if (history.get(iterations - 500) - llh > -1) {
                    System.out.print("break 500;\t");
                    broken = true;
                    break;
                }
            }
            log(llh);

            if (Globals.DEBUG) {
                if ((oldllh - llh) / llh == -1) {
                    System.out.print("0\t");
                } else {
                    System.out.print((oldllh - llh) / llh + "\t");
                }
                System.out.println(llh);
            }
            if (Globals.DEBUG && oldllh > llh && oldllh != Double.MIN_VALUE) {
                System.out.println("#" + K + "\t" + iterations);
            }
            jhmm.restart();
            iterations++;
        } while (((oldllh - llh) / llh) > this.delta || oldllh > llh);
        if (!broken) {
            System.out.print("\t\t");
        }

        System.out.println(
                (String.valueOf((oldllh - llh) / llh).contains("-") ? "dist: 1e-" + (String.valueOf((oldllh - llh) / llh).split("-")[1]) : String.valueOf((oldllh - llh) / llh)) + "(" + iterations + ")" + this.llh_opt + "\tthis: " + llh + "\topt:" + this.llh_opt + "\tmax:" + Globals.getMAX_LLH());

        this.calcBic(llh);

        if (Globals.DEBUG) {
            System.out.println("####");
        }
    }

    public void calcBic(double llh) {
        //overview
        double BIC_current = 0;

        // calculate loglikelihood from scaling factors
        for (ReadHMM r : jhmm.getReadHMMMap().keySet()) {
            int times = jhmm.getReadHMMMap().get(r);
            for (int j = 0; j < jhmm.getL(); j++) {
                BIC_current += Math.log(r.getC(j)) * times;
            }
        }

        // count free parameters
        int freeParameters = 0;
        double ERROR = 1e-8;
        
        double[][][] rho = jhmm.getRho();
        for (int i = 0; i < rho.length; i++) {
            for (int j = 0; j < rho[i].length; j++) {
                int f = 0;
                for (int k = 0; k < rho[i][j].length; k++) {
                    if (rho[i][j][k] > ERROR) {
                        f++;
                    }
                }
                freeParameters += f;
            }
        }
        double[][][] mu = jhmm.getMu();
        for (int i = 0; i < mu.length; i++) {
            for (int j = 0; j < mu[i].length; j++) {
                int f = 0;
                for (int k = 0; k < mu[i][j].length; k++) {
                    if (mu[i][j][k] > ERROR) {
                        f++;
                    }
                }
                freeParameters += f;
            }
        }
        // rho
//        freeParameters += L * K * (K - 1);
        // pi
        freeParameters += L - 1;
        // mu
//        freeParameters += L * K * (n - 1);

        BIC_current -= (freeParameters / 2d) * Math.log(N);

        Utils.appendFile(Globals.savePath + "BIC-" + K + ".txt", BIC_current + "\t"+freeParameters+ "\n");

        double[][][] mu_tmp = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                System.arraycopy(jhmm.getMu()[j][k], 0, mu_tmp[j][k], 0, n);
            }
        }

        this.or = new OptimalResult(N, K, L, n, reads, haplotypesArray,
                Arrays.copyOf(jhmm.getRho(), jhmm.getRho().length),
                Arrays.copyOf(jhmm.getPi(),
                jhmm.getPi().length),
                mu_tmp,
                llh,
                BIC_current, jhmm.getPrior_rho(),jhmm.getEps());
        if (llh > llh_opt) {
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
                    System.out.print(iterations + "\t" + t + "\t\t");
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
