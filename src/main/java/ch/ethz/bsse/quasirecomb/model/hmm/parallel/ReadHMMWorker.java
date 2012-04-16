package ch.ethz.bsse.quasirecomb.model.hmm.parallel;

import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.JHMM;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMM;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.RecursiveTask;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class ReadHMMWorker extends RecursiveTask<List<ReadHMM>> {

    private JHMM jhmm;
    private byte[][] reads;
    private double[][][] rho;
    private double[] pi;
    private double[][][] mu;
    private double[][] eps;
    private double[][] antieps;
    private int K;
    private int L;
    private int n;
    private int start;
    private int end;

    public ReadHMMWorker(JHMM jhmm, byte[][] reads, double[][][] rho, double[] pi, double[][][] mu, double[][] eps, double[][] antieps, int K, int L, int n, int start, int end) {
        this.jhmm = jhmm;
        this.reads = reads;
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.eps = eps;
        this.antieps = antieps;
        this.K = K;
        this.L = L;
        this.n = n;
        this.start = start;
        this.end = end;
    }

    @Override
    protected List<ReadHMM> compute() {
        if (end - start < Globals.STEPSIZE) {
            List<ReadHMM> list = new LinkedList<>();
//            long time = System.currentTimeMillis();
            for (int i = start; i < end; i++) {
                ReadHMM readHMM = new ReadHMM(L, K, n, reads[i], rho, pi, mu, eps, antieps);
                list.add(readHMM);
            }
//            System.out.println(start+"-"+end+"\t:"+(System.currentTimeMillis()-time));
            return list;
        } else {
            final int mid = start + (end - start) / 2;
            ReadHMMWorker left = new ReadHMMWorker(jhmm, reads, rho, pi, mu, eps, antieps, K, L, n, start, mid);
            ReadHMMWorker right = new ReadHMMWorker(jhmm, reads, rho, pi, mu, eps, antieps, K, L, n, mid, end);
            left.fork();
            List<ReadHMM> list = new LinkedList<>();
            list.addAll(right.compute());
            list.addAll(left.join());
            return list;
        }
    }
}
