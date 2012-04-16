package ch.ethz.bsse.quasirecomb.model.hmm.parallel;

import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMM;
import java.util.concurrent.RecursiveTask;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class ReadHMMWorkerRecalc extends RecursiveTask<Boolean> {

    private double[][] eps;
    private double[][][] rho;
    private double[] pi;
    private double[][][] mu;
    private int start;
    private int end;
    private ReadHMM[] readHMMArray;

    public ReadHMMWorkerRecalc(ReadHMM[] readHMMArray, double[][][] rho, double[] pi, double[][][] mu, double[][] eps, int start, int end) {
        this.readHMMArray = readHMMArray;
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.eps = eps;
        this.start = start;
        this.end = end;
    }

    @Override
    protected Boolean compute() {
        if (end - start < Globals.STEPSIZE) {
            long time = System.currentTimeMillis();
            for (int i = start; i < end; i++) {
                this.readHMMArray[i].recalc(rho, pi, mu, eps);
            }
//            System.out.println(start + "-" + end + "\t:" + (System.currentTimeMillis() - time));

            return true;
        } else {
            final int mid = start + (end - start) / 2;
            ReadHMMWorkerRecalc left = new ReadHMMWorkerRecalc(readHMMArray, rho, pi, mu, eps, start, mid);
            ReadHMMWorkerRecalc right = new ReadHMMWorkerRecalc(readHMMArray, rho, pi, mu, eps, mid, end);
            left.fork();
            right.compute();
            left.join();
            return true;
        }
    }
}
