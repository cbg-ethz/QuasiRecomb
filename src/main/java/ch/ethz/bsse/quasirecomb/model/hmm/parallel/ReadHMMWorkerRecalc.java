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
package ch.ethz.bsse.quasirecomb.model.hmm.parallel;

import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMM;
import java.util.concurrent.RecursiveTask;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadHMMWorkerRecalc extends RecursiveTask<Boolean> {

    private double[] eps;
    private double[][][] rho;
    private double[] pi;
    private double[][][] mu;
    private int start;
    private int end;
    private ReadHMM[] readHMMArray;

    public ReadHMMWorkerRecalc(ReadHMM[] readHMMArray, double[][][] rho, double[] pi, double[][][] mu, double[] eps, int start, int end) {
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
