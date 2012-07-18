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
public class ReadHMMWorkerRecalc extends RecursiveTask<EInfo> {

    private double[] eps;
    private double[] antieps;
    private double[][][] rho;
    private double[][] pi;
    private double[][][] mu;
    private int start;
    private int end;
    private ReadHMM[] readHMMArray;
    private int K;
    private int L;
    private int n;

    public ReadHMMWorkerRecalc(int K, int L, int n, ReadHMM[] readHMMArray, double[][][] rho, double[][] pi, double[][][] mu, double[] eps, double[] antieps, int start, int end) {
        this.readHMMArray = readHMMArray;
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.eps = eps;
        this.antieps = antieps;
        this.start = start;
        this.end = end;
        this.K = K;
        this.L = L;
        this.n = n;
    }

    @Override
    protected EInfo compute() {
        if (end - start < Globals.STEPSIZE) {
            EInfo einfo = new EInfo(K, L, n);
            for (int i = start; i < end; i++) {
                ReadHMM r = this.readHMMArray[i];
                r.recalc(rho, pi, mu, eps, antieps);
                int offset = r.getBegin();
                int times = r.getCount();
                //CONTINUE
                for (int j = 0; j < r.getLength(); j++) {
                    int jGlobal = offset + j;
                    einfo.coverage[jGlobal] += times;
                    for (int k = 0; k < K; k++) {
                        einfo.nJK[jGlobal][k] += r.gamma(j, k) * times;
                        if (j > 0) {
                            for (int l = 0; l < K; l++) {
                                einfo.nJKL[jGlobal][k][l] += r.xi(j, k, l) * times;
                                if (k == l) {
                                    einfo.nJeq[jGlobal] += r.xi(j, k, l) * times;
                                } else {
                                    einfo.nJneq[jGlobal] += r.xi(j, k, l) * times;
                                }
                            }
                        }
                        for (int v = 0; v < n; v++) {
                            einfo.nJKV[jGlobal][k][v] += r.gamma(j, k, v) * times;
                        }
                    }
                    for (int v = 0; v < n; v++) {
                        byte b = r.getSequence()[j];
                        for (int k = 0; k < K; k++) {
                            if (v != b) {
                                einfo.nneqPos[jGlobal] += r.gamma(j, k, v) * times;
                            } else {
                                einfo.neqPos[jGlobal] += r.gamma(j, k, v) * times;
                            }
                        }
                    }
                }
            }
            return einfo;
        } else {
            final int mid = start + (end - start) / 2;
            ReadHMMWorkerRecalc left = new ReadHMMWorkerRecalc(K, L, n, readHMMArray, rho, pi, mu, eps, antieps, start, mid);
            ReadHMMWorkerRecalc right = new ReadHMMWorkerRecalc(K, L, n, readHMMArray, rho, pi, mu, eps, antieps, mid, end);
            left.fork();
            EInfo einfo = right.compute();
            einfo.add(left.join());
            return einfo;
        }
    }
}
