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

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class EInfo {

    public double[][] nJK;
    public double[][][] nJKL;
    public double[][][] nJKV;
    public double[] neqPos;
    public double[] nneqPos;
    public double[] nJeq;
    public double[] nJneq;
    public int[] coverage;
    private int L;
    private int K;
    private int n;

    public EInfo(int K, int L, int n) {
        this.nJK = new double[L][K];
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nJeq = new double[L];
        this.nJneq = new double[L];
        this.neqPos = new double[L];
        this.nneqPos = new double[L];
        this.coverage = new int[L];
        this.L = L;
        this.K = K;
        this.n = n;
    }
    
    public void add(EInfo e) {
        for (int j = 0; j < L; j++) {
            this.nJeq[j] += e.nJeq[j];
            this.nJneq[j] += e.nJneq[j];
            this.neqPos[j] += e.neqPos[j];
            this.nneqPos[j] += e.nneqPos[j];
            this.coverage[j] += e.coverage[j];
            for (int k = 0; k < K; k++) {
                this.nJK[j][k] += e.nJK[j][k];
                for (int l = 0; l < K; l++) {
                    this.nJKL[j][k][l] += e.nJKL[j][k][l];
                }
                for (int v = 0; v < n; v++) {
                    this.nJKV[j][k][v] += e.nJKV[j][k][v];
                }
            }
        }
    }
}
