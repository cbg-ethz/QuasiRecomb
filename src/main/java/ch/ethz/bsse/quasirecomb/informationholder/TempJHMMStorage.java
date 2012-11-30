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
package ch.ethz.bsse.quasirecomb.informationholder;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class TempJHMMStorage {

    private double[][][] nJKL;
    private double[][][] nJKV;
    private double[] nneqPos;
    private int id;

    public TempJHMMStorage(int L, int K, int n, int id) {
        this.id = id;
        this.nJKL = new double[L][K][K];
        this.nJKV = new double[L][K][n];
        this.nneqPos = new double[L];
    }

    public void addnJKV(int j, int k, int v, double value) {
        this.nJKV[j][k][v] += value;
    }

    public void addnJKL(int j, int k, int l, double value) {
        this.nJKL[j][k][l] += value;
    }

    public void addnneqPos(int j, double value) {
        this.nneqPos[j] += value;
    }

    public int getId() {
        return id;
    }

    public TempJHMMStorage merge(TempJHMMStorage t) {
        for (int j = 0; j < nJKL.length; j++) {
            for (int k = 0; k < nJKL[j].length; k++) {
                for (int l = 0; l < nJKL[j][k].length; l++) {
                    this.nJKL[j][k][l] += t.nJKL[j][k][l];
                }
                for (int v = 0; v < nJKV[j][k].length; v++) {
                    this.nJKV[j][k][v] += t.nJKV[j][k][v];
                }
            }
            this.nneqPos[j] += t.nneqPos[j];
        }
        return this;
    }

    public double[][][] getnJKL() {
        return nJKL;
    }

    public double[][][] getnJKV() {
        return nJKV;
    }

    public double[] getNneqPos() {
        return nneqPos;
    }
}
