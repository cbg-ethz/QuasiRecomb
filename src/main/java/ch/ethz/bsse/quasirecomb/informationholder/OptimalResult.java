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

import java.io.Serializable;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class OptimalResult implements Serializable {

    private int N;
    private int K;
    private int L;
    private int n;
    private Read[] reads;
    private double[][][] rho;
    private double[][][] priorRho;
    private double[] pi;
    private double[][][] mu;
    private double llh;
    private double[] eps;
    private double BIC;
    private int restarts;

    public OptimalResult(int N, int K, int L, int n, Read[] reads,  double[][][] rho, double[] pi, double[][][] mu, double llh, double BIC, double[][][] priorRho, double[] eps, int restarts) {
        this.N = N;
        this.K = K;
        this.L = L;
        this.n = n;
        this.reads = reads;
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.llh = llh;
        this.BIC = BIC;
        this.priorRho = priorRho;
        this.eps = eps;
        this.restarts = restarts;
    }

    public int getK() {
        return K;
    }

    public int getL() {
        return L;
    }

    public int getN() {
        return N;
    }

    public double getLlh() {
        return llh;
    }

    public double[][][] getMu() {
        return mu;
    }

    public int getn() {
        return n;
    }

    public double[] getPi() {
        return pi;
    }

    public Read[] getReads() {
        return reads;
    }

    public double[][][] getRho() {
        return rho;
    }

    public double getBIC() {
        return BIC;
    }

    public double[][][] getPriorRho() {
        return priorRho;
    }

    public double[] getEps() {
        return eps;
    }

    public int getRestarts() {
        return restarts;
    }
}
