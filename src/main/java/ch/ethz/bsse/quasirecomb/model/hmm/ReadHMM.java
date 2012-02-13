package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.model.Globals;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class ReadHMM {

    private final int L;
    private final int K;
    private final int n;
    private final byte[] read;
    private double[][][] rho;
    private double[] pi;
    private double[][][] mu;
    private double eps;
    private double antieps;
    private double[][][] ufJKV;
    private double[][] fJK;
    private double[][] bJK;
    private double[] c;
    private double[][] gJK;
    private double[][] ugJK;
    private double[][][] gJKV;
    private double[][][] xJKL;
    private double[][][] uxJKL;
    private double[][] ubJK;
    private double[][][] ugJKV;
    private double[][] ufJK;
    private double[][][] fJKV;
    private double backLLH = 0d;

    public ReadHMM(int L, int K, int n, byte[] read, double[][][] rho, double[] pi, double[][][] mu, double eps, double antieps) {
        this.L = L;
        this.K = K;
        this.n = n;
        this.read = read;
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.eps = eps;
        this.antieps = antieps;
        calculate();
    }

    private void calculate() {
        this.forward();
        this.backward();
        this.gamma();
        this.xsi();
        if (Globals.TEST) {
            this.calculateUnscaled();
        }
    }

    public void recalc(double[][][] rho, double[] pi, double[][][] mu, double epsilon) {
        this.rho = rho;
        this.pi = pi;
        this.mu = mu;
        this.eps = epsilon;
        this.antieps = 1 - (n - 1) * epsilon;
        this.calculate();
    }

    private void forward() {
        this.c = new double[L];
        this.fJKV = new double[L][K][n];
        this.fJK = new double[L][K];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
//                double sum = 0d;
                for (int v = 0; v < n; v++) {
                    fJKV[j][k][v] = prRjHv(j, v) * mu[j][k][v];
                    if (j == 0) {
                        fJKV[j][k][v] *= pi[k];
                    } else {
                        double sumL = 0d;
                        for (int l = 0; l < K; l++) {
                            sumL += fJK[j - 1][l] * rho[j - 1][l][k];
                        }
                        fJKV[j][k][v] *= sumL;
                    }
                    c[j] += fJKV[j][k][v];
//                    sum += fJKV[j][k][v];
//                    if (Double.isNaN(fJKV[j][k][v])) {
//                        System.out.println("fJKV");
//                    }
                }
//                if (sum == 0d) {
//                    System.out.println("fJKVsum");
//                    System.err.println("grrrrrrrrrrr");
//                    System.exit(0);
//                }
            }

            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
//                    if (c[j] != 0d) {
                    fJKV[j][k][v] /= c[j];
//                    } else {
//                        System.out.println("c");
//                    }
                    fJK[j][k] += fJKV[j][k][v];
                }
            }
        }
    }

    private void backward() {
        this.bJK = new double[L][K];
        for (int j = L - 1; j >= 0; j--) {
            for (int k = 0; k < K; k++) {
                if (j == L - 1) {
                    bJK[j][k] = 1d / c[L - 1];
                } else {
                    for (int l = 0; l < K; l++) {
                        double sumV = 0d;
                        for (int v = 0; v < n; v++) {
                            sumV += prRjHv(j + 1, v) * mu[j + 1][l][v];
                        }
                        bJK[j][k] += sumV * rho[j][k][l] * bJK[j + 1][l];
                    }
                    if (c[j] != 0d) {
                        bJK[j][k] /= c[j];
                    }
//                    if (Double.isNaN(bJK[j][k])) {
//                        System.out.println("");
//                    }
                }
            }
        }
    }

    private void gamma() {
        this.gJK = new double[L][K];
        this.gJKV = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                this.gJK[j][k] = forward(j, k) * backward(j, k) * c[j];
                for (int v = 0; v < n; v++) {
                    this.gJKV[j][k][v] = forward(j, k, v) * backward(j, k) * c[j];
                }
            }
        }
    }

    private void xsi() {
        this.xJKL = new double[L][K][K];
        for (int j = 1; j < L; j++) {
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    double marginalV = 0d;
                    for (int v = 0; v < n; v++) {
                        marginalV += prRjHv(j, v) * mu[j][l][v];
                    }
                    this.xJKL[j][k][l] = forward(j - 1, k) * marginalV * rho[j - 1][k][l] * backward(j, l);
                }
            }
        }
    }

    public double gamma(int j, int k) {
        return this.gJK[j][k];
    }

    final public double gamma(int j, int k, int v) {
        return this.gJKV[j][k][v];
    }

    final public double xi(int j, int k, int l) {
        if (j < 1) {
            throw new IllegalStateException("J >= 1?");
        }
        return this.xJKL[j][k][l];
    }

    private double prRjHv(int j, int v) {
        return (read[j] == v) ? antieps : eps;
    }

    final public double forward(int j, int k) {
        return this.fJK[j][k];
    }

    final public double forward(int j, int k, int v) {
        return this.fJKV[j][k][v];
    }

    final public double backward(int j, int k) {
        return this.bJK[j][k];
    }

    final public double getC(int j) {
        return c[j];
    }

    final public byte[] getRead() {
        return read;
    }

    final public double[][] getfJK() {
        return fJK;
    }

    final public double[][][] getfJKV() {
        return fJKV;
    }

    final public double[][] getUfJK() {
        return ufJK;
    }

    final public double getBackLLH() {
        return backLLH;
    }

    final public double[][][] getUgJKV() {
        return ugJKV;
    }

    final public double[][][] getUxJKL() {
        return uxJKL;
    }

    private void calculateUnscaled() {
        // backward
        this.ubJK = new double[L][K];
        for (int j = L - 1; j >= 0; j--) {
            for (int k = 0; k < K; k++) {
                if (j == L - 1) {
                    ubJK[j][k] = 1d;
                } else {
                    for (int l = 0; l < K; l++) {
                        double sumV = 0d;
                        for (int v = 0; v < n; v++) {
                            sumV += prRjHv(j + 1, v) * mu[j + 1][l][v];
                        }
                        ubJK[j][k] += sumV * rho[j][l][k] * ubJK[j + 1][l];
                    }
                }
            }
        }
        int j = 0;
        double sum = 0d;
        for (int l = 0; l < K; l++) {
            double sumV = 0d;
            for (int v = 0; v < n; v++) {
                sumV += prRjHv(j, v) * mu[j][l][v];
            }
            sum += sumV * pi[l] * ubJK[j][l];
        }
        this.backLLH = sum;


        // forward
        this.ufJKV = new double[L][K][n];
        this.ufJK = new double[L][K];
        for (int jj = 0; jj < L; jj++) {
            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
                    ufJKV[jj][k][v] = prRjHv(jj, v) * mu[jj][k][v];
                    if (jj == 0) {
                        ufJKV[jj][k][v] *= pi[k];
                    } else {
                        double usumL = 0d;
                        for (int l = 0; l < K; l++) {
                            usumL += ufJK[jj - 1][l] * rho[jj - 1][k][l];
                        }
                        ufJKV[jj][k][v] *= usumL;
                    }
                    ufJK[jj][k] += ufJKV[jj][k][v];
                }
            }
        }

        // gamma
        this.ugJKV = new double[L][K][n];
        this.ugJK = new double[L][K];
        for (int jj = 0; jj < L; jj++) {
            double sumGamma = 0d;
            for (int k = 0; k < K; k++) {
                this.ugJK[jj][k] = this.ufJK[jj][k] * this.ubJK[jj][k];
                for (int v = 0; v < n; v++) {
                    ugJKV[jj][k][v] = this.ufJKV[jj][k][v] * this.ubJK[jj][k];
                    sumGamma += ugJKV[jj][k][v];
                }

            }
            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
                    ugJKV[jj][k][v] /= sumGamma;
                }
            }
        }
        // xsi
        this.uxJKL = new double[L][K][K];
        for (int jj = 1; jj < L; jj++) {
            double sumXsi = 0;
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    double marginalV = 0d;
                    for (int v = 0; v < n; v++) {
                        marginalV += prRjHv(jj, v) * mu[jj][l][v];
                    }
                    this.uxJKL[jj][k][l] = this.ufJK[jj - 1][k] * marginalV * rho[j][k][l] * this.ubJK[jj][l];
                    sumXsi += this.uxJKL[jj][k][l];
                }
            }
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    this.uxJKL[jj][k][l] /= sumXsi;
                }
            }
        }
    }
}
