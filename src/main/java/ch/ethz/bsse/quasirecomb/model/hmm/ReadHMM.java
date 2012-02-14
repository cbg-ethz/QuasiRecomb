package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.model.Globals;

/**
 * Calculates the forward / backward algorithm for a single read.
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
    private int begin;
    private int end;
    private int length;

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

        for (int j = 0; j < L; j++) {
            if (read[j] != 4) {
                begin = j;
                break;
            }
        }
        for (int j = L - 1; j >= 0; j--) {
            if (read[j] != 4) {
                end = j + 1;
                break;
            }
        }

        this.length = end - begin;

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
        this.c = new double[length];
        this.fJKV = new double[length][K][n];
        this.fJK = new double[length][K];
        for (int j = 0; j < length; j++) {
            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
                    fJKV[j][k][v] = prRjHv(begin + j, v) * mu[begin + j][k][v];
                    if (j == 0) {
                        fJKV[j][k][v] *= pi[k];
                    } else {
                        double sumL = 0d;
                        for (int l = 0; l < K; l++) {
                            sumL += fJK[j - 1][l] * rho[begin + j - 1][l][k];
                        }
                        fJKV[j][k][v] *= sumL;
                    }
                    c[j] += fJKV[j][k][v];
                }
            }

            for (int k = 0; k < K; k++) {
                for (int v = 0; v < n; v++) {
                    fJKV[j][k][v] /= c[j];
                    fJK[j][k] += fJKV[j][k][v];
                }
            }
        }
    }

    private void backward() {
        this.bJK = new double[length][K];
        for (int j = length - 1; j >= 0; j--) {
            for (int k = 0; k < K; k++) {
                if (j == length - 1) {
                    bJK[j][k] = 1d / c[length - 1];
                } else {
                    for (int l = 0; l < K; l++) {
                        double sumV = 0d;
                        for (int v = 0; v < n; v++) {
                            sumV += prRjHv(begin + j + 1, v) * mu[begin + j + 1][l][v];
                        }
                        bJK[j][k] += sumV * rho[begin + j][k][l] * bJK[j + 1][l];
                    }
                    if (c[j] != 0d) {
                        bJK[j][k] /= c[j];
                    }
                }
            }
        }
    }

    private void gamma() {
        this.gJK = new double[length][K];
        this.gJKV = new double[length][K][n];
        for (int j = 0; j < length; j++) {
            for (int k = 0; k < K; k++) {
                this.gJK[j][k] = this.fJK[j][k] * this.bJK[j][k] * c[j];
                for (int v = 0; v < n; v++) {
                    this.gJKV[j][k][v] = this.fJKV[j][k][v] * this.bJK[j][k] * c[j];
                }
            }
        }
    }

    private void xsi() {
        this.xJKL = new double[length][K][K];
        for (int j = 1; j < length; j++) {
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    double marginalV = 0d;
                    for (int v = 0; v < n; v++) {
                        marginalV += prRjHv(j, v) * mu[begin + j][l][v];
                    }
                    this.xJKL[j][k][l] = this.fJK[j - 1][k] * marginalV * rho[begin + j - 1][k][l] * this.bJK[j][l];
                }
            }
        }
    }

    public double gamma(int j, int k) {
        if (j >= begin && j - begin < length) {
            return this.gJK[j - begin][k];
        }
        return 0;
    }

    final public double gamma(int j, int k, int v) {
        if (j >= begin && j - begin < length) {
            return this.gJKV[j - begin][k][v];
        }
        return 0;
    }

    final public double xi(int j, int k, int l) {
//        if (j < 1) {
//            throw new IllegalStateException("J >= 1?");
//        }
        if (j >= begin && j - begin < length) {
            return this.xJKL[j - begin][k][l];
        }
        return 0;
    }

    private double prRjHv(int j, int v) {
        return (read[j] == v) ? antieps : eps;
    }

    final public double getC(int j) {
        if (j >= begin && j - begin < length) {
            return c[j - begin];
        }
        return 1;
    }

    final public byte[] getRead() {
        return read;
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
//        this.backLLH = sum;


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
