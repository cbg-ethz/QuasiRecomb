package ch.ethz.bsse.quasirecomb;

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.utils.BitMagic;
import ch.ethz.bsse.quasirecomb.utils.FastaParser;
import ch.ethz.bsse.quasirecomb.utils.Utils;

/**
 * Unit test for simple App.
 */
public class TestLLH {

    /**
     * Rigourous Test :-)
     */
    @org.junit.Test
    public void optimum() {
        int K = 2;
        int L = 19;
        int n = 4;
        double[][] tauOmega = null;
        double[] eps = new double[L];
        double[][][] rho = new double[L - 1][K][K];
        for (int j = 0; j < L - 1; j++) {
            eps[j] = 0;
            for (int k = 0; k < K; k++) {
                for (int l = 0; l < K; l++) {
                    rho[j][k][l] = k == l ? 1 : 0;
                }
            }
        }
        String[] haplotypes = FastaParser.parseFarFile("/Users/XLR/Dropbox/QuasiAsterisk/QuasiRecomb/src/main/resources/haplotypes/datasetLLH.fasta");
        byte[][] byteHaplotypes = new byte[K][L];
        for (int i = 0; i < haplotypes.length; i++) {
            byteHaplotypes[i] = Utils.splitReadIntoByteArray(haplotypes[i]);
        }
        double[][][] mu = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < K; k++) {
                byte b = byteHaplotypes[k][j];
                for (int v = 0; v < n; v++) {
                    mu[j][k][v] = v == b ? 1 : 0;
                }
            }
        }
        mu[1][1] = new double[]{0,0.9,0.1,0};
        mu[4][1] = new double[]{0.7,0,0.3,0};
        double[] pi = new double[]{0.8, 0.2};
        int N = 1000;
        OptimalResult optimalResult = new OptimalResult(N, K, L, n, rho, pi, mu, 0, 0, eps, 0, tauOmega);
        Utils.saveOptimum("/Users/XLR/Dropbox/QuasiAsterisk/QuasiRecomb/src/main/resources/haplotypes/LLH.optimum", optimalResult);
    }

    @org.junit.Test
    public void jumpy() {
        int K = 2;
        int L = 19;
        int n = 4;
        double[][] tauOmega = null;
        double[] eps = new double[L];
        double[][][] rho = new double[L][K][K];
        for (int j = 0; j < L; j++) {
            eps[j] = 0;
        }
        rho[0] = new double[][]{{1, 0}, {0, 1}};
        rho[1] = new double[][]{{1, 0}, {0, 1}};
        rho[2] = new double[][]{{1, 0}, {1, 0}};
        rho[3] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[4] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[5] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[6] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[7] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[8] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[9] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[10] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[11] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[12] = new double[][]{{1, 0}, {0.5, 0.5}};
        rho[13] = new double[][]{{0.8, 0.2}, {0.5, 0.5}};
        rho[14] = new double[][]{{1, 0}, {0, 1}};
        rho[15] = new double[][]{{1, 0}, {0, 1}};
        rho[16] = new double[][]{{1, 0}, {0, 1}};
        rho[17] = new double[][]{{1, 0}, {0, 1}};
        String[] haplotypes = FastaParser.parseFarFile("/Users/XLR/Dropbox/QuasiAsterisk/QuasiRecomb/src/main/resources/haplotypes/datasetLLH.fasta");
        byte[][] byteHaplotypes = new byte[K][L];
        for (int i = 0; i < haplotypes.length; i++) {
            byteHaplotypes[i] = Utils.splitReadIntoByteArray(haplotypes[i]);
        }
        double[][][] mu = new double[L][K][n];
        for (int j = 0; j < L; j++) {
            int k = 0;
            byte b = byteHaplotypes[k][j];
            for (int v = 0; v < n; v++) {
                mu[j][k][v] = v == b ? 1 : 0;
            }

        }
        for (int j = 0; j < L; j++) {
            int k = 1;
            byte b = byteHaplotypes[k][j];
            for (int v = 0; v < n; v++) {
                if (j < 2 || j > 13) {
                    mu[j][k][v] = v == b ? 1 : 0;
                } else {
                    mu[j][k][v] = 1d/n;
                }
            }

        }
        for (int j = 0; j < L; j++) {
            int k = 0;
            byte b = byteHaplotypes[k][j];
            for (int v = 0; v < n; v++) {
                mu[j][k][v] = v == b ? 1 : 0;
            }

        }
        double[] pi = new double[]{0.8, 0.2};
        int N = 1000;
        OptimalResult optimalResult = new OptimalResult(N, K, L, n, rho, pi, mu, 0, 0, eps, 0, tauOmega);
        Utils.saveOptimum("/Users/XLR/Dropbox/QuasiAsterisk/QuasiRecomb/src/main/resources/haplotypes/LLH_jump.optimum", optimalResult);
    }
}
