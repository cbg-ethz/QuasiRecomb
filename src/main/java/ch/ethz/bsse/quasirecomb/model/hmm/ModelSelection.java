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
package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Map;

/**
 * Selects the best model among the specified range of generators.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ModelSelection {

    private int Kmin;
    private int Kmax;
    private int N;
    private int L;
    private int n;
    private int bestK;
    private Map<byte[], Integer> clusterReads;
    private byte[][] haplotypesArray;
    private double[][][] mu = null;
    private double[][][] rho = null;
    private double[] pi = null;

    public ModelSelection(Map<byte[], Integer> clusterReads, int Kmin, int Kmax, int N, int L, int n, byte[][] haplotypesArray) {
        this.Kmax = Kmax;
        this.Kmin = Kmin;
        this.N = N;
        this.L = L;
        this.n = n;
        this.clusterReads = clusterReads;
        this.haplotypesArray = haplotypesArray;
        this.start();
    }

    private void start() {
        double optBIC = 0;

//        Globals.REPEATS = 5;
        System.out.println("Model selection (" + Globals.REPEATS + " iterations):");
        OptimalResult or = null;
        if (Kmin == 0) {
            for (int k = 1;; k++) {
                if (!Globals.rho0force || k == 1) {
                    checkRho0(k);
                }
                EM em = new EM(this.N, this.L, k, this.n, this.clusterReads, this.haplotypesArray);

                if (em.getOr().getBIC() > optBIC || optBIC == 0) {
                    or = em.getOr();
                    optBIC = em.getOr().getBIC();
                    this.bestK = k;
                    this.mu = em.getMu_opt();
                    this.pi = em.getPi_opt();
                    this.rho = em.getRho_opt();
                } else {
                    break;
                }
                Globals.PERCENTAGE = 0;
            }
        } else if (Kmin != Kmax) {
            for (int k = Kmin; k <= Kmax; k++) {
                if (!Globals.rho0force || k == 1) {
                    checkRho0(k);
                }
                EM em = new EM(this.N, this.L, k, this.n, this.clusterReads, this.haplotypesArray);

                if (em.getOr().getBIC() > optBIC || optBIC == 0) {
                    or = em.getOr();
                    optBIC = em.getOr().getBIC();
                    this.bestK = k;
                    this.mu = em.getMu_opt();
                    this.pi = em.getPi_opt();
                    this.rho = em.getRho_opt();
                } else {
                    break;
                }
                Globals.PERCENTAGE = 0;
            }
        } else {
            bestK = Kmin;
        }
        System.out.println("\nBest model: " + bestK);
        Globals.REPEATS = Globals.DESIRED_REPEATS;
        Globals.PERCENTAGE = 0;
        System.out.println("Model training (" + Globals.REPEATS + " iterations):");
        EM em = new EM(this.N, this.L, bestK, this.n, this.clusterReads, this.haplotypesArray);
        if (em.getOr().getBIC() > optBIC || optBIC == 0) {
            or = em.getOr();
            this.mu = em.getMu_opt();
            this.pi = em.getPi_opt();
            this.rho = em.getRho_opt();
        }
        //save optimumJava
        StringBuilder sb = new StringBuilder();
        sb.append(new Summary().print(or));
        if (!new File(Globals.savePath + "support").exists()) {
            new File(Globals.savePath + "support").mkdirs();
        }
        Utils.saveFile(Globals.savePath + "support" + File.separator + "K" + or.getK() + "-result.txt", sb.toString());
        try {
            String s = Globals.savePath + "support" + File.separator + "optimumJava";// + (bestK ? "" : K);
            FileOutputStream fos = new FileOutputStream(s);
            try (ObjectOutputStream out = new ObjectOutputStream(fos)) {
                out.writeObject(or);
            }
        } catch (IOException ex) {
            System.out.println("Optimum Java saving\n" + ex.getMessage());
        }
        ModelSampling modelSampling = new ModelSampling(L, n, or.getK(), or.getRho(), or.getPi(), or.getMu(), Globals.savePath);
        modelSampling.save();
        System.out.println("Quasispecies saved: " + Globals.savePath + "quasispecies.fasta");
    }

    private static void checkRho0(int K) {
        if (K == 1) {
            Globals.rho0 = true;
        } else {
            Globals.rho0 = false;
        }
    }

    public int getBestK() {
        return bestK;
    }

    public double[][][] getMu() {
        return mu;
    }

    public double[] getPi() {
        return pi;
    }

    public double[][][] getRho() {
        return rho;
    }
}
