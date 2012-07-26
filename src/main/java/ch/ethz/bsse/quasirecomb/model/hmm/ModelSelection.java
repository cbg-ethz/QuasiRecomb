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

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;

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
    private double[][][] mu = null;
    private double[][][] rho = null;
    private double[] pi = null;

    public ModelSelection(Read[] reads, int Kmin, int Kmax, int N, int L, int n) {
        this.Kmax = Kmax;
        this.Kmin = Kmin;
        this.N = N;
        this.L = L;
        this.n = n;
        this.start(reads);
    }

    private void start(Read[] reads) {
        double optBIC = 0;
        if (!new File(Globals.SAVEPATH + "support").exists()) {
            new File(Globals.SAVEPATH + "support").mkdirs();
        }
//        Globals.REPEATS = 5;
//        System.out.println("Model selection (" + Globals.REPEATS + " iterations):");
        
        OptimalResult or = null;
        if (Kmin == 0) {
            Globals.MODELSELECTION = true;
            for (int k = 1;; k++) {
                if (!Globals.NO_RECOMB || k == 1) {
                    checkRho0(k);
                }
                EM em = new EM(this.N, this.L, k, this.n, reads);
                if (Globals.LOG_BIC) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(new Summary().print(em.getOr()));
                    Utils.saveFile(Globals.SAVEPATH + "support" + File.separator + "K" + em.getOr().getK() + "-result.txt", sb.toString());
                }
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
            Globals.MODELSELECTION = true;
            for (int k = Kmin; k <= Kmax; k++) {
                if (!Globals.NO_RECOMB || k == 1) {
                    checkRho0(k);
                }
                EM em = new EM(this.N, this.L, k, this.n, reads);
                if (Globals.LOG_BIC) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(new Summary().print(em.getOr()));
                    Utils.saveFile(Globals.SAVEPATH + "support" + File.separator + "K" + em.getOr().getK() + "-result.txt", sb.toString());
                }
                if (em.getOr().getBIC() > optBIC || optBIC == 0) {
                    or = em.getOr();
                    optBIC = em.getOr().getBIC();
                    this.bestK = k;
                    this.mu = em.getMu_opt();
                    this.pi = em.getPi_opt();
                    this.rho = em.getRho_opt();
                }
                Globals.PERCENTAGE = 0;
            }
        } else {
            bestK = Kmin;
        }
        Globals.MODELSELECTION = false;
//        System.out.println("\nBest model: " + bestK);
        Globals.REPEATS = Globals.DESIRED_REPEATS;
        Globals.PERCENTAGE = 0;
//        System.out.println("Model training (" + Globals.REPEATS + " iterations):");
        EM em = new EM(this.N, this.L, bestK, this.n, reads);
        if (em.getOr().getLlh() > optBIC || optBIC == 0) {
            or = em.getOr();
            this.mu = em.getMu_opt();
            this.pi = em.getPi_opt();
            this.rho = em.getRho_opt();
        }
        //save optimumJava
        StringBuilder sb = new StringBuilder();
        sb.append(new Summary().print(or));

        Utils.saveFile(Globals.SAVEPATH + "support" + File.separator + "K" + or.getK() + "-result.txt", sb.toString());
        try {
            String s = Globals.SAVEPATH + "support" + File.separator + "optimumJava";// + (bestK ? "" : K);
            FileOutputStream fos = new FileOutputStream(s);
            try (ObjectOutputStream out = new ObjectOutputStream(fos)) {
                out.writeObject(or);
            }
        } catch (IOException ex) {
            System.out.println("Optimum Java saving\n" + ex.getMessage());
        }

        ModelSampling modelSampling = new ModelSampling(L, n, or.getK(), or.getRho(), or.getPi(), or.getMu(), Globals.SAVEPATH);
        modelSampling.save();
        System.out.println("Quasispecies saved: " + Globals.SAVEPATH + "quasispecies.fasta");
    }

    private static void checkRho0(int K) {
        if (K == 1) {
            Globals.NO_RECOMB = true;
        } else {
            Globals.NO_RECOMB = false;
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
