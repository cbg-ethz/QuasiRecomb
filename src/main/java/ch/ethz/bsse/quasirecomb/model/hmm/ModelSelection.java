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

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;

/**
 * Selects the best model among the specified range of generators.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ModelSelection {

    private static void checkRho0(int K) {
        if (K == 1) {
            Globals.getINSTANCE().setNO_RECOMB(true);
        } else {
            Globals.getINSTANCE().setNO_RECOMB(false);
        }
    }

    private int kMin;
    private int kMax;
    private int N;
    private int L;
    private int n;
    private int bestK;
    private OptimalResult or;

    public ModelSelection(Read[] reads, int Kmin, int Kmax, int N, int L, int n) {
        this.kMax = Kmax;
        this.kMin = Kmin;
        this.N = N;
        this.L = L;
        this.n = n;
        this.start(reads);
    }

    private void start(Read[] reads) {
        double optBIC = 0;
        String save = Globals.getINSTANCE().getSAVEPATH() + "support";
        Utils.mkdir(Globals.getINSTANCE().getSnapshotDir());

        if (kMin != kMax) {
            Globals.getINSTANCE().setMODELSELECTION(true);
            Utils.mkdir(Globals.getINSTANCE().getSnapshotDir() + File.separator + "modelselection");
            for (int k = kMin; k <= kMax; k++) {
                if (!Globals.getINSTANCE().isFORCE_NO_RECOMB()) {
                    checkRho0(k);
                }
                EM em = new EM(this.N, this.L, k, this.n, reads);
                System.out.println("");
                if (Globals.getINSTANCE().isLOG_BIC()) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(new Summary().print(em.getOr()));
                    Utils.saveFile(save + File.separator + "K" + em.getOr().getK() + "-result.txt", sb.toString());
                }
                if (em.getOr().getBIC() > optBIC || optBIC == 0) {
                    or = em.getOr();
                    optBIC = em.getOr().getBIC();
                    this.bestK = k;
                } else {
                    break;
                }
                Globals.getINSTANCE().setPERCENTAGE(0);
            }
        } else {
            bestK = kMin;
        }
        Globals.getINSTANCE().setMODELSELECTION(false);
        Utils.mkdir(Globals.getINSTANCE().getSnapshotDir() + File.separator + "training");
        Globals.getINSTANCE().setREPEATS(Globals.getINSTANCE().getDESIRED_REPEATS());
        Globals.getINSTANCE().setPERCENTAGE(0);
        EM em = new EM(this.N, this.L, bestK, this.n, reads);
        if (em.getOr().getLlh() > optBIC || optBIC == 0) {
            or = em.getOr();
        }

        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + File.separator + "K" + or.getK() + "-result.txt", new Summary().print(or));
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + File.separator + "K" + or.getK() + "-minimal.txt", new Summary().minimal(or));
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + File.separator + "K" + or.getK() + "-summary.html", new Summary().html(or));
        //save optimumJava
        Utils.saveOptimum(save + File.separator + "best.optimum", or);
    }

    public int getBestK() {
        return bestK;
    }

    public OptimalResult getOptimalResult() {
        return or;
    }
}