/**
 * Copyright (c) 2011-2013 Armin Töpfer
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
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math.stat.descriptive.rank.Median;

/**
 * Responsible for the start of the repeats within the EM algorithm.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class EM extends Utils {

    private OptimalResult or;
    private double medianBIC;
    private double lowerBoundBIC;

    protected EM(int N, int L, int K, int n, Read[] reads) {
        this.blackbox(reads, N, L, K, n);
    }

    private void blackbox(Read[] reads, int N, int L, int K, int n) {
        Globals.getINSTANCE().setLOG(new StringBuilder());
        Globals.getINSTANCE().setMAX_LLH(-1);
        Globals.getINSTANCE().setMIN_BIC(Double.MAX_VALUE);
        String pathOptimum = null;
        if (K == 1 || Globals.getINSTANCE().isFORCE_NO_RECOMB()) {
            Globals.getINSTANCE().setNO_RECOMB(true);
        } else {
            Globals.getINSTANCE().setNO_RECOMB(false);
        }
        if (Globals.getINSTANCE().getOPTIMUM() == null) {
            double maxLLH = Double.NEGATIVE_INFINITY;
            double bics[] = new double[Globals.getINSTANCE().getREPEATS()];
            for (int i = 0; i < Globals.getINSTANCE().getREPEATS(); i++) {
                SingleEM sem = new SingleEM(N, K, L, n, reads, Globals.getINSTANCE().getDELTA_LLH(), i);
                bics[i] = sem.getOptimalResult().getBIC();
                if (sem.getLoglikelihood() > maxLLH) {
                    maxLLH = sem.getLoglikelihood();
                    pathOptimum = sem.getOptimumPath();
                }
            }
            medianBIC = new Median().evaluate(bics);
            lowerBoundBIC = medianBIC - new StandardDeviation().evaluate(bics) * Math.sqrt(1 + 1d / bics.length);
        } else {
            pathOptimum = Globals.getINSTANCE().getOPTIMUM();
        }
        if (Globals.getINSTANCE().isMODELSELECTION()) {
            try {
                FileInputStream fis = new FileInputStream(pathOptimum);
                try (ObjectInputStream in = new ObjectInputStream(fis)) {
                    or = (OptimalResult) in.readObject();
                }
            } catch (IOException | ClassNotFoundException ex) {
                System.err.println(ex);
            }
            Globals.getINSTANCE().printBIC(K, (int) or.getBIC());
        } else {
            long time = System.currentTimeMillis();

            try {
                FileInputStream fis = new FileInputStream(pathOptimum);
                try (ObjectInputStream in = new ObjectInputStream(fis)) {
                    or = (OptimalResult) in.readObject();
                }
            } catch (IOException | ClassNotFoundException ex) {
                System.err.println(ex);
            }

            Globals.getINSTANCE().printBIC(K, (int) or.getBIC());
            System.out.print("\n");
            SingleEM bestEM = new SingleEM(or, Globals.getINSTANCE().getDELTA_REFINE_LLH(), reads);
            this.or = bestEM.getOptimalResult();
            if (Globals.getINSTANCE().isLOGGING()) {
                Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "log_K" + K, Globals.getINSTANCE().getLOG().toString());
            }
        }
    }

    /**
     * The double 3D array of mu for the best result.
     *
     * @return mu
     */
    public double[][][] getMu_opt() {
        return or.getMu();
    }

    /**
     * The double array of pi for the best result.
     *
     * @return pi
     */
    public double[] getPi_opt() {
        return or.getPi();
    }

    /**
     * The double 3D array of rho for the best result.
     *
     * @return rho
     */
    public double[][][] getRho_opt() {
        return or.getRho();
    }

    public OptimalResult getOr() {
        return or;
    }

    public double getMedianBIC() {
        return medianBIC;
    }

    public double getLowerBoundBIC() {
        return lowerBoundBIC;
    }
}
