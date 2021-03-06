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
import ch.ethz.bsse.quasirecomb.informationholder.ModelSelectionStorage;
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.utils.StatusUpdate;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

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
    private ModelSelectionStorage msTemp = new ModelSelectionStorage();
    private SortedMap<Integer, Double[]> bics = new TreeMap<>();

    public ModelSelection(Read[] reads, int Kmin, int Kmax, int L, int n) {
        this.kMax = Kmax;
        this.kMin = Kmin;
        this.L = L;
        this.n = n;
        Globals.getINSTANCE().setREPEATS(Globals.getINSTANCE().getMS_REPEATS());
        this.start(reads);
    }

    private void start(Read[] reads) {
        double optBIC = 0;
        String save = Globals.getINSTANCE().getSAVEPATH() + "support";
        Utils.mkdir(Globals.getINSTANCE().getSnapshotDir());

        if (Globals.getINSTANCE().getSUB_SAMPLE_PERC() > 0d) {
            shuffleArray(reads);
            reads = Arrays.copyOfRange(reads, 0, (int) Math.round(Globals.getINSTANCE().getSUB_SAMPLE_PERC() * reads.length));
        }

        for (Read r : reads) {
            this.N += r.getCount();
        }
        StatusUpdate.getINSTANCE().println("Unique reads\t" + reads.length);
        StatusUpdate.getINSTANCE().println("Total reads\t" + this.N + "\n");
        if (!Globals.getINSTANCE().isUSER_OPTIMUM()) {
            if (kMin != kMax) {
                Globals.getINSTANCE().setMODELSELECTION(true);
                Utils.mkdir(Globals.getINSTANCE().getSnapshotDir() + File.separator + "modelselection");
                select(reads, save);
                saveBics();
                Utils.saveR();
            } else {
                bestK = kMin;
            }
        }
        Globals.getINSTANCE().setMODELSELECTION(false);
//        if (!Globals.getINSTANCE().isBOOTSTRAP()) {
        Utils.mkdir(Globals.getINSTANCE().getSnapshotDir() + File.separator + "training");
        Globals.getINSTANCE().setREPEATS(Globals.getINSTANCE().getDESIRED_REPEATS());

        if (Globals.getINSTANCE().getDESIRED_REPEATS() > 0) {
            StatusUpdate.getINSTANCE().setPERCENTAGE(0);
            if (Globals.getINSTANCE().isSUBSAMPLE()) {
                shuffleArray(reads);
                List<Read> subsample = new LinkedList<>();
                Map<String, Integer> generators = new HashMap<>();
                Globals.getINSTANCE().setDESIRED_REPEATS(0);
                double mult_mu = Globals.getINSTANCE().getMULT_MU();
                double mult_rho = Globals.getINSTANCE().getMULT_RHO();
                Globals.getINSTANCE().setMULT_MU(1);
                Globals.getINSTANCE().setMULT_RHO(1);
                for (int i = 0; i < reads.length; i++) {
                    if (subsample.size() < reads.length / Globals.getINSTANCE().getSUBSAMPLE_COUNT() || i + reads.length / Globals.getINSTANCE().getSUBSAMPLE_COUNT() > reads.length) {
                        subsample.add(reads[i]);
                    } else {
                        Read[] readsSubsample = subsample.toArray(new Read[subsample.size()]);
                        Globals.getINSTANCE().setREPEATS(Globals.getINSTANCE().getMS_REPEATS());
                        EM em = new EM(this.N, this.L, bestK, this.n, readsSubsample);
                        subsample.clear();
                        saveSubSample(em.getOr().getMu(), generators, L);
                    }
                }
                for (Map.Entry<String, Integer> e : generators.entrySet()) {
                    System.out.println(e.getValue() + "\t" + e.getKey());
                }
                Globals.getINSTANCE().setREPEATS(Globals.getINSTANCE().getDESIRED_REPEATS());
                System.exit(0);
            }

            EM em = new EM(this.N, this.L, bestK, this.n, reads);
            if (em.getOr().getLlh() > optBIC || optBIC == 0) {
                or = em.getOr();
            }

            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "K" + or.getK() + "-result.txt", new Summary().print(or));
            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "K" + or.getK() + "-minimal.txt", new Summary().minimal(or));
            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "K" + or.getK() + "-summary.html", new Summary().html(or));
            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "K" + or.getK() + "-snvs.txt", new Summary().snvs(or));
            //save optimumJava
            Utils.saveOptimum(save + File.separator + "best.optimum", or);
        }
//        }
    }

    private static void saveSubSample(double[][][] mu, Map<String, Integer> generators, int L) {
        int K = mu[0].length;
        for (int k = 0; k < K; k++) {
            StringBuilder gen = new StringBuilder();
            for (int j = 0; j < L; j++) {
                gen.append(max(mu[j][k]));
            }
            String generator = gen.toString();
            int count = generators.containsKey(generator) ? generators.get(generator) : 0;
            generators.put(generator, count + 1);
        }
    }

    private static String max(double[] array) {
        int index = -1;
        double max = -1;
        for (int i = 0; i < array.length; i++) {
            if (array[i] > max) {
                index = i;
                max = array[i];
            }
        }
        return Utils.reverse(index);
    }

    private static void shuffleArray(Read[] ar) {
        Random rnd = new Random();
        for (int i = ar.length - 1; i >= 0; i--) {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            Read a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }

    public int getBestK() {
        return bestK;
    }

    public OptimalResult getOptimalResult() {
        return or;
    }

    private void select(Read[] reads, String save) {
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
            bics.put(k, em.getBics());
            msTemp.add(em, k);
            StatusUpdate.getINSTANCE().setPERCENTAGE(0);
            or = msTemp.getBestOR();
            if (!Globals.getINSTANCE().isBOOTSTRAP() && or.getK() < k && or.getK() + 1 == k) {
                break;
            }
        }
        bestK = or.getK();
    }

    public void saveBics() {
        StringBuilder sb = new StringBuilder();
        int x = bics.values().iterator().next().length;
        Set<Integer> keySet = bics.keySet();
        for (int i : keySet) {
            sb.append(i).append("\t");
        }
        sb.setLength(sb.length() - 1);
        sb.append("\n");

        for (int l = 0; l < x; l++) {
            for (int i : keySet) {
                sb.append(bics.get(i)[l]).append("\t");
            }
            sb.setLength(sb.length() - 1);
            sb.append("\n");
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "bics.txt", sb.toString());
    }

    public ModelSelectionStorage getMsTemp() {
        return msTemp;
    }
}
