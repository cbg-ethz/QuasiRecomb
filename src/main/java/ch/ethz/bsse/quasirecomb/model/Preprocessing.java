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
package ch.ethz.bsse.quasirecomb.model;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.ModelSelection;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.utils.BitMagic;
import ch.ethz.bsse.quasirecomb.utils.Plot;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.util.HashMap;
import java.util.Map;
import org.apache.commons.math.stat.descriptive.moment.Mean;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;

/**
 * Responsible for orchestrating parsing, proper read placement in alignment and
 * forwards parameters to ModelSelection.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Preprocessing {

    /**
     * Entry point. Forwards invokes of the specified workflow.
     *
     * @param L length of the reads
     * @param exp is it an experimental dataset
     * @param input path to the fasta file
     * @param Kmin minimal amount of generators
     * @param Kmax minimal amount of generators
     * @param n size of the alphabet
     * @param f distribution of the haplotypes if sampling has to be done
     * @param N amount of reads in case exp is false
     */
    public static void workflow(String input, int Kmin, int Kmax) {
        Utils.mkdir(Globals.getINSTANCE().getSAVEPATH() + "support");
        Globals.getINSTANCE().print("Parsing");
        Read[] reads = Utils.parseInput(input);

        for (Read r : reads) {
            Globals.getINSTANCE().setALIGNMENT_BEGIN(Math.min(r.getBegin(), Globals.getINSTANCE().getALIGNMENT_BEGIN()));
            Globals.getINSTANCE().setALIGNMENT_END(Math.max(r.getEnd(), Globals.getINSTANCE().getALIGNMENT_END()));
        }
        int L = Globals.getINSTANCE().getALIGNMENT_END() - Globals.getINSTANCE().getALIGNMENT_BEGIN();
        Globals.getINSTANCE().setALIGNMENT_END(L);
        for (Read r : reads) {
            r.shrink();
        }

        Globals.getINSTANCE().print("Parsing\t25%");
        int[][] alignment = countPos(reads, L);


        Globals.getINSTANCE().print("Parsing\t50%");
        StringBuilder sb = new StringBuilder();
        sb.append("Start: ").append(Globals.getINSTANCE().getALIGNMENT_BEGIN()).append("\n");
        for (int i = 0; i < L; i++) {
            int hits = 0;
            sb.append(i);
            for (int v = 0; v < alignment[i].length; v++) {
                sb.append("\t").append(alignment[i][v]);
                if (alignment[i][v] != 0) {
                    hits++;
                }
            }
            sb.append("\n");
            if (hits == 0) {
                System.out.println("Position " + i + " is not covered.");
            }
        }
        Globals.getINSTANCE().setMU_PRIOR(alignment);
        Globals.getINSTANCE().print("Parsing\t75%");
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "hit_dist.txt", sb.toString());
        sb = null;
        System.gc();
        System.gc();
        int n = countChars(reads);
        Globals.getINSTANCE().print("Parsing\t100%");
        Globals.getINSTANCE().println("Unique reads\t" + reads.length);
        Globals.getINSTANCE().println("Paired reads\t" + Globals.getINSTANCE().getPAIRED_COUNT());
        if (Globals.getINSTANCE().getPAIRED_COUNT() > 0) {
            double[] inserts = new double[Globals.getINSTANCE().getPAIRED_COUNT()];
            int x = 0;
            for (Read r : reads) {
                if (r.isPaired()) {
                    inserts[x++] = r.getCrickBegin() - r.getWatsonEnd();
                }
            }
            Globals.getINSTANCE().println("Insert size\t" + Math.round((new Mean().evaluate(inserts))*10)/10 + " (±" + Math.round(new StandardDeviation().evaluate(inserts)*10)/10 + ")");
        }
        Globals.getINSTANCE().println("Merged reads\t" + Globals.getINSTANCE().getMERGED() + "\n");
        if (Globals.getINSTANCE().isPLOT()) {
            Globals.getINSTANCE().println("Plotting\t");
            Plot.plotCoverage(alignment);
        }
        if (Globals.getINSTANCE().isPRINT_ALIGNMENT()) {
            Globals.getINSTANCE().println("Saving alignment\t");
            new Summary().printAlignment(reads);
        }
        if (Globals.getINSTANCE().isCIRCOS()) {
            new Summary().circos(L, alignment);
            System.exit(0);
        }
        if (Globals.getINSTANCE().isDEBUG()) {
            new File(Globals.getINSTANCE().getSAVEPATH() + "support/log/").mkdirs();
        }

        ModelSelection ms = new ModelSelection(reads, Kmin, Kmax, reads.length, L, n);
        if (!Globals.getINSTANCE().isNOSAMPLE()) {
            ModelSampling modelSampling = new ModelSampling(ms.getOptimalResult(), Globals.getINSTANCE().getSAVEPATH());
            modelSampling.save();
            System.out.println("\nQuasispecies saved: " + Globals.getINSTANCE().getSAVEPATH() + "quasispecies.fasta");
        }
    }

    private static int countChars(Read[] rs) {
        Map<Byte, Boolean> map = new HashMap<>();
        for (Read r : rs) {
            if (r.isPaired()) {
                for (int i = 0; i < r.getLength(); i++) {
                    if (r.isHit(i)) {
                        map.put(r.getBase(i), Boolean.TRUE);
                    }
                }
            } else {
                for (int i = 0; i < r.getWatsonLength(); i++) {
                    map.put(r.getBase(i), Boolean.TRUE);
                }
            }
        }
        return map.keySet().size();
    }

    private static void saveUnique(Read[] reads) {
        if (Globals.getINSTANCE().isDEBUG()) {
            StringBuilder sb = new StringBuilder();
            for (Read r : reads) {
                sb.append(r.getCount()).append("\t");
                if (r.getCount() < 1000) {
                    sb.append("\t");
                }
                for (int i = 0; i < r.getBegin(); i++) {
                    sb.append(" ");
                }
                sb.append(Utils.reverse(r.getSequence())).append("\n");
            }
            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + File.separator + "in.fasta", sb.toString());
        }
    }

    public static int[][] countPos(Read[] reads, int L) {
        int[][] alignment = new int[L][5];
        for (Read r : reads) {
            int begin = r.getWatsonBegin();
            for (int i = 0; i < r.getWatsonLength(); i++) {
                try {
                    alignment[i + begin][BitMagic.getPosition(r.getSequence(), i)] += r.getCount();
                } catch (ArrayIndexOutOfBoundsException e) {
                    System.out.println(e);
                }
            }
            if (r.isPaired()) {
                begin = r.getCrickBegin();
                for (int i = 0; i < r.getCrickLength(); i++) {
                    alignment[i + begin][BitMagic.getPosition(r.getCrickSequence(), i)] += r.getCount();
                }
            }
        }
        return alignment;
    }
}
