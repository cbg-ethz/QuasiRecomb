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
import ch.ethz.bsse.quasirecomb.informationholder.MSBTemp;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.ModelSelection;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.utils.BitMagic;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
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
        //parse file
        Globals.getINSTANCE().print("Parsing");
        Read[] reads = Utils.parseInput(input);
        int L = fixAlignment(reads);
        int[][] alignment = computeAlignment(reads, L);
        Globals.getINSTANCE().println("Unique reads\t" + reads.length);
        Globals.getINSTANCE().println("Paired reads\t" + Globals.getINSTANCE().getPAIRED_COUNT());
        computeInsertDist(reads);
        Globals.getINSTANCE().println("Merged reads\t" + Globals.getINSTANCE().getMERGED() + "\n");
        printAlignment(reads);
        circos(L, alignment);
        if (Globals.getINSTANCE().isDEBUG()) {
            new File(Globals.getINSTANCE().getSAVEPATH() + "support/log/").mkdirs();
        }
        int n = countChars(reads);
        Globals.getINSTANCE().setTAU_OMEGA(reads, L);
        plot();

        double N = 0;
        for (Read r : reads) {
            N += r.getCount();
        }

        if (Globals.getINSTANCE().isBOOTSTRAP()) {
            Multimap<Integer, Double> bics = ArrayListMultimap.create();
            Map<Read, Double> piMap = new HashMap<>();
            for (Read r : reads) {
                piMap.put(r, r.getCount() / N);
            }
            Frequency<Read> readDist = new Frequency<>(piMap);

            for (int i = 0; i < 10; i++) {
                Globals.getINSTANCE().println("Bootstrap " + i + "\n");
                Map<Integer, Read> hashed = new HashMap<>();

                for (int x = 0; x < N; x++) {
                    Read r = readDist.roll();
                    int hash = r.hashCode();
                    if (hashed.containsKey(hash)) {
                        hashed.get(hash).incCount();
                    } else {
                        hashed.put(hash, new Read(r));
                    }
                }

                Read[] rs = hashed.values().toArray(new Read[hashed.values().size()]);
                ModelSelection ms = new ModelSelection(rs, Kmin, Kmax, rs.length, L, n);
                bics.putAll(ms.getMsTemp().getMaxBICs());
            }
            StringBuilder sb = new StringBuilder();

            int x = bics.asMap().values().iterator().next().size();
            Set<Integer> keySet = bics.keySet();
            for (int i : keySet) {
                sb.append(i).append("\t");
            }
            sb.setLength(sb.length() - 1);
            sb.append("\n");



            for (int l = 0; l < x; l++) {
                for (int i : keySet) {
                    ArrayList arrayList = new ArrayList(bics.get(i));
                    if (l < arrayList.size()) {
                        sb.append(arrayList.get(l));
                    }
                    sb.append("\t");
                }
                sb.setLength(sb.length() - 1);
                sb.append("\n");
            }
            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "bics.txt", sb.toString());
            Globals.getINSTANCE().println("");
            MSBTemp msbt = new MSBTemp(bics);
            Kmin = msbt.getBestK();
            Kmax = Kmin;
            Globals.getINSTANCE().setBOOTSTRAP(false);
        }

        ModelSelection ms = new ModelSelection(reads, Kmin, Kmax, reads.length, L, n);

        if (!Globals.getINSTANCE().isNOSAMPLE()) {
            ModelSampling modelSampling = new ModelSampling(ms.getOptimalResult(), Globals.getINSTANCE().getSAVEPATH());
            modelSampling.save();
            System.out.println("\nQuasispecies saved: " + Globals.getINSTANCE().getSAVEPATH() + "quasispecies.fasta");
        }
        if (!Globals.getINSTANCE().isDEBUG()) {
            deleteDirectory(new File(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "snapshots"));
        }
    }

    static public boolean deleteDirectory(File path) {
        if (path.exists()) {
            File[] files = path.listFiles();
            for (int i = 0; i < files.length; i++) {
                if (files[i].isDirectory()) {
                    deleteDirectory(files[i]);
                } else {
                    files[i].delete();
                }
            }
        }
        return (path.delete());
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

//    private static void saveUnique(Read[] reads) {
//        if (Globals.getINSTANCE().isDEBUG()) {
//            StringBuilder sb = new StringBuilder();
//            for (Read r : reads) {
//                sb.append(r.getCount()).append("\t");
//                if (r.getCount() < 1000) {
//                    sb.append("\t");
//                }
//                for (int i = 0; i < r.getBegin(); i++) {
//                    sb.append(" ");
//                }
//                sb.append(Utils.reverse(r.getSequence())).append("\n");
//            }
//            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + File.separator + "in.fasta", sb.toString());
//        }
//    }
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

    public static double[][] countPosWeighted(Read[] reads, int L) {
        double[][] alignment = new double[L][5];
        for (Read r : reads) {
            int begin = r.getWatsonBegin();
            for (int i = 0; i < r.getWatsonLength(); i++) {
                try {
                    if (r.getWatsonQuality() == null || r.getWatsonQuality().length > 0) {
                        alignment[i + begin][BitMagic.getPosition(r.getSequence(), i)] += r.getCount();
                    } else {
                        alignment[i + begin][BitMagic.getPosition(r.getSequence(), i)] += r.getCount() * r.getWatsonQuality()[i];
                    }
                } catch (ArrayIndexOutOfBoundsException e) {
                    System.out.println(e);
                }
            }
            if (r.isPaired()) {
                begin = r.getCrickBegin();
                for (int i = 0; i < r.getCrickLength(); i++) {
                    if (r.getCrickQuality() == null || r.getCrickQuality().length > 0) {
                        alignment[i + begin][BitMagic.getPosition(r.getCrickSequence(), i)] += r.getCount();
                    } else {
                        alignment[i + begin][BitMagic.getPosition(r.getCrickSequence(), i)] += r.getCount() * r.getCrickQuality()[i];
                    }
                }
            }
        }

        return alignment;
    }

    private static int fixAlignment(Read[] reads) {
        //fix alignment to position 0
        int N = 0;
        for (Read r : reads) {
            if (r.isPaired()) {
                Globals.getINSTANCE().setPAIRED(true);
            }
            Globals.getINSTANCE().setALIGNMENT_BEGIN(Math.min(r.getBegin(), Globals.getINSTANCE().getALIGNMENT_BEGIN()));
            Globals.getINSTANCE().setALIGNMENT_END(Math.max(r.getEnd(), Globals.getINSTANCE().getALIGNMENT_END()));
            N += r.getCount();
        }
        Globals.getINSTANCE().setNREAL(N);
        int L = Globals.getINSTANCE().getALIGNMENT_END() - Globals.getINSTANCE().getALIGNMENT_BEGIN();
        Globals.getINSTANCE().setALIGNMENT_END(L);
        Globals.getINSTANCE().println("Modifying reads\t");
        double shrinkCounter = 0;
        for (Read r : reads) {
            r.shrink();
            Globals.getINSTANCE().print("Modifying reads\t" + (Math.round((shrinkCounter++ / reads.length) * 100)) + "%");
        }
        Globals.getINSTANCE().print("Modifying reads\t100%");
        return L;
    }

    private static void computeAllelFrequencies(int L, int[][] alignment, double[][] alignmentWeighted) {
        Globals.getINSTANCE().println("Allel frequencies\t");
        double allelCounter = 0;
        StringBuilder sb = new StringBuilder();
        StringBuilder sbw = new StringBuilder();
        char[] alphabet = new char[]{'A', 'C', 'G', 'T', '-'};
        sb.append("#Offset: ").append(Globals.getINSTANCE().getALIGNMENT_BEGIN()).append("\n");
        sbw.append("#Offset: ").append(Globals.getINSTANCE().getALIGNMENT_BEGIN()).append("\n");
        sb.append("Pos");
        sbw.append("Pos");
        for (int i = 0; i < alignment[0].length; i++) {
            sb.append("\t").append(alphabet[i]);
            sbw.append("\t").append(alphabet[i]);
        }
        sb.append("\n");
        sbw.append("\n");
        for (int i = 0; i < L; i++) {
            int hits = 0;
            sb.append(i);
            sbw.append(i);
            for (int v = 0; v < alignment[i].length; v++) {
                sb.append("\t").append(alignment[i][v]);
                sbw.append("\t").append(Summary.shorten(alignmentWeighted[i][v]));
                if (alignment[i][v] != 0) {
                    hits++;
                }
            }
            sb.append("\n");
            sbw.append("\n");
            if (hits == 0) {
                System.out.println("Position " + i + " is not covered.");
            }
            Globals.getINSTANCE().print("Allel frequencies\t" + (Math.round((allelCounter++ / L) * 100)) + "%");
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "allel_distribution.txt", sb.toString());
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "allel_distribution_phred_weighted.txt", sbw.toString());
        sb.setLength(0);
        sb = null;
    }

    private static int[][] computeAlignment(Read[] reads, int L) {
        Globals.getINSTANCE().println("Computing entropy\t");
        double entropyCounter = 0;
        int[][] alignment = countPos(reads, L);
        double[][] alignmentWeighted = countPosWeighted(reads, L);
        double alignmentEntropy = 0;
        for (int i = 0; i < L; i++) {
            double sum = 0;
            for (int j = 0; j < 5; j++) {
                sum += alignmentWeighted[i][j];
            }
            double shannonEntropy_pos = 0d;
            for (int j = 0; j < 5; j++) {
                alignmentWeighted[i][j] /= sum;
                if (alignmentWeighted[i][j] > 0) {
                    shannonEntropy_pos -= alignmentWeighted[i][j] * Math.log(alignmentWeighted[i][j]) / Math.log(5);
                }
            }
            alignmentEntropy += shannonEntropy_pos;
            Globals.getINSTANCE().print("Computing entropy\t" + (Math.round((entropyCounter++ / L) * 100)) + "%");
        }
        Globals.getINSTANCE().print("Computing entropy\t100%");
        alignmentEntropy /= L;
        computeAllelFrequencies(L, alignment, alignmentWeighted);
        Globals.getINSTANCE().print("Allel frequencies\t100%");
        Globals.getINSTANCE().println("Alignment entropy\t" + Math.round(alignmentEntropy*1000)/1000d);
        return alignment;
    }

    private static void computeInsertDist(Read[] reads) {
        if (Globals.getINSTANCE().getPAIRED_COUNT() > 0) {
            double[] inserts = new double[Globals.getINSTANCE().getPAIRED_COUNT()];
            int x = 0;
            for (Read r : reads) {
                if (r.isPaired()) {
                    inserts[x++] = r.getCrickBegin() - r.getWatsonEnd();
                }
            }
            Globals.getINSTANCE().println("Insert size\t" + Math.round((new Mean().evaluate(inserts)) * 10) / 10 + " (±" + Math.round(new StandardDeviation().evaluate(inserts) * 10) / 10 + ")");
        }
    }

    private static void plot() {
        Globals.getINSTANCE().println("Compute coverage\t");
        StringBuilder sb = new StringBuilder();
        int[] coverage = Globals.getINSTANCE().getTAU_OMEGA().getCoverage();
        {
            int start = Globals.getINSTANCE().getALIGNMENT_BEGIN();
            for (int i = 0; i < coverage.length; i++) {
                sb.append(String.valueOf(start++)).append("\t").append(coverage[i]).append("\n");
            }
        }

        int start = Globals.getINSTANCE().getALIGNMENT_BEGIN();
        int begin_H = -1;
        int end_H = -1;
        int begin_FH = -1;
        int end_FH = -1;
        int begin_T = -1;
        int end_T = -1;
        int begin_TT = -1;
        int end_TT = -1;
        for (int i = 0; i < coverage.length; i++) {
            if (begin_H == -1 && coverage[i] >= 100) {
                begin_H = start + i;
            }
            if (begin_FH == -1 && coverage[i] >= 500) {
                begin_FH = start + i;
            }
            if (begin_T == -1 && coverage[i] >= 1000) {
                begin_T = start + i;
            }
            if (begin_TT == -1 && coverage[i] >= 10000) {
                begin_TT = start + i;
            }
        }
        for (int i = coverage.length - 1; i >= begin_H - start; i--) {
            if (end_H == -1 && coverage[i] >= 100) {
                end_H = start + i;
            }
            if (end_FH == -1 && coverage[i] >= 500) {
                end_FH = start + i;
            }
            if (end_T == -1 && coverage[i] >= 1000) {
                end_T = start + i;
            }
            if (end_TT == -1 && coverage[i] >= 10000) {
                end_TT = start + i;
            }
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "coverage.txt", sb.toString());
        if (Globals.getINSTANCE().isCOVERAGE()) {
            Utils.saveCoveragePlot();
            Globals.getINSTANCE().println("To create a coverage plot, please execute: R CMD BATCH support/coverage.R");
            if (begin_H == -1 || end_H == -1) {
                Globals.getINSTANCE().println("There is no region with a sufficient coverage of >100x");
            } else {
                Globals.getINSTANCE().println("A coverage >100x is in region " + begin_H + "-" + end_H + "");
                if (begin_FH == -1 || end_FH == -1) {
                    Globals.getINSTANCE().println("There is no region with a sufficient coverage of >500x");
                } else {
                    Globals.getINSTANCE().println("A coverage >500x is in region " + begin_FH + "-" + end_FH + "");
                    if (begin_T == -1 || end_T == -1) {
                        Globals.getINSTANCE().println("There is no region with a sufficient coverage of >1000x");
                    } else {
                        Globals.getINSTANCE().println("A coverage >1000x is in region " + begin_T + "-" + end_T + "");
                        if (begin_TT == -1 || end_TT == -1) {
                            Globals.getINSTANCE().println("There is no region with a sufficient coverage of >10000x");
                        } else {
                            Globals.getINSTANCE().println("A coverage >10000x is in region " + begin_TT + "-" + end_TT + "");
                        }
                    }
                }
            }
            System.out.println("");
            System.out.println("Exiting.");
            System.exit(9);
        }
    }

    private static void printAlignment(Read[] reads) {
        if (Globals.getINSTANCE().isPRINT_ALIGNMENT()) {
            Globals.getINSTANCE().println("Saving alignment\t");
            new Summary().printAlignment(reads);
        }
    }

    private static void circos(int L, int[][] alignment) {
        if (Globals.getINSTANCE().isCIRCOS()) {
            new Summary().circos(L, alignment);
            System.exit(0);
        }
    }
}
