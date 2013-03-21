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
import com.google.common.collect.Multimaps;
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
            Multimap<Integer,Double> bics = ArrayListMultimap.create();
            Map<Read, Double> piMap = new HashMap<>();
            for (Read r : reads) {
                piMap.put(r, r.getCount() / N);
            }
            Frequency<Read> readDist = new Frequency<>(piMap);

            for (int i = 0; i < 10; i++) {
                Map<Integer, Read> hashed = new HashMap<>();

                for (int x = 0; x < N; x++) {
                    Read r = readDist.roll();
                    int hash = r.hashCode();
                    if (hashed.containsKey(hash)) {
                        hashed.get(hash).incCount();
                    } else {
                        hashed.put(hash, r);
                    }
                }
                
                Read[] rs = hashed.values().toArray(new Read[hashed.values().size()]);
                ModelSelection ms = new ModelSelection(rs, Kmin, Kmax, rs.length, L, n);
                bics.putAll(ms.getMsTemp().getMaxBICs());
            }
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
        if (Globals.getINSTANCE().isDEBUG()) {
            if (new File(Globals.getINSTANCE().getSAVEPATH() + "gap.txt").exists()) {
                new File(Globals.getINSTANCE().getSAVEPATH() + "gap.txt").delete();
            }
        }
        StringBuilder indelSB = new StringBuilder();
        for (Read r : reads) {
            r.shrink();
            indelSB.append(r.getInsertion()).append(" ");
            Globals.getINSTANCE().print("Modifying reads\t" + (Math.round((shrinkCounter++ / reads.length) * 100)) + "%");
        }
        if (Globals.getINSTANCE().isDEBUG()) {
            Utils.appendFile(Globals.getINSTANCE().getSAVEPATH() + "gap.txt", indelSB.toString());
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
        Globals.getINSTANCE().println("Alignment entropy\t" + alignmentEntropy);
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
        Globals.getINSTANCE().println("Plotting\t");
        StringBuilder sb = new StringBuilder();
        int start = Globals.getINSTANCE().getALIGNMENT_BEGIN();
        int[] coverage = Globals.getINSTANCE().getTAU_OMEGA().getCoverage();
        for (int i = 0; i < coverage.length; i++) {
            sb.append(String.valueOf(start++)).append("\t").append(coverage[i]).append("\n");
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "coverage.txt", sb.toString());
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
