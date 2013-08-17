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
package ch.ethz.bsse.quasirecomb;

import ch.ethz.bsse.quasirecomb.distance.DistanceUtils;
import ch.ethz.bsse.quasirecomb.distance.HammerWorker;
import ch.ethz.bsse.quasirecomb.distance.IntersectQuasispecies;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.Preprocessing;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.simulation.Recombinator;
import ch.ethz.bsse.quasirecomb.simulation.Simulator;
import ch.ethz.bsse.quasirecomb.utils.CutNHam;
import ch.ethz.bsse.quasirecomb.utils.Cutter;
import ch.ethz.bsse.quasirecomb.utils.FastaParser;
import ch.ethz.bsse.quasirecomb.utils.StatusUpdate;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import com.google.common.collect.Collections2;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import net.sf.samtools.SAMFormatException;
import org.javatuples.Pair;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Startup {

    public static void main(String[] args) throws IOException {
        new Startup().doMain(args);
        System.exit(0);
    }
    //GENERAL
    @Option(name = "-i")
    private String input;
    @Option(name = "-print")
    private boolean print;
    @Option(name = "-o", usage = "Path to the output directory (default: current directory)", metaVar = "PATH")
    private String output;
    @Option(name = "-log")
    private boolean log;
    @Option(name = "-verbose")
    private boolean verbose;
    //TRAIN
    @Option(name = "-K")
    private String K = "1-5";
    @Option(name = "-prune")
    private boolean prune;
    @Option(name = "-multMu")
    private double multMu = 0;
    @Option(name = "-multMuMin")
    private double multMuMin = 0;
    @Option(name = "-multRhoMin")
    private double multRhoMin = 0;
    @Option(name = "-multRho")
    private double multRho = 0;
    @Option(name = "-nosample")
    private boolean nosample;
    @Option(name = "-m")
    private int m = 5;
    @Option(name = "-t")
    private int t = 50;
    @Option(name = "-N")
    private int N = 2000;
    @Option(name = "-e")
    private double e = .001;
    @Option(name = "-noInfoEps")
    private boolean noInfoEps;
    @Option(name = "-d")
    private double d = 1e-4;
    @Option(name = "-dd")
    private double dd = 1e-8;
    @Option(name = "-alphah")
    private double alphah = 1e-4;
    @Option(name = "-alphaz")
    private double alphaz = 1e-4;
    @Option(name = "-p")
    private double p = 1e-10;
    @Option(name = "-noRecomb")
    private boolean noRecomb;
    @Option(name = "-betaz")
    private double betaz = 0.1;
    @Option(name = "-logBic")
    private boolean logBIC;
    @Option(name = "--recombine")
    private boolean recombine;
    @Option(name = "-spots")
    private String spots;
    @Option(name = "--sample", usage = "Sample from given trained model", metaVar = "OPTIMUMJAVA", multiValued = true)
    private boolean sample;
    @Option(name = "-cutoff")
    private double cutoff = 0;
    @Option(name = "--summary")
    private boolean summary;
    @Option(name = "--cut")
    private boolean cut;
    @Option(name = "-begin")
    private int begin;
    @Option(name = "-end")
    private int end;
    @Option(name = "-size")
    private int size;
    @Option(name = "--hamming")
    private boolean hamming;
    @Option(name = "--distance")
    private boolean distance;
    @Option(name = "--distanceDetail")
    private boolean distanceDetail;
    @Option(name = "-h")
    private String haplotypes;
    @Option(name = "--simulate")
    private boolean simulate;
    @Option(name = "-paired")
    private boolean paired;
    @Option(name = "-f")
    private String f;
    @Option(name = "-L")
    private int L;
    @Option(name = "-snapshots")
    private boolean snapshots;
    @Option(name = "-minmem")
    private boolean minmem;
    @Option(name = "-plot")
    private boolean plot;
    @Option(name = "-debug")
    private boolean debug;
    @Option(name = "--html")
    private boolean html;
    @Option(name = "-pdelta")
    private boolean pdelta;
    @Option(name = "-overlap")
    private boolean overlap;
    @Option(name = "-optimum")
    private String optimum;
    @Option(name = "--kl")
    private boolean kl;
    @Option(name = "--intersect")
    private boolean intersect;
    @Option(name = "--circos")
    private boolean circos;
    @Option(name = "-g")
    private String genome;
    @Option(name = "-printAlignment")
    private boolean printAlignment;
    @Option(name = "-amplicons")
    private String amplicons;
    @Option(name = "-ampliconDist")
    private String ampliconDist;
    @Option(name = "-length")
    private int length;
    @Option(name = "-muPrior")
    private boolean muPrior;
    @Option(name = "-stopQuick")
    private boolean stopQuick;
    @Option(name = "-spikeRho")
    private boolean spikeRho = true;
    @Option(name = "-conservative")
    private boolean conservative;
    @Option(name = "-global")
    private boolean global;
    @Option(name = "-silent")
    private boolean silent;
    @Option(name = "-unpaired")
    private boolean unpaired;
    @Option(name = "-steps")
    private int steps = 2;
    @Option(name = "-interpolateMu")
    private double interpolateMu = 1;
    @Option(name = "-interpolateRho")
    private double interpolateRho = 1;
    @Option(name = "-refine")
    private boolean refine;
    @Option(name = "-r")
    private String region;
    @Option(name = "-quality")
    private boolean quality;
    @Option(name = "--cutnham")
    private boolean cutnham;
    @Option(name = "--annotate")
    private boolean annotate;
    @Option(name = "-sampleReads")
    private boolean sampleReads;
    @Option(name = "-sampleProteins")
    private boolean sampleProteins;
    @Option(name = "-bootstrap")
    private boolean bootstrap;
    @Option(name = "-coverage")
    private boolean coverage;
    @Option(name = "-max")
    private boolean max;
    @Option(name = "--extended")
    private boolean extended;
    @Option(name = "-noGaps")
    private boolean noGaps;
    @Option(name = "-annealing")
    private boolean annealing;
    @Option(name = "-noGradient")
    private boolean noGradient;
    @Option(name = "-maxDel")
    private double maxDel = Integer.MAX_VALUE;
    @Option(name = "-maxPercDel")
    private double maxPercDel = 1;
    @Option(name = "-HIV")
    private String hiv;
    @Option(name = "-prior")
    private String prior;
    @Option(name = "--mix")
    private boolean mix;
    @Option(name = "-onlyPaired")
    private boolean onlyPaired;

    private void setInputOutput() {
        if (output == null) {
            this.output = System.getProperty("user.dir") + File.separator;
        } else {
            Globals.getINSTANCE().setSAVEPATH(this.output);
        }
        if (output.endsWith("/") || output.endsWith("\\")) {
            if (!new File(this.output).exists()) {
                if (!new File(this.output).mkdirs()) {
                    System.out.println("Cannot create directory: " + this.output);
                }
            }
        }
    }

    private void setMainParameters() {
        Globals.getINSTANCE().setSILENT(this.silent);
        Globals.getINSTANCE().setSTORAGE(!this.minmem);
        Globals.getINSTANCE().setSNAPSHOTS(this.snapshots);
        Globals.getINSTANCE().setDEBUG(this.verbose || this.debug);
        Globals.getINSTANCE().setPRINT(this.print || this.debug);
        Globals.getINSTANCE().setLOGGING(this.log);
        Globals.getINSTANCE().setPAIRED(this.paired);
        Globals.getINSTANCE().setPLOT(this.plot);
        Globals.getINSTANCE().setCUTOFF(this.cutoff);
        Globals.getINSTANCE().setSTEPS(this.steps);
        Globals.getINSTANCE().setSAMPLE_READS(this.sampleReads);
        Globals.getINSTANCE().setSAMPLE_PROTEINS(this.sampleProteins);
        Globals.getINSTANCE().setCOVERAGE(this.coverage);
        Globals.getINSTANCE().setMAX_DEL(this.maxDel);
        Globals.getINSTANCE().setMAX_OVERALL_DEL(this.maxPercDel);
        Globals.getINSTANCE().setNO_GAPS(this.noGaps);
    }

    private void sample() {
        Globals.getINSTANCE().setNREAL(N);
        ModelSampling simulation = new ModelSampling(input, output);
        simulation.save();
    }

    private void simulate() throws RuntimeException, NumberFormatException {
        Globals.getINSTANCE().setOVERLAP(this.overlap);
        double fArray[];
        String[] split = f.split(",");
        fArray = new double[split.length];
        int i = 0;
        double sum = 0d;
        for (String s : split) {
            fArray[i++] = Double.parseDouble(s);
            sum += fArray[i - 1];
        }
        if (sum != 1d && Math.abs(sum - 1d) > 1e-6) {
            throw new RuntimeException("Frequencies do not add up to 1, instead to " + sum);
        }
        if (this.output.endsWith(File.separator)) {
            this.output += "reads";
        }
        if (this.global) {
            Simulator.fromHaplotypesGlobal(FastaParser.parseFarFile(input), N, L, this.e, fArray, this.output);
        } else if (paired) {
            Simulator.fromHaplotypesGlobalPaired(FastaParser.parseFarFile(input), N, L, this.e, fArray, this.output);
        } else {
            if (this.amplicons != null) {
                Simulator.fromHaplotypesGlobalAmplicon(FastaParser.parseFarFile(input), N, L, this.e, fArray, this.output, this.amplicons, this.ampliconDist, this.length);
            } else {
                Simulator.fromHaplotypes(FastaParser.parseFarFile(input), N, L, this.e, fArray, 5, this.output);
            }
        }
    }

    private void recombine() throws NumberFormatException {
        if (this.spots != null) {
            String[] split = this.spots.split(",");
            int[] spotsArray = new int[split.length];
            int i = 0;
            for (String s : split) {
                spotsArray[i++] = Integer.parseInt(s);
            }
            Recombinator.recombine(input, spotsArray, output);
        } else {
            System.out.println("Please provide -spots, i.e. -spots 50,140,321");
        }
    }

    private void cutnham() {
        new CutNHam(input, begin, end, size);
    }

    private void hamming() {
        Map<String, String> hapMap = FastaParser.parseHaplotypeFile(input);
        for (String hap : hapMap.keySet()) {
            hapMap.put(hap, hapMap.get(hap).substring(1));
        }
        String[] fasta = hapMap.keySet().toArray(new String[hapMap.size()]);
        Globals.getINSTANCE().setHammingMax((int) Math.pow(fasta.length, 2));
        Map<String, Map<String, Integer>> invoke = Globals.getINSTANCE().getFjPool().invoke(new HammerWorker(fasta, 0, fasta.length));
        StatusUpdate.getINSTANCE().println("done");
        int i = 0;
        Utils.saveFile(output + "dist.txt", "");
        StringBuilder sb = new StringBuilder();
        sb.append("PAIRS");
        for (String f : fasta) {
            sb.append("\t").append(hapMap.get(f));
        }
        sb.append("\n");
        for (String f : fasta) {
            Map<String, Integer> entry = invoke.get(f);
            sb.append(hapMap.get(f));
            for (String f2 : fasta) {
                Integer dist = entry.get(f2);
                sb.append("\t ");
                sb.append(String.valueOf(dist));
            }
            sb.append("\n");
            StatusUpdate.getINSTANCE().print("SB\t" + i++);
        }
        Utils.saveFile(output + "dist.txt", sb.toString());
        System.out.println(sb.toString());
    }

    private void intersect() {
        if (N == 2000) {
            N = 2;
        }
        Map<String, Double> map = new ModelSampling(input, output).getMap();
        for (int i = 1; i < N; i++) {
            Map<String, Double> map2 = new ModelSampling(input, output).getMap();
            map = IntersectQuasispecies.compute(map, map2);
        }

        StringBuilder sb = new StringBuilder();
        int i = 0;
        for (Object o : ModelSampling.sortMapByValue(map).entrySet()) {
            Entry<String, Double> e = (Entry) o;
            sb.append(">INTERSECT").append(i++).append("_").append(IntersectQuasispecies.shorten(e.getValue())).append("\n");
            sb.append(e.getKey()).append("\n");
        }
        System.out.println(sb.toString());
    }

    private void kullbackLeibler() {
        if (this.input.equals("+")) {
            List<String> fileList = new ArrayList<>(Arrays.asList(new File(System.getProperty("user.dir") + File.separator).list()));
            fileList.remove("support");
            String[] files = fileList.toArray(new String[fileList.size()]);
            Arrays.sort(files);
            boolean startHeader = true;
            for (int i = 0; i < files.length; i++) {
                if (files[i].equals("support")) {
                    continue;
                }
                if (!startHeader) {
                    System.out.print("\t");
                } else {
                    startHeader = false;
                }
                System.out.print(files[i].split("\\.")[0]);
            }
            System.out.println("");
            for (int i = 0; i < files.length; i++) {
                if (files[i].equals("support")) {
                    continue;
                }
                boolean start = true;
                for (int j = 0; j < files.length; j++) {
                    if (files[i].equals("support")) {
                        continue;
                    }
                    Map<String, Double> x = FastaParser.parseQuasispeciesFile(files[i]);
                    Map<String, Double> x2 = FastaParser.parseQuasispeciesFile(files[j]);
                    Double calculateKLD2 = DistanceUtils.calculateKLD2(x, x2);
                    if (!start) {
                        System.out.print("\t");
                    } else {
                        start = false;
                    }
                    System.out.print(calculateKLD2);
                }
                System.out.println("");
            }
        } else {
            Map<String, Double> i = FastaParser.parseQuasispeciesFile(this.input);
            Map<String, Double> i2 = FastaParser.parseQuasispeciesFile(this.haplotypes);
            Double calculateKLD2 = DistanceUtils.calculateKLD2(i, i2);
            System.out.println("Jensen-Shannon divergence: " + calculateKLD2);
        }
    }

    private void html() {
        OptimalResult or = null;
        try {
            FileInputStream fis = new FileInputStream(input);
            try (ObjectInputStream in = new ObjectInputStream(fis)) {
                or = (OptimalResult) in.readObject();
            }
        } catch (IOException | ClassNotFoundException ex) {
            System.err.println(ex);
        }
        if (or != null) {
            Summary s = new Summary();
            System.out.println(s.html(or));
        }
    }

    private void summary() {
        OptimalResult or = null;
        try {
            FileInputStream fis = new FileInputStream(input);
            try (ObjectInputStream in = new ObjectInputStream(fis)) {
                or = (OptimalResult) in.readObject();
            }
        } catch (IOException | ClassNotFoundException ex) {
            System.err.println(ex);
        }
        if (or != null) {
            Summary s = new Summary();
            System.out.println(s.print(or));
        }
        if (this.haplotypes != null) {
            ModelSampling ms = new ModelSampling(input, "");
            System.out.println("\n#Quasispecies:");
            ms.printQuasispecies();
            Pair[] phi = DistanceUtils.calculatePhi(FastaParser.parseHaplotypeFile(haplotypes), ms.getReadsReversed());
            System.out.println("\n#Phi distance:");
            System.out.println("q\tphi");
            int i = 0;
            for (Pair phiLocal : phi) {
                System.out.println(i++ + "\t" + phiLocal.getValue0());
                if (((double) phiLocal.getValue0()) == 1d) {
                    break;
                }
            }
        }
    }

    private void distance() {
        Map<String, Double> quasiDouble = FastaParser.parseQuasispeciesFile(input);
        Map<String, String> haps = FastaParser.parseHaplotypeFile(haplotypes);
        double[] precision = DistanceUtils.calculatePhi2(haps, quasiDouble);

        String[] quasispecies = quasiDouble.keySet().toArray(new String[quasiDouble.size()]);
        quasiDouble.clear();
        Map<String, Pair<String, Double>> quasiHeadMap = new HashMap<>();
        for (String h : haps.keySet()) {
            quasiHeadMap.put(h, Pair.with(haps.get(h), 1d / haps.size()));
        }
        haps.clear();
        for (String q : quasispecies) {
            haps.put(q, "");
        }
        String[] recall = DistanceUtils.calculatePhi3(haps, quasiHeadMap);
        System.out.println("q\tprecision\trecall");
        int i = 0;
        int max = Math.max(precision.length, recall.length);
        for (int j = 0; j < max; j++) {
            System.out.println(j + "\t" + (j < precision.length ? precision[j] : "1.0") + "\t" + (j < recall.length ? recall[j] : "1.0"));
        }
    }

    private void distanceDetail() {
        Map<String, Double> quasiDouble = FastaParser.parseQuasispeciesFile(input);
        Map<String, String> haps = FastaParser.parseHaplotypeFile(haplotypes);
        Map<String, String> strain2Hap = new HashMap<>();
        double[][] precision = new double[haps.size()][300];
        String[] head = new String[haps.size()];
        int x = 0;
        for (Entry<String, String> entry : haps.entrySet()) {
            head[x++] = entry.getValue().replaceAll(">", "");
            strain2Hap.put(entry.getValue().replaceAll(">", ""), entry.getKey());
        }
        Arrays.sort(head);
        x = 0;
        for (String strain : head) {
            Map<String, String> h = new HashMap<>();
            h.put(strain2Hap.get(strain), strain);
            precision[x++] = DistanceUtils.calculatePhi2(h, quasiDouble);
        }
        for (int i = 0; i < head.length; i++) {
            System.out.print("\t" + head[i]);
        }
        System.out.println("\tFP");
        for (int j = 0; j < 200; j++) {
            System.out.print(j);
            double sum = 0d;
            for (int i = 0; i < precision.length; i++) {
                System.out.print("\t" + precision[i][j]);
                sum += precision[i][j];
            }
            System.out.println("\t" + (1d - sum));
        }
    }

    private void cut() {
        String[] split = this.region.split(";");
        for (String s : split) {
            String[] r = s.split("-");
            int start = Integer.parseInt(r[0]);
            int stop = Integer.parseInt(r[1]);
            Cutter.cut(input, start + "-" + stop + ".fasta", start - 1, stop - 1);
        }
    }

    private void circos() {
//        for (int i = 0; i < 10; i++) {
//
//        }

        Globals.getINSTANCE().setCIRCOS(this.circos);
        Globals.getINSTANCE().setGENOME(this.genome);

        if (this.hiv != null) {
            int begin_hiv = -1;
            int end_hiv = -1;
            switch (this.hiv) {
                case "1":
                case "p17":
                    begin_hiv = 790;
                    end_hiv = 1186;
                    break;
                case "2":
                case "p24":
                    begin_hiv = 1186;
                    end_hiv = 1879;
                    break;
                case "3":
                case "p2p6":
                    begin_hiv = 1879;
                    end_hiv = 2292;
                    break;
                case "4":
                case "prot":
                    begin_hiv = 2253;
                    end_hiv = 2550;
                    break;
                case "5":
                case "RT":
                    begin_hiv = 2550;
                    end_hiv = 3870;
                    break;
                case "6":
                case "RNase":
                    begin_hiv = 3870;
                    end_hiv = 4230;
                    break;
                case "7":
                case "int":
                    begin_hiv = 4230;
                    end_hiv = 5096;
                    break;
                case "8":
                case "vif":
                    begin_hiv = 5041;
                    end_hiv = 5619;
                    break;
                case "9":
                case "vpr":
                    begin_hiv = 5559;
                    end_hiv = 5850;
                    break;
                case "10":
                case "vpu":
                    begin_hiv = 6062;
                    end_hiv = 6310;
                    break;
                case "11":
                case "gp120":
                    begin_hiv = 6225;
                    end_hiv = 7758;
                    break;
                case "12":
                case "gp41":
                    begin_hiv = 7758;
                    end_hiv = 8795;
                    break;
                case "13":
                case "nef":
                    begin_hiv = 8797;
                    end_hiv = 9417;
                    break;
                case "14":
                case "gag":
                    begin_hiv = 790;
                    end_hiv = 2292;
                    break;
                case "15":
                case "pol":
                    begin_hiv = 2085;
                    end_hiv = 5096;
                    break;
                case "16":
                case "env":
                    begin_hiv = 6225;
                    end_hiv = 8795;
                    break;
                case "17":
                case "cg":
                    begin_hiv = 490;
                    end_hiv = 9540;
                    break;
                case "LTRGAG":
                    begin_hiv = 490;
                    end_hiv = 2292;
                    break;
                case "5LTR":
                    begin_hiv = 490;
                    end_hiv = 790;
                    break;
                case "ENVGP120":
                    begin_hiv = 5041;
                    end_hiv = 7758;
                    break;
                case "GP41LTR":
                    begin_hiv = 7758;
                    end_hiv = 9540;
                    break;
                case "VVV":
                    begin_hiv = 5096;
                    end_hiv = 6310;
                    break;
                case "POL":
                    begin_hiv = 2253;
                    end_hiv = 5096;
                    break;
            }
            Globals.getINSTANCE().setWINDOW_BEGIN(begin_hiv - 1);
            Globals.getINSTANCE().setWINDOW_END(end_hiv - 1);
            Globals.getINSTANCE().setWINDOW(true);
        }
        Preprocessing.workflow(this.input, 0, 0);
    }

    private void train() throws NumberFormatException, CmdLineException {
        if (this.input == null) {
//            System.out.println("No input given");
//            System.exit(0);
            throw new CmdLineException("No input given");
        }
        int Kmin, Kmax;
        if (K.contains("-")) {
            Kmin = Integer.parseInt(K.split("-")[0]);
            Kmax = Integer.parseInt(K.split("-")[1]);
        } else {
            Kmin = Integer.parseInt(K);
            Kmax = Integer.parseInt(K);
        }

        if (this.hiv != null) {
            int begin_hiv = -1;
            int end_hiv = -1;
            switch (this.hiv) {
                case "1":
                case "p17":
                    begin_hiv = 790;
                    end_hiv = 1186;
                    break;
                case "2":
                case "p24":
                    begin_hiv = 1186;
                    end_hiv = 1879;
                    break;
                case "3":
                case "p2p6":
                    begin_hiv = 1879;
                    end_hiv = 2292;
                    break;
                case "4":
                case "prot":
                    begin_hiv = 2253;
                    end_hiv = 2550;
                    break;
                case "5":
                case "RT":
                    begin_hiv = 2550;
                    end_hiv = 3870;
                    break;
                case "6":
                case "RNase":
                    begin_hiv = 3870;
                    end_hiv = 4230;
                    break;
                case "7":
                case "int":
                    begin_hiv = 4230;
                    end_hiv = 5096;
                    break;
                case "8":
                case "vif":
                    begin_hiv = 5041;
                    end_hiv = 5619;
                    break;
                case "9":
                case "vpr":
                    begin_hiv = 5559;
                    end_hiv = 5850;
                    break;
                case "10":
                case "vpu":
                    begin_hiv = 6062;
                    end_hiv = 6310;
                    break;
                case "11":
                case "gp120":
                    begin_hiv = 6225;
                    end_hiv = 7758;
                    break;
                case "12":
                case "gp41":
                    begin_hiv = 7758;
                    end_hiv = 8795;
                    break;
                case "13":
                case "nef":
                    begin_hiv = 8797;
                    end_hiv = 9417;
                    break;
                case "14":
                case "gag":
                    begin_hiv = 790;
                    end_hiv = 2292;
                    break;
                case "15":
                case "pol":
                    begin_hiv = 2085;
                    end_hiv = 5096;
                    break;
                case "16":
                case "env":
                    begin_hiv = 6225;
                    end_hiv = 8795;
                    break;
            }
            Globals.getINSTANCE().setWINDOW_BEGIN(begin_hiv - 1);
            Globals.getINSTANCE().setWINDOW_END(end_hiv - 1);
            Globals.getINSTANCE().setWINDOW(true);
        } else if (this.region != null && !this.region.isEmpty()) {
            String[] r = this.region.split("-");
            Globals.getINSTANCE().setWINDOW_BEGIN(Integer.parseInt(r[0]) - 1);
            Globals.getINSTANCE().setWINDOW_END(Integer.parseInt(r[1]) - 1);
            Globals.getINSTANCE().setWINDOW(true);
        }

        Globals.getINSTANCE().setGRADIENT(!this.noGradient);
        if (!this.noGradient) {
            Globals.getINSTANCE().setMULT_MU(100);
            Globals.getINSTANCE().setMULT_RHO(1000);
            Globals.getINSTANCE().setMULT_RHO_MIN(100);
            Globals.getINSTANCE().setMULT_MU_MIN(10);
        }
        if (this.multRhoMin > 0) {
            Globals.getINSTANCE().setMULT_RHO_MIN(this.multRhoMin);
        }
        if (this.multMuMin > 0) {
            Globals.getINSTANCE().setMULT_MU_MIN(this.multMuMin);
        }
        if (this.conservative) {
            Globals.getINSTANCE().setALPHA_H(1e-6);
            Globals.getINSTANCE().setALPHA_Z(1e-6);
            Globals.getINSTANCE().setINTERPOLATE_MU(1);
            Globals.getINSTANCE().setINTERPOLATE_RHO(1);
            Globals.getINSTANCE().setMULT_RHO_MIN(1);
            Globals.getINSTANCE().setMULT_MU_MIN(1);
        } else {
            Globals.getINSTANCE().setALPHA_H(this.alphah);
            Globals.getINSTANCE().setALPHA_Z(this.alphaz);
            if (this.multMu > 0) {
                Globals.getINSTANCE().setMULT_MU(this.multMu);
            }
            if (this.multRho > 0) {
                Globals.getINSTANCE().setMULT_RHO(this.multRho);
            }
            Globals.getINSTANCE().setINTERPOLATE_MU(this.interpolateMu);
            Globals.getINSTANCE().setINTERPOLATE_RHO(this.interpolateRho);
        }

        Globals.getINSTANCE().setNO_QUALITY(!this.quality);
        Globals.getINSTANCE().setSPIKERHO(this.spikeRho);
        Globals.getINSTANCE().setUNPAIRED(this.unpaired);
        Globals.getINSTANCE().setSTOP_QUICK(this.stopQuick);
        Globals.getINSTANCE().setPRINT_ALIGNMENT(this.printAlignment);
        Globals.getINSTANCE().setPRIORMU(this.muPrior);
        Globals.getINSTANCE().setNOSAMPLE(this.nosample);
        Globals.getINSTANCE().setPDELTA(this.pdelta);
        Globals.getINSTANCE().setPRUNE(this.prune);
        Globals.getINSTANCE().setPCHANGE(this.p);
        Globals.getINSTANCE().setBETA_Z(this.betaz);
        Globals.getINSTANCE().setLOG_BIC(this.logBIC);
        if (this.e != .001) {
            Globals.getINSTANCE().setFLAT_EPSILON_PRIOR(true);
        }
        Globals.getINSTANCE().setUNINFORMATIVE_EPSILON_PRIOR(this.noInfoEps);
        Globals.getINSTANCE().setESTIMATION_EPSILON(this.e);
        Globals.getINSTANCE().setDELTA_LLH(this.d);
        Globals.getINSTANCE().setDELTA_REFINE_LLH(this.dd);
        Globals.getINSTANCE().setREPEATS(this.m);
        Globals.getINSTANCE().setDESIRED_REPEATS(this.t);
        Globals.getINSTANCE().setDEBUG(this.verbose);
        Globals.getINSTANCE().setSAVEPATH(output + File.separator);
        Globals.getINSTANCE().setNO_RECOMB(this.noRecomb);
        Globals.getINSTANCE().setFORCE_NO_RECOMB(this.noRecomb);
        Globals.getINSTANCE().setOPTIMUM(this.optimum);
        Globals.getINSTANCE().setUSER_OPTIMUM(this.optimum != null);
        if (this.refine) {
            Globals.getINSTANCE().setALPHA_H(1e-6);
            Globals.getINSTANCE().setALPHA_Z(1e-6);
            Globals.getINSTANCE().setMULT_MU(1);
            Globals.getINSTANCE().setMULT_RHO(1);
            Globals.getINSTANCE().setINTERPOLATE_MU(1);
            Globals.getINSTANCE().setINTERPOLATE_RHO(1);
            Globals.getINSTANCE().setOPTIMUM("support/best.optimum");
            Globals.getINSTANCE().setUSER_OPTIMUM(true);
            if (!new File(this.output + File.separator + "support/best.optimum").exists()) {
                System.err.println("QuasiRecomb needs to be executed once without -refine before it can be used with -refine.");
            }
        }
        Globals.getINSTANCE().setBOOTSTRAP(this.bootstrap);
        Globals.getINSTANCE().setMAX(this.max);
        Globals.getINSTANCE().setANNEALING(this.annealing);
        if (this.global) {
            System.err.println("");
            System.err.println("Parameter -global is not supported anymore.");
            System.err.println("Please use -conservative if the quasispecies should be peaked.");
            System.err.println("");
        }
        Globals.getINSTANCE().setK_MIN(Kmin);
        Globals.getINSTANCE().setPRIOR(this.prior);
        Globals.getINSTANCE().setONLY_PAIRED(this.onlyPaired);
        Preprocessing.workflow(this.input, Kmin, Kmax);
    }

    public void doMain(String[] args) throws IOException {
        CmdLineParser parser = new CmdLineParser(this);

        parser.setUsageWidth(80);
        try {
            parser.parseArgument(args);
            setInputOutput();
            StringBuilder sb = new StringBuilder();
            for (String arg : args) {
                sb.append(arg).append(" ");
            }
            new File(this.output + File.separator + "support/").mkdirs();
            Utils.appendFile(this.output + File.separator + "support/CMD", sb.toString());
            setMainParameters();

            if (this.sample) {
                sample();
            } else if (this.simulate) {
                simulate();
            } else if (this.recombine) {
                recombine();
            } else if (this.hamming) {
                hamming();
            } else if (this.intersect) {
                intersect();
            } else if (this.kl) {
                kullbackLeibler();
            } else if (this.html) {
                html();
            } else if (this.summary) {
                summary();
            } else if (this.distance) {
                distance();
            } else if (this.distanceDetail) {
                distanceDetail();
            } else if (this.circos) {
                circos();
            } else if (cut) {
                cut();
            } else if (cutnham) {
                cutnham();
            } else if (annotate) {
                annotate();
            } else if (mix) {
                mix();
            } else {
                train();
            }

        } catch (SAMFormatException e) {
            System.err.println("");
            System.err.println("Input file is not in SAM or BAM format.");
            System.err.println(e);
        } catch (CmdLineException cmderror) {
            System.err.println(cmderror.getMessage());
            System.err.println("");
            System.err.println("QuasiRecomb version: " + Startup.class.getPackage().getImplementationVersion());
            System.err.println("Get latest version from http://bit.ly/QuasiRecomb");
            System.err.println("");
            System.err.println("USAGE: java -jar QuasiRecomb.jar options...\n");
            System.err.println(" -------------------------");
            System.err.println(" === GENERAL options ===");
            System.err.println("  -i INPUT\t\t: Alignment file in BAM or SAM format.");
            System.err.println("  -o PATH\t\t: Path to the output directory (default: current directory).");
            System.err.println("");
            System.err.println("  -K INT or INT-INT\t: The interval or fixed number of sequence generators, i.e. 1-4 or 2\n\t\t\t  In a grid enviroment the $SGE_TASK_ID."
                    + "\n\t\t\t  In case of no input, K will be incremented as long as max BIC has not been reached, but will stop at K=5.");
            System.err.println("  -m INT\t\t: The number of EM restarts during model selection (default: 5).");
            System.err.println("  -t INT\t\t: The number of EM restarts for best K to find optimum (default: 50).");
            System.err.println("  -r INT-INT\t\t: Only reconstruct a specific region.");
            System.err.println("  -noRecomb\t\t: Do not allow recombination.");
            System.err.println("  -quality\t\t: Account phred quality scores (slower runtime).");
            System.err.println("  -printAlignment\t: Save alignment.txt in a human readable format.");
//            System.err.println("  -sampleReads\t\t: Sample reads in addition to haplotypes");
            System.err.println("  -sampleProteins\t: Sample full-length protein sequences in three reading frames.");
            System.err.println("  -coverage\t\t: If your dataset only contains a single region of interest, "
                    + "\n\t\t\t  regions with a minimum coverage of 100x, 500x, 1,000x and 10,000x are reported.");
            System.err.println("  -bootstrap\t\t: Model-selection is performed on 10 bootstrapped datasets. Very time consuming, but robust.");
            System.err.println("  -refine\t\t: Can only be used after QuasiRecomb has been executed once before on the same dataset in the same directory."
                    + "\n\t\t\t  Thins the number of haplotypes.");
            System.err.println("  -noGaps\t\t: Ignore gaps; useful if data is 454 and gaps are only technical errors.");
            System.err.println("  -conservative\t\t: Use this if the major haplotypes are only of interest.");
            System.err.println("  -maxDel INT\t\t: Remove reads with more than INT consecutive deletions.");
            System.err.println("  -maxPercDel DOUBLE\t: Remove reads with more than DOUBLE ratio of deletions, between 0.0 - 1.0");
            System.err.println("  -unpaired\t\t: If read names are not unique and reads are single-end, prevent pairing and merging.");
            System.err.println(" -------------------------");
            System.err.println(" === Technical options ===");
            System.err.println("  -XX:NewRatio=9\t: Reduces the memory consumption (RECOMMENDED to use).");
            System.err.println("  -Xms2G -Xmx10G\t: Increase heap space.");
            System.err.println("  -XX:+UseParallelGC\t: Enhances performance on multicore systems.");
            System.err.println("  -XX:+UseNUMA\t\t: Enhances performance on multi-CPU systems.");
            System.err.println(" -------------------------");
            System.err.println(" === EXAMPLES ===");
            System.err.println("   java -XX:NewRatio=9 -jar QuasiRecomb.jar -i alignment.bam");
            System.err.println("   java -XX:NewRatio=9 -jar QuasiRecomb.jar -i alignment.bam -conservative ");
            System.err.println("   java -XX:NewRatio=9 -jar QuasiRecomb.jar -i alignment.bam -K 1:10");
            System.err.println("   java -XX:NewRatio=9 -jar QuasiRecomb.jar -i alignment.bam -noRecomb -r 790-2292");
            System.err.println("   java -XX:+UseParallelGC -Xms2g -Xmx10g -XX:+UseNUMA -XX:NewRatio=9 -jar QuasiRecomb.jar -i alignment.bam");
            System.err.println(" -------------------------");
            System.err.println("  For further information, see http://bit.ly/QuasiRecomb-howto");
            System.err.println(" -------------------------");
//            System.err.println("  -d DOUBLE\t\t: Relative likehood threshold (default: 1e-8)");
//            System.err.println("  -pdelta\t\t: Stop if there is no change of parameters, convergence criterium");
//            System.err.println("  -e DOUBLE\t\t: Fix error rate of the sequencing machine");
//            System.err.println("  -noInfoEps\t\t: Do not use the error rate of 0.8% as an informative prior");
//            System.err.println("");
//            System.err.println("");
            if (this.extended) {
                System.err.println(" === SAMPLE from model === ");
                System.err.println("  --sample \t\t: Sample from given trained model");
                System.err.println("  -i FILE\t\t: Path to best.optimum file");
                System.err.println("");
                System.err.println("  Example for sampling:\n   java -jar QuasiRecomb.jar --sample -i support/best.optimum");
//            System.err.println(" -------------------------");
//            System.err.println(" === DISTANCE === ");
//            System.err.println("  --distance ");
//            System.err.println("  -i FILE\t\t: Multiple fasta file with quasispecies incl. frequencies"
//                    + "\n\t\t\t  The corresponding frequency has to be the suffix in the fasta description delimited by an underscore, i.e. >seq1231_0.4212");
//            System.err.println("  -h FILE\t\t: Multiple fasta file with original haplotypes.");
//            System.err.println("");
//            System.err.println("  Example for distance:\n   java -jar QuasiRecomb.jar --distance -i quasispecies.fasta -h dataset.fasta");
                System.err.println(" -------------------------");
                System.err.println(" === DISTANCE === ");
                System.err.println("  --distanceDetail\t: Reports frequencies of the original haplotypes that are present in the quasispecies and false-positive rate, allowing q mismatches");
                System.err.println("  -i FILE\t\t: Multiple fasta file with quasispecies incl. frequencies"
                        + "\n\t\t\t  The corresponding frequency has to be the suffix in the fasta description delimited by an underscore, i.e. >seq1231_0.4212");
                System.err.println("  -h FILE\t\t: Multiple fasta file with original haplotypes.");
                System.err.println("");
                System.err.println("  Example for distance:\n   java -jar QuasiRecomb.jar --distanceDetail -i quasispecies.fasta -h dataset.fasta");
                System.err.println(" -------------------------");
//            System.err.println(" === SIMULATE === ");
//            System.err.println("  --simulate ");
//            System.err.println("  -i FILE\t\t: Multiple fasta file with haplotypes");
//            System.err.println("  -f INT_ARRAY\t\t: Array with frequencies for haplotypes."
//                    + "\n\t\t\t  The number of frequencies has to be equal the number of haplotypes.");
//            System.err.println("  -e DOUBLE\t\t: Error-rate per base, per position (default: 0.003)");
//            System.err.println("  -N int\t\t: Number of reads");
//            System.err.println("  -paired\t\t: Paired-end reads 2x250bp");
//            System.err.println("");
//            System.err.println("  Example for distance:\n   java -jar QuasiRecomb.jar --simulate -i quasispecies.fasta -h dataset.fasta");
//            System.err.println(" -------------------------");
            }
        }
    }

    private void posterior(OptimalResult or, String haplotypeFile) {
        Map<String, Double> haps = FastaParser.parseQuasispeciesFile(haplotypeFile);

    }

//    private void viterbi(OptimalResult or, String haplotypeFile) {
//        int K = or.getK();
//        int L = or.getL();
//
//        double[] pi = new double[K];
//        double piSum = 0;
//        for (int j = 0; j < or.getPi().length; j++) {
//            for (int k = 0; k < or.getPi()[j].length; k++) {
//                pi[k] += or.getPi()[j][k];
//                piSum += or.getPi()[j][k];
//            }
//        }
//        for (int k = 0; k < K; k++) {
//            pi[k] /= piSum;
//        }
//
//        Map<String, Double> haps = FastaParser.parseQuasispeciesFile(haplotypeFile);
//
//        for (String hap : haps.keySet()) {
//            byte[] read = Utils.splitReadIntoByteArray(hap);
//            double[][] f = new double[K][L];
//            double[][] f2 = new double[K][L];
//            for (int k = 0; k < K; k++) {
//                f[k][0] = pi[k] * or.getMu()[0][k][read[0]];
//                f2[k][0] = 0;
//            }
//            for (int j = 1; j < L; j++) {
//                for (int k = 0; k < K; k++) {
//                    double max = -1;
//                    int argmax = -1;
//                    for (int l = 0; l < K; l++) {
//                        double x = f[l][j - 1] * or.getRho()[j - 1][l][k] * or.getMu()[j - 1][k][read[j]];
//                        if (x > max) {
//                            max = x;
//                            argmax = l;
//                        }
//                    }
//                    f[k][j] = max;
//                    f2[k][j] = argmax;
//                }
//            }
//            double[] z = new double[L];
//            double[] x = new double[L];
//            for (int k = 0; k < K; k++) {
//                if (f[k][L - 1] > z[L - 1]) {
//                    z[L - 1] = k;
//                }
//            }
//            x[L-1] =
////            z[L-1] =
//        }
//    }
    private void annotate() {
        StringBuilder fasta = new StringBuilder();
        StringBuilder protein = new StringBuilder();
        Map<String, Double> parseQuasispeciesFile = FastaParser.parseQuasispeciesFile(this.input);
        Map<String, Double> mutationMap = new HashMap<>();
        Map<String, Double> proteins = new LinkedHashMap<>();
        Map<String, String> proteinAnnotation = new LinkedHashMap<>();
        int i = 0;
        for (Object o : ModelSampling.sortMapByValue(parseQuasispeciesFile).entrySet()) {
            Entry<String, Double> e = (Entry<String, Double>) o;
            StringBuilder sb = new StringBuilder();
            String c = e.getKey();
            String s = ModelSampling.dna2protein(c, 0);
            char[] cs = s.toCharArray();
            StringBuilder mutations = new StringBuilder();
            int secondary = 0;
            if (cs[21] == 'M') {
                secondary++;
                mutations.append("M");
            } else {
                mutations.append(" ");
            }
            if (cs[44] == 'A') {
                secondary++;
                mutations.append("A");
            } else {
                mutations.append(" ");
            }
            if (cs[90] == 'H') {
                secondary++;
                mutations.append("H");
            } else {
                mutations.append(" ");
            }
            if (cs[98] == 'I') {
                secondary++;
                mutations.append("I");
            } else {
                mutations.append(" ");
            }
            if (cs[110] == 'R') {
                secondary++;
                mutations.append("R");
            } else {
                mutations.append(" ");
            }
            if ((cs[90] == 'R' || cs[90] == 'C') && cs[102] == 'H') {
                sb.append(17 + secondary);
                mutations.append(cs[90]);
                mutations.append(cs[102]);
            } else if (cs[90] == 'R' || cs[90] == 'C') {
                sb.append(6 + secondary);
                mutations.append(cs[90]);
                mutations.append(" ");
            } else if (cs[102] == 'H') {
                sb.append(11 + secondary);
                mutations.append(" ");
                mutations.append(cs[102]);
            } else {
                mutations.append("  ");
                sb.append(secondary);
            }
            String mutationSummary = mutations.toString();
            if (mutationMap.containsKey(mutationSummary)) {
                mutationMap.put(mutationSummary, mutationMap.get(mutationSummary) + e.getValue());
            } else {
                mutationMap.put(mutationSummary, e.getValue());
            }

            fasta.append(">read").append(String.valueOf(i)).append("_").append(parseQuasispeciesFile.get(c)).append("_0_").append(sb).append("\n").append(c).append("\n");
            if (proteins.containsKey(s)) {
                proteins.put(s, proteins.get(s) + e.getValue());
            } else {
                proteins.put(s, e.getValue());
                proteinAnnotation.put(s, sb.toString());
            }
            i++;
        }
        StringBuilder mutationSummary = new StringBuilder();
        int maxLength = 0;
        for (Object o : ModelSampling.sortMapByValue(mutationMap).entrySet()) {
            Entry<String, Double> e = (Entry<String, Double>) o;
            maxLength = Math.max(maxLength, e.getKey().length());
        }
        for (Object o : ModelSampling.sortMapByValue(mutationMap).entrySet()) {
            Entry<String, Double> e = (Entry<String, Double>) o;
            mutationSummary.append(e.getKey());
            for (int j = e.getKey().length() - 1; j < maxLength; j++) {
                mutationSummary.append(" ");
            }
            mutationSummary.append(Math.round(e.getValue() * 10000) / 10000d).append("\n");
        }
        i = 0;
        for (Object o : ModelSampling.sortMapByValue(proteins).entrySet()) {
            Entry<String, Double> e = (Entry<String, Double>) o;
            protein.append(">read").append(String.valueOf(i++)).append("_").append(e.getValue()).append("_0_").append(proteinAnnotation.get(e.getKey())).append("\n").append(e.getKey()).append("\n");
        }
        Utils.saveFile(this.input + "_annotated", fasta.toString());
        Utils.saveFile(this.input + "_protein0_annotated", protein.toString());
        Utils.saveFile(this.input + "_summary", mutationSummary.toString());

        char[] as = {' ', 'M'};
        char[] bs = {' ', 'A'};
        char[] cs = {' ', 'H'};
        char[] ds = {' ', 'I'};
        char[] es = {' ', 'R'};
        char[] fs = {' ', 'H'};
        char[] gs = {' ', 'R', 'C'};
        StringBuilder table = new StringBuilder();
        StringBuilder table2 = new StringBuilder();
        for (char a : as) {
            for (char b : bs) {
                for (char c : cs) {
                    for (char d : ds) {
                        for (char e : es) {
                            for (char f : fs) {
                                for (char g : gs) {
                                    String x = "" + a + b + c + d + e + g + f;
                                    table2.append(x).append("\n");
                                    table.append(mutationMap.get(x) == null ? "0.0" : Math.round(mutationMap.get(x) * 10000) / 10000d).append("\n");
                                }
                            }
                        }
                    }
                }
            }
        }
        Utils.saveFile(this.input + "_tableOverview", table2.toString());
        Utils.saveFile(this.input + "_table", table.toString());
    }
//    private void annotate() {
//        StringBuilder fasta = new StringBuilder();
//        StringBuilder hesamFasta = new StringBuilder();
//        Map<String, Double> parseQuasispeciesFile = FastaParser.parseQuasispeciesFile(this.input);
//
//        int i = 0;
//        for (Object o : ModelSampling.sortMapByValue(parseQuasispeciesFile).keySet()) {
//            StringBuilder hesam = new StringBuilder();
//            StringBuilder sb = new StringBuilder();
//            String s = (String) o;
//            char[] cs = s.toCharArray();
//            int secondary = 0;
//            if (cs[21] == 'M') {
//                secondary++;
//                hesam.append("1 ");
//            } else {
//                hesam.append("0 ");
//            }
//            if (cs[44] == 'A') {
//                secondary++;
//                hesam.append("1 ");
//            } else {
//                hesam.append("0 ");
//            }
//            if (cs[98] == 'I') {
//                secondary++;
//                hesam.append("1 ");
//            } else {
//                hesam.append("0 ");
//            }
//            if ((cs[90] == 'R' || cs[90] == 'C') && cs[102] == 'H') {
//                sb.append(12 + secondary);
//            } else if (cs[90] == 'R' || cs[90] == 'C') {
//                sb.append(4 + secondary);
//            } else if (cs[102] == 'H') {
//                sb.append(8 + secondary);
//            } else {
//                sb.append(secondary);
//            }
//
//            if (cs[90] == 'R' || cs[90] == 'C') {
//                hesam.append("1 ");
//            } else {
//                hesam.append("0 ");
//            }
//            if (cs[102] == 'H') {
//                hesam.append("1");
//            } else {
//                hesam.append("0");
//            }
//
////            sb.append(cs[21]);
////            sb.append(cs[44]);
////            sb.append(cs[90]);
////            sb.append(cs[98]);
////            sb.append(cs[102]);
//            fasta.append(">read").append(String.valueOf(i)).append("_").append(parseQuasispeciesFile.get(s)).append("_0_").append(sb).append("\n").append(s).append("\n");
//            hesamFasta.append("read").append(String.valueOf(i)).append("_").append(parseQuasispeciesFile.get(s)).append(" ").append(hesam).append("\n");
//            i++;
//        }
//        Utils.saveFile(this.input + "_annotated", fasta.toString());
//        Utils.saveFile(this.input + "_hesam", hesamFasta.toString());
//    }

    private void mix() {
        String[] f1 = FastaParser.parseFarFile(this.input);
        String[] f2 = FastaParser.parseFarFile(this.haplotypes);
        if (f1.length != f2.length) {
            System.err.println("Files are of different sizes.");
            System.exit(0);
        }
        List<Integer> order = new LinkedList<>();
        for (int i = 0; i < f2.length; i++) {
            order.add(i);
        }
        int x = 0;
        Collection<List<Integer>> permutations = Collections2.permutations(order);
        for (List<Integer> l : permutations) {
            StringBuilder sb = new StringBuilder();
            int i = 0;
            for (int j : l) {
                sb.append(">").append(i).append(j).append("\n");
                sb.append(f1[i++]).append(f2[j]).append("\n");
            }
            Utils.saveFile("permutation_" + x++ + ".fasta", sb.toString());
        }
    }
}
