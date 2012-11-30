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
import ch.ethz.bsse.quasirecomb.utils.Cutter;
import ch.ethz.bsse.quasirecomb.utils.FastaParser;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import org.javatuples.Pair;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Startup {

    //GENERAL
    @Option(name = "-i")
    private String input;
    @Option(name = "-i2")
    private String input2;
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
    private String K = "1:5";
    @Option(name = "-prune")
    private boolean prune;
    @Option(name = "-spikeRho")
    private boolean spikeRho;
    @Option(name = "-multMu")
    private double multMu = 100;
    @Option(name = "-multRho")
    private double multRho = 100;
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
    private double d = 1e-6;
    @Option(name = "-dd")
    private double dd = 1e-8;
    @Option(name = "-alphah")
    private double alphah = 1e-3;
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
    @Option(name = "--summary")
    private boolean summary;
    @Option(name = "--cut")
    private boolean cut;
    @Option(name = "-begin")
    private int begin;
    @Option(name = "-end")
    private int end;
    @Option(name = "--hamming")
    private boolean hamming;
    @Option(name = "--distance")
    private boolean distance;
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
    @Option(name = "-perturb")
    private int perturb = 0;
    @Option(name = "--kl")
    private boolean kl;
    @Option(name = "--intersect")
    private boolean intersect;
    @Option(name = "-circos")
    private boolean circos;

    public static void main(String[] args) throws IOException {
        new Startup().doMain(args);
        System.exit(0);
    }

    public void doMain(String[] args) throws IOException {
        CmdLineParser parser = new CmdLineParser(this);

        parser.setUsageWidth(80);
        try {
            parser.parseArgument(args);
            setInputOutput();
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
            } else if (cut) {
                cut();
            } else {
                train();
            }

        } catch (CmdLineException cmderror) {
            System.err.println(cmderror.getMessage());
            System.err.println("java -jar QuasiRecomb.jar options...\n");
            System.err.println(" ------------------------");
            System.err.println("  -i INPUT\t\t: Multiple fasta file");
            System.err.println("  -o PATH\t\t: Path to the output directory (default: current directory)");
            System.err.println("  -paired\t\t: Reads are paired");
            System.err.println("");
            System.err.println("  -K INT or INT:INT\t: The interval or fixed number of sequence generators, i.e. 1:4 or 2\n\t\t\t  In a grid enviroment the $SGE_TASK_ID."
                    + "\n\t\t\t  In case of no input, K will be incremented as long as max BIC has not been reached, but will stop at K=5.");
            System.err.println("  -m INT\t\t: The number of EM restarts during model selection (default: 5)");
            System.err.println("  -t INT\t\t: The number of EM restarts for best K to find optimum (default: 50)");
            System.err.println("  -d DOUBLE\t\t: Relative likehood threshold (default: 1e-8)");
            System.err.println("  -pdelta\t\t: Stop if no parameters change, convergence criterium");
            System.err.println("  -e DOUBLE\t\t: Fix error rate of the sequencing machine");
            System.err.println("  -noInfoEps\t\t: Do not use the error rate of 0.8% as an informative prior");
            System.err.println("  -noRecomb\t\t: Do not allow recombination");
            System.err.println("  -plot\t\t: Plot coverage");
            System.err.println("");
            System.err.println("  Example for training:\n   java -jar QuasiRecomb.jar -i input.fasta");
            System.err.println(" ------------------------");
//            System.err.println("");
//            System.err.println("");
//            System.err.println(" ------------------------");
//            System.err.println(" === SAMPLE from model === ");
//            System.err.println("  --sample ");
//            System.err.println("  -i FILE\t\t: Sample from given trained model");
//            System.err.println("");
//            System.err.println("  Example for sampling:\n   java -jar QuasiRecomb.jar --sample -i path/to/optimumJava");
//            System.err.println("");
//            System.err.println(" ------------------------");
//            System.err.println("");
//            System.err.println(" ------------------------");
//            System.err.println(" === SUMMARY of model === ");
//            System.err.println("  --sample ");
//            System.err.println("  -i FILE\t\t: Summary of given trained model");
//            System.err.println("  -h FILE\t\t: Calculates phi distance to this true haplotypes");
//            System.err.println("");
//            System.err.println("  Example for summary:\n   java -jar QuasiRecomb.jar --summary -i path/to/optimumJava");
//            System.err.println("");
//            System.err.println(" ------------------------");
//            System.err.println("");
//            System.err.println(" ------------------------");
//            System.err.println(" === DISTANCE (phi) === ");
//            System.err.println("  --distance ");
//            System.err.println("  -i FILE\t\t: Multiple fasta file with quasispecies incl. frequencies\n\t\t\t  The corresponding frequency has to be the suffix in the fasta description delimited by an underscore, i.e. >seq1231_0.4212");
//            System.err.println("  -h FILE\t\t: Multiple fasta file with original haplotypes sampled from");
//            System.err.println("");
//            System.err.println("  Example for distance:\n   java -jar QuasiRecomb.jar --distance -i quasiespecies.fasta -h dataset.fasta");
//            System.err.println(" ------------------------");
        }
    }

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
        Globals.getINSTANCE().setSTORAGE(!this.minmem);
        Globals.getINSTANCE().setSNAPSHOTS(this.snapshots);
        Globals.getINSTANCE().setDEBUG(this.verbose || this.debug);
        Globals.getINSTANCE().setPRINT(this.print || this.debug);
        Globals.getINSTANCE().setLOGGING(this.log);
        Globals.getINSTANCE().setPAIRED(this.paired);
        Globals.getINSTANCE().setPLOT(this.plot);
    }

    private void sample() {
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
        if (paired) {
            Simulator.fromHaplotypesGlobalPaired(FastaParser.parseFarFile(input), N, L, this.e, fArray, this.output);
        } else {
            Simulator.fromHaplotypes(FastaParser.parseFarFile(input), N, L, this.e, fArray, 4, this.output);
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

    private void hamming() {
        Map<String, String> hapMap = FastaParser.parseHaplotypeFile(input);
        for (String hap : hapMap.keySet()) {
            hapMap.put(hap, hapMap.get(hap).substring(1));
        }
        String[] fasta = hapMap.keySet().toArray(new String[hapMap.size()]);
        Globals.getINSTANCE().setHammingMax((int) Math.pow(fasta.length, 2));
        Map<String, Map<String, Integer>> invoke = Globals.getINSTANCE().getFjPool().invoke(new HammerWorker(fasta, 0, fasta.length));
        Globals.getINSTANCE().println("done");
        int i = 0;
        Utils.saveFile(output + "dist.txt", "");
        StringBuilder sb = new StringBuilder();
        sb.append("PAIRS");
        for (Entry<String, Map<String, Integer>> entry : invoke.entrySet()) {
            sb.append("\t").append(hapMap.get(entry.getKey()));
        }
        sb.append("\n");
        for (Entry<String, Map<String, Integer>> entry : invoke.entrySet()) {
            sb.append(hapMap.get(entry.getKey()));
            for (Entry<String, Integer> inner : entry.getValue().entrySet()) {
                sb.append("\t ");
                sb.append(String.valueOf(inner.getValue()));
            }
            sb.append("\n");
            Globals.getINSTANCE().print("SB\t" + i++);
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
            System.out.println(s.kl(or));
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

    private void cut() {
        Cutter.cut(input, output, begin, end);
    }

    private void train() throws NumberFormatException, CmdLineException {
        if (this.input == null) {
            throw new CmdLineException("No input given");
        }
        int Kmin, Kmax;
        if (K.contains(":")) {
            Kmin = Integer.parseInt(K.split(":")[0]);
            Kmax = Integer.parseInt(K.split(":")[1]);
        } else {
            Kmin = Integer.parseInt(K);
            Kmax = Integer.parseInt(K);
        }

        Globals.getINSTANCE().setMULT_MU(this.multMu);
        Globals.getINSTANCE().setMULT_RHO(this.multRho);
        Globals.getINSTANCE().setCIRCOS(this.circos);
        Globals.getINSTANCE().setNOSAMPLE(this.nosample);
        Globals.getINSTANCE().setALPHA_H(this.alphah);
        Globals.getINSTANCE().setPDELTA(this.pdelta);
        Globals.getINSTANCE().setPRUNE(this.prune);
        Globals.getINSTANCE().setPERTURB(this.perturb);
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
        Preprocessing.workflow(this.input, Kmin, Kmax);
    }
}
