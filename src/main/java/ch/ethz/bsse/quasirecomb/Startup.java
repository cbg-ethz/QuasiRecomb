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

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.model.Preprocessing;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.simulation.Recombinator;
import ch.ethz.bsse.quasirecomb.simulation.Sampling;
import ch.ethz.bsse.quasirecomb.utils.Cutter;
import ch.ethz.bsse.quasirecomb.utils.DistanceUtils;
import ch.ethz.bsse.quasirecomb.utils.FastaParser;
import ch.ethz.bsse.quasirecomb.utils.Summary;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.HashMap;
import java.util.Map;
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
    @Option(name = "-m")
    private int m = 5;
    @Option(name = "-t")
    private int t = 50;
    @Option(name = "-N")
    private int N = 2000;
    @Option(name = "-samplingNumber")
    private int samplingNumber = 10000;
    @Option(name = "-e")
    private double e = .001;
    @Option(name = "-ee")
    private double ee = .001;
    @Option(name = "-d")
    private double d = 1e-8;
    @Option(name = "-p")
    private double p = 1e-4;
    @Option(name = "-parallelRestarts")
    private boolean parallelRestarts;
    @Option(name = "-singleCore")
    private boolean singleCore;
    @Option(name = "-flatEpsPrior")
    private boolean flatEpsilonPrior;
    @Option(name = "-noRecomb")
    private boolean noRecomb;
    @Option(name = "-alphah")
    private double alphah = 0.00001;
    @Option(name = "-alphaz")
    private double alphaz = 0.001;
    @Option(name = "-betaz")
    private double betaz = 0.01;
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
    @Option(name = "-afile")
    private String afile;
    @Option(name = "-bfile")
    private String bfile;
    @Option(name = "-h")
    private String haplotypes;
    @Option(name = "--simulate")
    private boolean simulate;
    @Option(name = "-f")
    private String f;
    @Option(name = "-L")
    private int L;

    public static void main(String[] args) throws IOException {
        new Startup().doMain(args);
    }

    public void doMain(String[] args) throws IOException {
        CmdLineParser parser = new CmdLineParser(this);

        parser.setUsageWidth(80);
        try {
            parser.parseArgument(args);
            if (output == null) {
                this.output = System.getProperty("user.dir") + File.separator;
            } else {
                Globals.getINSTANCE().setSAVEPATH(this.output);
            }
            if (!new File(this.output).exists()) {
                if (!new File(this.output).mkdirs()) {
                    System.out.println("Cannot create directory: " + this.output);
                }
            }

            Globals.getINSTANCE().setDEBUG(this.verbose);
            Globals.getINSTANCE().setLOGGING(this.log);
            Globals.getINSTANCE().setLOG_BIC(this.logBIC);
            Globals.getINSTANCE().setSAMPLING_NUMBER(this.samplingNumber);
            Globals.getINSTANCE().setPRINT(this.print);

            if (this.sample) {
                ModelSampling simulation = new ModelSampling(input, output);
                simulation.save();
            } else if (this.simulate) {
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
                    this.output += "reads.fasta";
                }
                Sampling.fromHaplotypes(FastaParser.parseFarFile(input), N, L, this.ee, fArray, 4, this.output);
//                Sampling.fromHaplotypesGlobal(FastaParser.parseFarFile(input), N, L, this.ee, fArray, 4, this.output);
            } else if (this.recombine) {
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
            } else if (this.hamming) {
                System.out.println(DistanceUtils.calcHamming(afile, bfile));
            } else if (this.summary) {
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
            } else if (this.distance) {
                Map<String, Double> quasiDouble = FastaParser.parseQuasispeciesFile(input);
                Map<String, Integer> quasiInt = new HashMap<>();
                for (Map.Entry<String, Double> s : quasiDouble.entrySet()) {
                    quasiInt.put(s.getKey(), (int) (s.getValue().doubleValue() * 10000));
                }
                Pair[] phi = DistanceUtils.calculatePhi(FastaParser.parseHaplotypeFile(haplotypes), quasiInt);
                System.out.println("q\tphi");
                int i = 0;
                for (Pair phiLocal : phi) {
                    System.out.println(i++ + "\t" + phiLocal.getValue0());
                    if (((double) phiLocal.getValue0()) == 1d) {
                        break;
                    }
                }

            } else if (cut) {
                Cutter.cut(input, output, begin, end);
            } else {
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

                Globals.getINSTANCE().setFLAT_EPSILON_PRIOR(this.flatEpsilonPrior);
                Globals.getINSTANCE().setPCHANGE(this.p);
                Globals.getINSTANCE().setALPHA_Z(this.alphaz);
                Globals.getINSTANCE().setALPHA_H(this.alphah);
                Globals.getINSTANCE().setBETA_Z(this.betaz);
                Globals.getINSTANCE().setPARALLEL_JHMM(!this.singleCore);
                Globals.getINSTANCE().setPARALLEL_RESTARTS(this.parallelRestarts);
                Globals.getINSTANCE().setESTIMATION_EPSILON(this.e);
                Globals.getINSTANCE().setSAMPLING_EPSILON(this.ee);
                Globals.getINSTANCE().setDELTA_LLH(this.d);
                Globals.getINSTANCE().setREPEATS(this.m);
                Globals.getINSTANCE().setDESIRED_REPEATS(this.t);
                Globals.getINSTANCE().setDEBUG(this.verbose);
                Globals.getINSTANCE().setSAVEPATH(output + File.separator);
                Globals.getINSTANCE().setNO_RECOMB(this.noRecomb);
                Globals.getINSTANCE().setFORCE_NO_RECOMB(this.noRecomb);
                Preprocessing.workflow(this.input, Kmin, Kmax, N);
            }

        } catch (CmdLineException cmderror) {
            System.err.println(cmderror.getMessage());
            System.err.println("java -jar QuasiRecomb.jar options...\n");
            System.err.println(" ------------------------");
//            System.err.println(" === GENERAL options ===");
//            System.err.println("  -i PATH\t\t: Path to the input file [REQUIRED]\n\t\t\t  More information about the input type in sections below");
//            System.err.println("  -verbose\t\t: Print debug information");
//            System.err.println(" ------------------------");
//            System.err.println("");
//            System.err.println(" ------------------------");
//            System.err.println(" === TRAINING options ===");
            System.err.println("  -i INPUT\t\t: Multiple fasta file");
            System.err.println("  -o PATH\t\t: Path to the output directory (default: current directory)");
            System.err.println("");
            System.err.println("  -K INT or INT:INT\t: The interval or fixed number of sequence generators, i.e. 1:4 or 2\n\t\t\t  In a grid enviroment the $SGE_TASK_ID."
                    + "\n\t\t\t  In case of no input, K will be incremented as long as max BIC has not been reached, but will stop at K=5.");
            System.err.println("  -m INT\t\t: The number of EM restarts during model selection (default: 5)");
            System.err.println("  -t INT\t\t: The number of EM restarts for best K to find optimum (default: 50)");
            System.err.println("  -d DOUBLE\t\t: Relative likehood threshold (default: 1e-8)");
            System.err.println("  -flatEpsPrior\t: Do not train epsilon, in combination with -e");
            System.err.println("  -e DOUBLE\t\t: Error rate of the sequencing machine (default: 0.0001)");
            System.err.println("  -noRecomb\t: Do not allow recombination");
            System.err.println("  -parallelRestarts\t: Parallelize the EM restarts, use this only on machines with 10+ cores!");
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
}
