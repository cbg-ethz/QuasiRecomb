package ch.ethz.bsse.quasirecomb;

import ch.ethz.bsse.quasirecomb.distance.DistanceUtils;
import ch.ethz.bsse.quasirecomb.diversity.Diversity;
import ch.ethz.bsse.quasirecomb.diversity.PairwiseEntropyComparison;
import ch.ethz.bsse.quasirecomb.filter.Cutter;
import ch.ethz.bsse.quasirecomb.filter.MAExtract;
import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.ArtificialExperimentalForwarder;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelEntropy;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.quasiviz.QuasiViz;
import ch.ethz.bsse.quasirecomb.simulation.Recombinator;
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
 * Hello world! --train -i
 * C:/Users/XLR/Dropbox/QuasiRecomb/src/main/resources/haplotypes/dataset_3.fasta
 * -o C:/Users/XLR/Dropbox/QuasiRecomb/src/main/resources/d3Std -f
 * .28,.28,.26,.03,.03,.03,.03,.03,.03 -K 3 -verbose
 *
 * --trainHarder -i
 * C:/Users/XLR/Dropbox/QuasiRecomb/src/main/resources/d3Std/optimumJavaK3
 * -verbose -o C:/Users/XLR/Dropbox/QuasiRecomb/src/main/resources/d3Harder/
 */
public class Startup {

    @Option(name = "-i")
    private String input;
    @Option(name = "--recombine")
    private boolean recombine;
    @Option(name = "-spots")
    private String spots;
    @Option(name = "--sample", usage = "Sample from given trained model", metaVar = "OPTIMUMJAVA", multiValued = true)
    private boolean sample;
    @Option(name = "-o", usage = "Path to the output directory (default: current directory)", metaVar = "PATH")
    private String output;
    @Option(name = "--muentropy")
    private boolean muentropy;
    @Option(name = "--train", usage = "Train model for given multiple alignment")
    private boolean train;
    @Option(name = "--summary")
    private boolean summary;
//    @Option(name = "-noSample", usage = "Do not infer haplotypes from best model")
//    private boolean noSample;
    @Option(name = "-c10")
    private boolean crossvalidation = false;
    @Option(name = "-b")
    private boolean bootstrap = false;
    @Option(name = "-K")
    private String K = "1:5";
    @Option(name = "-t")
    private int t = 50;
    @Option(name = "-a")
    private double a = 0.01;
    @Option(name = "-N")
    private int N = 2000;
    @Option(name = "-chunk")
    private int chunk = 20;
    @Option(name = "-maxsim")
    private int maxsim = 1;
    @Option(name = "-f")
    private String f;
    @Option(name = "-e")
    private double e = 0.0003;
    @Option(name = "-d")
    private double d = 1e-8;
    @Option(name = "-min")
    private double minLLH = Double.NEGATIVE_INFINITY;
    @Option(name = "-parallelRestarts")
    private boolean parallelRestarts;
    @Option(name = "-singleCore")
    private boolean singleCore;
    @Option(name = "-verbose")
    private boolean verbose;
    @Option(name = "-noRecomb")
    private boolean noRecomb;
    @Option(name = "-alphah")
//    private double alphah = 0.0001;
    private double alphah = 0.01;
    @Option(name = "-betah")
    private double betah = 2;
    @Option(name = "-alphaz")
//    private double alphaz = 0.0015;
    private double alphaz = 0.01;
    @Option(name = "-betaz")
//    private double betaz = 0.0025;
    private double betaz = 0.005;
    @Option(name = "--filter")
    private boolean filter;
    @Option(name = "-c")
    private double gapc = 0.01;
    @Option(name = "--viz")
    private boolean viz;
    @Option(name = "-exec")
    private String exec;
    @Option(name = "--cut")
    private boolean cut;
    @Option(name = "-begin")
    private int begin;
    @Option(name = "-end")
    private int end;
    @Option(name = "-entropy")
    private boolean entropy;
    @Option(name = "-plot")
    private boolean plot;
    @Option(name = "-shannonindex")
    private boolean shannonindex;
    @Option(name = "--diversity")
    private boolean diversity;
    @Option(name = "-desc")
    private String description;
    @Option(name = "--pairEntropyTest")
    private boolean pairEntropyTest;
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
    private int SAMPLING_AMOUNT = 10000;

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
            } else if (this.output.endsWith(File.separator)
                    && !new File(this.output).exists()) {
                new File(this.output).mkdirs();
            }
            Globals.DEBUG = this.verbose;

            if (this.sample) {
                if (input.contains("#")) {
                    String[] splitBracket = input.split("#");
                    String[] split = splitBracket[1].split("-");

                    for (int i = Integer.parseInt(split[0]); i <= Integer.parseInt(split[1]); i++) {
                        System.out.println("Sampling " + splitBracket[0] + i);
                        ModelSampling simulation = new ModelSampling(splitBracket[0] + i, output, SAMPLING_AMOUNT);
                        simulation.save();
                    }
                } else {
                    System.out.println("Sampling " + input);
                    ModelSampling simulation = new ModelSampling(input, output, SAMPLING_AMOUNT);
                    simulation.save();
                }
            } else if (this.recombine) {
                if (this.spots != null) {
                    String[] split = this.spots.split(",");
                    int[] spots = new int[split.length];
                    int i = 0;
                    for (String s : split) {
                        spots[i++] = Integer.parseInt(s);
                    }
                    Recombinator.recombine(input, spots, output);
                } else {
                    System.out.println("Please provide -spots, i.e. -spots 50,140,321");
                }
            } else if (this.muentropy) {
                new ModelEntropy(this.input);
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

                Summary s = new Summary();
                System.out.println(s.print(or));
                if (this.haplotypes != null) {
                    ModelSampling ms = new ModelSampling(input, "", SAMPLING_AMOUNT);
                    System.out.println("\n#Quasispecies:");
                    ms.printQuasispecies();
                    Pair[] phi = DistanceUtils.calculatePhi(FastaParser.parseHaplotypeFile(haplotypes), ms.getReadsReversed());
                    System.out.println("\n#Phi distance:");
                    for (Pair p : phi) {
                        System.out.println(p.getValue0() + "\t" + p.getValue1());
                    }
                }
            } else if (this.distance) {
                Map<String, Double> quasiDouble = FastaParser.parseQuasispeciesFile(input);
                Map<String, Integer> quasiInt = new HashMap<>();
                for (String s : quasiDouble.keySet()) {
                    quasiInt.put(s, (int) (quasiDouble.get(s).doubleValue() * 10000));
                }
                Pair[] phi = DistanceUtils.calculatePhi(FastaParser.parseHaplotypeFile(haplotypes), quasiInt);
                System.out.println("\n#Phi distance:");
                for (Pair p : phi) {
                    System.out.println(p.getValue0() + "\t" + p.getValue1());
                }
            } else if (this.train) {
                Globals.CROSSVALIDATION = this.crossvalidation;
                Globals.BOOTSTRAP = this.bootstrap;
                int Kmin, Kmax;
                if (K.contains(":")) {
                    Kmin = Integer.parseInt(K.split(":")[0]);
                    Kmax = Integer.parseInt(K.split(":")[1]);
                } else {
                    Kmin = Integer.parseInt(K);
                    Kmax = Integer.parseInt(K);
                }
                double fArray[] = null;
                boolean exp = true;
                if (f != null) {
                    exp = false;
                    String[] split = f.split(",");
                    fArray = new double[split.length];
                    int i = 0;
                    double sum = 0d;
                    for (String s : split) {
                        fArray[i++] = Double.parseDouble(s);
                        sum += fArray[i - 1];
                    }
                    if (sum != 1d && Math.abs(sum) - 1d > 1e-6) {
                        System.err.println("Frequencies do not add up to 1, instead to " + sum);
                        System.exit(0);
                    }
                }

                Globals.ALPHA_Z = this.alphaz;
                Globals.ALPHA_H = this.alphah;
                Globals.BETA_Z = this.betaz;
                Globals.BETA_H = this.betah;
                Globals.MIN_LLH = this.minLLH;
                Globals.PARALLEL_JHMM = !this.singleCore;
                Globals.PARALLEL_RESTARTS = this.parallelRestarts;
                Globals.ESTIMATION_EPSILON = this.e;
//            Globals.SAMPLING_EPSILON = Globals.ESTIMATION_EPSILON;
                Globals.DELTA_LLH = this.d;
                Globals.REPEATS = this.t;
                Globals.DEBUG = this.verbose;
//            Globals.SIMULATION = config.getBoolean("Simulation");
//            Globals.DISTANCE = config.getBoolean("Distance");
                Globals.PARALLEL_RESTARTS_UPPER_BOUND = maxsim;
                Globals.STEPSIZE = chunk;
                Globals.savePath = output + File.separator;
                new File(Globals.savePath).mkdirs();
                if (this.noRecomb) {
                    Globals.rho0 = true;
                    Globals.rho0force = true;
                }
                ArtificialExperimentalForwarder.forward(exp, this.input, Kmin, Kmax, fArray, N);
            } else if (viz) {
                QuasiViz.paint(input, output);
            } else if (cut) {
                Cutter.cut(input, output, begin, end);
            } else if (diversity) {
                Diversity diversity = new Diversity(input, output, description, plot);
                if (entropy) {
                    diversity.entropy();
                }
                if (shannonindex) {
                    diversity.shannonIndex();
                }
            } else if (pairEntropyTest) {
//                File dir;
//                FileFilter fileFilter;
//                if (input.contains(File.separator)) {
//                    String[] split = input.split(File.separator);
//                    String path = "";
//                    for (int i = 0; i < split.length - 1; i++) {
//                        path += split[i] + File.separator;
//                    }
//                    dir = new File(path);
//                    fileFilter = new WildcardFileFilter(split[split.length - 1].replace("#", "*"));
//                } else {
//                    dir = new File(".");
//                    fileFilter = new WildcardFileFilter(input.replace("#", "*"));
//                }
//                File[] files = dir.listFiles(fileFilter);
//                for (int i = 0; i < files.length; i++) {
//                    System.out.println(files[i]);
//                }
                new PairwiseEntropyComparison(afile, bfile);
            } else if (filter) {
                if (input.contains("#")) {
                    String[] splitBracket = input.split("#");
                    String[] split = splitBracket[1].split("-");

                    for (int i = Integer.parseInt(split[0]); i <= Integer.parseInt(split[1]); i++) {
                        System.out.println("Filtering " + splitBracket[0] + i);
                        MAExtract.calc(splitBracket[0] + i, output + i, gapc);
                    }
                } else {
                    System.out.println("Filtering " + input);
                    MAExtract.calc(input, output, gapc);
                }
            } else {
                throw new CmdLineException("No input given");
            }

        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            System.err.println("java -jar QuasiRecomb.jar [options...] arguments...\n");
            System.err.println(" ------------------------");
            System.err.println(" === GENERAL options ===");
            System.err.println("  -o PATH\t\t: Path to the output directory (default: current directory)");
            System.err.println("  -i PATH\t\t: Path to the input file [REQUIRED]\n\t\t\t More information about the input type in sections below");
            System.err.println("  -verbose\t\t: Print debug information");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" ------------------------");
            System.err.println(" === TRAINING options ===");
            System.err.println("  --train");
            System.err.println("  -i INPUT\t\t: Multiple fasta file");
            System.err.println("");
            System.err.println("  -K INT or INT:INT\t: The interval or fixed number of sequence generators, i.e. 1:4 or 2\n\t\t\t  In a grid enviroment the $SGE_TASK_ID");
            System.err.println("  -t INT\t\t: The number of EM restarts to find optimum (default: 50)");
            System.err.println("  -e DOUBLE\t\t: Error rate of the sequencing machine (default: 0.0001)");
            System.err.println("  -min DOUBLE\t\t: Minimal likelihood which has to be reached after 50 iterations");
            System.err.println("  -d DOUBLE\t\t: Relative likehood change cut off (default: 1e-8)");
            System.err.println("  -f DOUBLE,DOUBLE,...\t: The distribution of the original haplotypes\n"
                    + "\t\t\t  Comma seperated format, i.e. 0.8,0.1,0.1\n"
                    + "\t\t\t  If not specified, the input is treated as experimental dataset");
            System.err.println("  -parallelRestarts\t: Parallelize the EM restarts, use this only on machines with 10+ cores!");
            System.err.println("");
            System.err.println("  Example for training:\n   java -jar QuasiRecomb.jar --train -i input.fasta");
//            System.err.println("  -noSample\t\t: Do not infer haplotypes, only model training");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" ------------------------");
            System.err.println(" === FILTER alignment ===");
            System.err.println("  --filter");
            System.err.println("  -i INPUT\t\t: Multiple alignment in fasta format");
            System.err.println("");
            System.err.println("  -c DOUBLE\t\t: Percentage of gaps allowed (default: 0.01)");
            System.err.println("");
            System.err.println("  Example for filtering:\n   java -jar QuasiRecomb.jar --filter -i input.fasta -o output_filtered.fasta");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" ------------------------");
            System.err.println(" === Cut alignment ===");
            System.err.println("  --cut");
            System.err.println("  -i INPUT\t\t: Multiple alignment in fasta format");
            System.err.println("");
            System.err.println("  -begin INT\t\t: Beginning position of the window");
            System.err.println("  -end INT\t\t: Ending position of the window");
            System.err.println("");
            System.err.println("  Example for cutting:\n   java -jar QuasiRecomb.jar --cut -i input.fasta -o output_w1200-2420.fasta -begin 1200 -end 2420");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" ------------------------");
            System.err.println(" === Visualize Quasispezies ===");
            System.err.println("  --viz");
            System.err.println("  -i INPUT\t\t: Multiple alignment in fasta format");
            System.err.println("");
            System.err.println("  -exec PATH\t\t: Executable dot file of graphviz");
//            System.err.println("  -c DOUBLE\t\t\t: Percentage of gaps allowed (default: 0.01)");
            System.err.println("");
            System.err.println("  Example for vizualization:\n   java -jar QuasiRecomb.jar --viz -i hapDist.fasta -o hapViz -exec /PATH/TO/dot");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" ------------------------");
            System.err.println(" === SAMPLE from model === ");
            System.err.println("  --sample ");
            System.err.println("  -i FILE\t\t: Sample from given trained model");
            System.err.println("");
            System.err.println("  Example for sampling:\n   java -jar QuasiRecomb.jar --sample -i path/to/optimumJava");
            System.err.println("");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" ------------------------");
            System.err.println(" === DIVERSITY measurements === ");
            System.err.println("  --diversity ");
            System.err.println("  -i FILE\t\t: Multiple fasta file");
            System.err.println("  -desc STRING\t\t: Short description (i.e. name of strain)");
            System.err.println("  -shannonindex\t\t\t: Calculate Shannon-index");
            System.err.println("  -entropy\t\t\t: Plot site-wise Shannon-Entropy");
            System.err.println("  -coverage\t\t\t: Plot site-wise Coverage");
            System.err.println("");
            System.err.println("  Example for entropy:\n   java -jar QuasiRecomb.jar --entropy -i input.far");
            System.err.println(" ------------------------");
            System.err.println("");
            System.err.println(" ------------------------");
            System.err.println(" === DISTANCE (phi) === ");
            System.err.println("  --distance ");
            System.err.println("  -i FILE\t\t: Multiple fasta file with quasispecies incl. frequencies");
            System.err.println("  -h FILE\t\t: Multiple fasta file with original haplotypes sampled from");
            System.err.println("");
            System.err.println("  Example for distance:\n   java -jar QuasiRecomb.jar --distance -i quasiespecies.fasta -h dataset.fasta");
            System.err.println(" ------------------------");
        }
    }
}
