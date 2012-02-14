package ch.ethz.bsse.quasirecomb;

import ch.ethz.bsse.quasirecomb.entropy.ShannonEntropy;
import ch.ethz.bsse.quasirecomb.filter.Cutter;
import ch.ethz.bsse.quasirecomb.filter.MAExtract;
import ch.ethz.bsse.quasirecomb.model.ArtificialExperimentalForwarder;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.quasiviz.QuasiViz;
import java.io.File;
import java.io.IOException;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * Hello world!
 *
 */
public class Startup {

    @Option(name = "-i")
    private String input;
    @Option(name = "--sample", usage = "Sample from given trained model", metaVar = "OPTIMUMJAVA", multiValued = true)
    private boolean sample;
    @Option(name = "-o", usage = "Path to the output directory (default: current directory)", metaVar = "PATH")
    private String output;
    @Option(name = "--train", usage = "Train model for given multiple alignment")
    private boolean train;
    @Option(name = "-noSample", usage = "Do not infer haplotypes from best model")
    private boolean noSample;
    @Option(name = "-K")
    private String K = "1:5";
    @Option(name = "-t")
    private int t = 50;
    @Option(name = "-N")
    private int N = 2000;
    @Option(name = "-chunk")
    private int chunk = 20;
    @Option(name = "-maxsim")
    private int maxsim = 1;
    @Option(name = "-f")
    private String f;
    @Option(name = "-e")
    private double e = 0.0001;
    @Option(name = "-d")
    private double d = 1e-8;
    @Option(name = "-min")
    private double minLLH = Double.NEGATIVE_INFINITY;
    @Option(name = "-parallelRestarts")
    private boolean parallelRestarts;
    @Option(name = "-verbose")
    private boolean verbose;
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
    @Option(name = "--entropy")
    private boolean entropy;
    

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
            }

            if (this.sample) {
                int AMOUNT = 10000;
                if (input.contains("#")) {
                    String[] splitBracket = input.split("#");
                    String[] split = splitBracket[1].split("-");

                    for (int i = Integer.parseInt(split[0]); i <= Integer.parseInt(split[1]); i++) {
                        System.out.println("Sampling " + splitBracket[0] + i);
                        ModelSampling simulation = new ModelSampling(splitBracket[0] + i, output, AMOUNT);
                    }
                } else {
                    System.out.println("Sampling " + input);
                    ModelSampling simulation = new ModelSampling(input, output, AMOUNT);
                }
            } else if (this.train) {
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
                Globals.MIN_LLH = this.minLLH;

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
                if (Globals.PRIOR_ALPHA == 0.0) {
                    Globals.rho0 = true;
                    Globals.rho0force = true;
                }
                ArtificialExperimentalForwarder.forward(exp, this.input, Kmin, Kmax, fArray, N);
            } else if (viz) {
                QuasiViz.paint(input, output);
            } else if (cut) {
                Cutter.cut(input, output, begin, end);
            } else if (entropy) {
                ShannonEntropy.entropy(input);
            } else if (filter) {
                if (input.contains("#")) {
                    String[] splitBracket = input.split("#");
                    String[] split = splitBracket[1].split("-");

                    for (int i = Integer.parseInt(split[0]); i <= Integer.parseInt(split[1]); i++) {
                        System.out.println("Filtering " + splitBracket[0] + i);
                        MAExtract.calc(splitBracket[0]+i, output+i, gapc);
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
            System.err.println(" === SHANNON ENTROPY === ");
            System.err.println("  --entropy ");
            System.err.println("  -i FILE\t\t: Multiple fasta file");
            System.err.println("");
            System.err.println("  Example for entropy:\n   java -jar QuasiRecomb.jar --entropy -i input.far");
            System.err.println("");
        }
    }
}
