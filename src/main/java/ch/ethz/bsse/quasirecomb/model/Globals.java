package ch.ethz.bsse.quasirecomb.model;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;

/** 
 * Information holder for all necessary given and inferred parameters.
 * 
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class Globals {

    public static boolean SAVE = true;
    public static String filePrefix = "";
    public static String savePath = "C:\\Users\\XLR\\Dropbox\\JData\\";
    public static double ESTIMATION_EPSILON = .0003;
    public static double SAMPLING_EPSILON = .0001;
    public static double FILTER_LLH = 1e-4;
    public static double DELTA_LLH = 1e-8;
    public static boolean MASK_RHO = false;
    public static int STEPSIZE = 20;
    public static final ForkJoinPool fjPool = new ForkJoinPool();
    public static boolean DEBUG = false;
    public static boolean SIMULATION = false;
    public static boolean DISTANCE = false;
    public static List<Integer> runtime = new LinkedList<>();
    public static boolean rho0 = false;
    public static int RUNS = 20;
    public static int REPEATS = 10;
    public static boolean rho0force = false;
    public static boolean TEST = false;
    public static int N;
    public static double MASK_RHO_THRESHOLD;
    public static double PRIOR_ALPHA = 0.01;
    public static boolean ADD_ALPHA = false;
    public static int MAX_PRE_BREAK = 50;
    public static int PARALLEL_RESTARTS_UPPER_BOUND = 10;
    public static boolean PARALLEL_JHMM = true;
    public static boolean PARALLEL_RESTARTS = false;
    public static boolean NO_REFINE = false;
    public static boolean NO_BREAK_THRESHOLD = false;
    public static double MIN_LLH = Double.NEGATIVE_INFINITY;
    public static double BIAS = Double.NEGATIVE_INFINITY;
    private static double MAX_LLH = Double.NEGATIVE_INFINITY;

    public static synchronized double getMAX_LLH() {
        return MAX_LLH;
    }

    public static synchronized void maxMAX_LLH(double llh) {
        MAX_LLH = Math.max(MAX_LLH, llh);
    }
}
