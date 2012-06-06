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
package ch.ethz.bsse.quasirecomb.model;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;

/**
 * Information holder for all necessary given and inferred parameters.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Globals {

    public static boolean FIX_EPSILON = false;
    public static double BETA_Z = 0.01;
    public static double BETA_H = 2;
    public static double ALPHA_Z = 0.001;
    public static boolean ALPHA_H_USER = false;
    public static double[][] ALPHA_H;
    public static double ALPHAH = 0.01;
    public static String[] HAPLOTYPE_ARRAY_EMPIRICAL;
    public static boolean GLOBAL = false;
    public static boolean BOOTSTRAP = false;
    public static boolean CROSSVALIDATION = false;
    public static boolean SAVE = true;
    public static String filePrefix = "";
    public static String savePath = "C:\\Users\\XLR\\Dropbox\\JData\\";
    public static double ESTIMATION_EPSILON = .0003;
    public static double SAMPLING_EPSILON = .001;
    public static double FILTER_LLH = 1e-4;
    public static double DELTA_LLH_HARDER = 1e-50;
    public static double DELTA_LLH = 1e-8;
    public static boolean MASK_RHO = false;
    public static int STEPSIZE = 50;
    public static final ForkJoinPool fjPool = new ForkJoinPool();
    public static boolean DEBUG = false;
    public static boolean SIMULATION = false;
    public static boolean DISTANCE = false;
    public static List<Integer> runtime = new LinkedList<>();
    public static boolean rho0 = false;
    public static int RUNS = 20;
    public static int REPEATS = 10;
    public static int DESIRED_REPEATS = 50;
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
    public static double BIAS = 0.1;
    private static double MAX_LLH = Double.NEGATIVE_INFINITY;
    public static StringBuilder LOG = new StringBuilder();
    public static boolean LOG_BIC = false;
    public static boolean LOGGING = false;
    public static boolean PRINT = true;
    public static int SAMPLING_NUMBER = 10000;

    public static void log(Object o) {
        if (PRINT) {
            System.out.print(o);
        } else {
            if (LOGGING) {
                LOG.append(o);
            }
        }
    }
    public static double PERCENTAGE = 0;

    public static void printPercentage(int K) {
        PERCENTAGE += 100d / Globals.REPEATS;
        System.out.print("\r\tK " + K + ":\t" + Math.round(PERCENTAGE*1000)/1000 + "%");
    }

    public static synchronized double getMAX_LLH() {
        return MAX_LLH;
    }

    public static synchronized void maxMAX_LLH(double llh) {
        MAX_LLH = Math.max(MAX_LLH, llh);
    }
}
