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

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;

/**
 * Information holder for all necessary given and inferred parameters.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Globals {

    public static int ALIGNMENT_BEGIN = Integer.MAX_VALUE;
    public static int ALIGNMENT_END = Integer.MIN_VALUE;
    public static boolean FLAT_EPSILON_PRIOR;
    public static double PCHANGE;
    public static double BETA_Z;
    public static double ALPHA_Z;
    public static double ALPHA_H;
    public static String[] HAPLOTYPE_ARRAY_EMPIRICAL;
    public static String SAVEPATH;
    public static double ESTIMATION_EPSILON;
    public static double SAMPLING_EPSILON;
    public static double DELTA_LLH = 1e-8;
    public static final ForkJoinPool fjPool = new ForkJoinPool();
    public static boolean DEBUG;
    public static List<Integer> runtime = new LinkedList<>();
    public static boolean NO_RECOMB = false;
    public static int REPEATS;
    public static int DESIRED_REPEATS;
    public static int STEPSIZE = 50;
    public static int PARALLEL_RESTARTS_UPPER_BOUND = 10;
    public static boolean PARALLEL_JHMM = true;
    public static boolean PARALLEL_RESTARTS = false;
    public static double MAX_LLH = Double.NEGATIVE_INFINITY;
    public static boolean LOG_BIC = false;
    public static boolean LOGGING = false;
    public static boolean PRINT = true;
    public static boolean MODELSELECTION;
    public static int SAMPLING_NUMBER;
    public static StringBuilder LOG = new StringBuilder();
    private static final long start = System.currentTimeMillis();

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

    public static void incPercentage() {
        PERCENTAGE += 100d / Globals.REPEATS;
    }

    public static void printBIC(int K, int bic) {
        if (MODELSELECTION) {
            System.out.print("\r" + time() + " Model selection [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[BIC: " + (int) bic + "]                 ");
        } else {
            System.out.print("\r" + time() + " Model training  [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[BIC: " + (int) bic + "]                 ");
        }
    }

    public static void print(String s) {
        System.out.print("\r" + time() + " " + s);
    }
    public static void printPercentage(int K) {
        if (!DEBUG)
        if (MODELSELECTION) {
            System.out.print("\r" + time() + " Model selection [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[LLH: " +  MAX_LLH + "]                 ");
        } else {
            System.out.print("\r" + time() + " Model training  [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[LLH: " +  MAX_LLH + "]                 ");
        }
    }

    private static String time() {
        return new SimpleDateFormat("hh:mm:ss:SSS").format(new Date(System.currentTimeMillis() - start));
    }

    public static synchronized double getMAX_LLH() {
        return MAX_LLH;
    }

    public static synchronized void maxMAX_LLH(double llh) {
        MAX_LLH = Math.max(MAX_LLH, llh);
    }
}
