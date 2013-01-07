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
package ch.ethz.bsse.quasirecomb.informationholder;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;
import java.util.TimeZone;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Information holder for all necessary given and inferred parameters.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Globals {

    private static final Globals INSTANCE = new Globals();
    private static final BlockingQueue<Runnable> blockingQueue = new ArrayBlockingQueue<>(Runtime.getRuntime().availableProcessors() - 1);
    private static final RejectedExecutionHandler rejectedExecutionHandler = new ThreadPoolExecutor.CallerRunsPolicy();
    private static ExecutorService executor = refreshExecutor();

    private static ExecutorService refreshExecutor() {
//        return Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1);
        return new ThreadPoolExecutor(Runtime.getRuntime().availableProcessors() - 1, Runtime.getRuntime().availableProcessors() - 1, 0L, TimeUnit.MILLISECONDS, blockingQueue, rejectedExecutionHandler);
    }

    public static void renewExecutor() {
        executor = refreshExecutor();
    }

    public static Globals getINSTANCE() {
        return INSTANCE;
    }
    private boolean BIAS_MU;
    private boolean SILENT;
    private boolean ML;
    private boolean STOP_QUICK;
    private boolean PRINT_ALIGNMENT;
    private boolean CIRCOS;
    private boolean NOSAMPLE;
    private boolean OVERLAP;
    private boolean UNINFORMATIVE_EPSILON_PRIOR;
    private boolean PDELTA;
    private boolean PLOT;
    private boolean STORAGE;
    private boolean SNAPSHOTS;
    private boolean FLAT_EPSILON_PRIOR;
    private boolean DEBUG;
    private boolean NO_RECOMB = false;
    private boolean FORCE_NO_RECOMB = false;
    private boolean PARALLEL_RESTARTS = false;
    private boolean LOG_BIC = false;
    private boolean LOGGING = false;
    private boolean PRINT = true;
    private boolean MODELSELECTION;
    private boolean PAIRED = false;
    private boolean PRIORMU;
    private boolean SPIKERHO;
    private double PCHANGE;
    private double MULT_RHO;
    private double MULT_MU;
    private double BETA_Z;
    private double ALPHA_Z;
    private double ALPHA_H;
    private double ESTIMATION_EPSILON;
    private double SAMPLING_EPSILON;
    private double DELTA_LLH;
    private double DELTA_REFINE_LLH;
    private boolean PRUNE;
    private double INTERPOLATE_MU;
    private double INTERPOLATE_RHO;
    private double CURRENT_DELTA_LLH = 0;
    private double MAX_LLH = -1;
    private double MIN_BIC = Double.MIN_VALUE;
    private int ALIGNMENT_BEGIN = Integer.MAX_VALUE;
    private int ALIGNMENT_END = Integer.MIN_VALUE;
    private int STEPS;
    private int REPEATS;
    private int DESIRED_REPEATS;
    private int SAMPLING_NUMBER;
    private int PERTURB;
    private final int cpus = Runtime.getRuntime().availableProcessors();
    private List<Integer> runtime = new LinkedList<>();
    private long start = System.currentTimeMillis();
    private String GENOME;
    private String OPTIMUM;
    private String[] HAPLOTYPE_ARRAY_EMPIRICAL;
    private int[][] MU_PRIOR;
    private String SAVEPATH;
    private StringBuilder LOG = new StringBuilder();
    private final DateFormat df = new SimpleDateFormat("HH:mm:ss:SSS");
//    private final ExecutorService executor = new ThreadPoolExecutor(Runtime.getRuntime().availableProcessors() - 1, Runtime.getRuntime().availableProcessors() - 1, 0L, TimeUnit.MILLISECONDS, blockingQueue, rejectedExecutionHandler);
    private final ForkJoinPool fjPool = new ForkJoinPool();
    private final AtomicInteger MERGED = new AtomicInteger(0);
    private double PERCENTAGE = 0;
    private long oldTime = 0;
    private String oldOut = "";
    private double hammingCount = 0;
    private int hammingMax = 0;
    
    private Globals() {
        df.setTimeZone(TimeZone.getTimeZone("GMT"));
    }
    public Globals getInstance() {
        return INSTANCE;
    }

    public synchronized void log(Object o) {
        if (PRINT) {
            System.out.print(o);
        } else {
            if (LOGGING) {
                LOG.append(o);
            }
        }
    }

    public void incPercentage() {
        PERCENTAGE += 100d / REPEATS;
    }

    public void printBIC(int K, int bic) {
        System.out.print("\r                                                                                                                                                   ");
        if (MODELSELECTION) {
            System.out.print("\r" + time() + " Model selection [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[BIC: " + (int) bic + "]                 ");
        } else {
            System.out.print("\r" + time() + " Model training  [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[BIC: " + (int) bic + "]                 ");
        }
    }

    public void print(String s) {
        if (!oldOut.equals(s)) {
            this.oldOut = s;
            System.out.print("\r" + time() + " " + s);
        }
    }

    public void println(String s) {
        System.out.print("\n" + time() + " " + s);
    }

    public void resetTimer() {
        this.oldTime = 0;
    }

    public void printPercentage(int K, double read, double Kmin) {
        if (!SILENT) {
            if (!oldOut.equals(time())) {
                this.oldOut = time();
                if (oldTime == 0) {
                    oldTime = System.currentTimeMillis();
                }
                long time = System.currentTimeMillis() - oldTime;
                System.out.print("\r                                                                                                                                                   ");
//            System.out.print("\r" + time() + " Model " + (MODELSELECTION ? "selection" : "training") + " [K " + (int)Kmin + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "% [ETA:" + df.format(new Date((long) ((1 - read) * time / read))) + "]" + "[cK " + K + "]" + "[LLH " + ((int)MAX_LLH *1000)/100d + "]" + "[BIC " + ((int)MIN_BIC *1000)/100d + "]" + "[D-LLH " + Summary.shorten(CURRENT_DELTA_LLH) + "]");
                System.out.print("\r" + time() + " Model " + (MODELSELECTION ? "selection" : "training") + " [K " + (int) Kmin + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "% [ETA:" + df.format(new Date((long) ((1 - read) * time / read))) + "]" + "[cK " + K + "]");
            }
        }
    }

    public synchronized void incHamming(int inc) {
        hammingCount += inc * (100d / hammingMax);
    }

    public synchronized void printHamming(int inc) {
        incHamming(inc);
        System.out.print("\r" + time() + " Computing:\t" + Math.round(hammingCount * 1000) / 1000 + "%");
    }

    public String time() {
        return df.format(new Date(System.currentTimeMillis() - start));
    }

    public void setSPIKERHO(boolean SPIKERHO) {
        this.SPIKERHO = SPIKERHO;
    }

    public boolean isSPIKERHO() {
        return SPIKERHO;
    }

    public boolean isPRIORMU() {
        return PRIORMU;
    }

    public void setPRIORMU(boolean PRIORMU) {
        this.PRIORMU = PRIORMU;
    }

    public int[][] getMU_PRIOR() {
        return MU_PRIOR;
    }

    public void setMU_PRIOR(int[][] MU_PRIOR) {
        this.MU_PRIOR = MU_PRIOR;
    }

    public double getCURRENT_DELTA_LLH() {
        return CURRENT_DELTA_LLH;
    }

    public void setCURRENT_DELTA_LLH(double CURRENT_DELTA_LLH) {
        this.CURRENT_DELTA_LLH = CURRENT_DELTA_LLH;
    }

    public double getMIN_BIC() {
        return MIN_BIC;
    }

    public void setPRUNE(boolean PRUNE) {
        this.PRUNE = PRUNE;
    }

    public boolean isPRUNE() {
        return PRUNE;
    }

    public int getMERGED() {
        return MERGED.get();
    }

    public void incMERGED() {
        MERGED.getAndIncrement();
    }

    public void setOVERLAP(boolean OVERLAP) {
        this.OVERLAP = OVERLAP;
    }

    public boolean isOVERLAP() {
        return OVERLAP;
    }

    public void setUNINFORMATIVE_EPSILON_PRIOR(boolean UNINFORMATIVE_EPSILON_PRIOR) {
        this.UNINFORMATIVE_EPSILON_PRIOR = UNINFORMATIVE_EPSILON_PRIOR;
    }

    public boolean isUNINFORMATIVE_EPSILON_PRIOR() {
        return UNINFORMATIVE_EPSILON_PRIOR;
    }

    public boolean isPDELTA() {
        return PDELTA;
    }

    public void setPDELTA(boolean PDELTA) {
        this.PDELTA = PDELTA;
    }

    public void setDELTA_REFINE_LLH(double DELTA_REFINE_LLH) {
        this.DELTA_REFINE_LLH = DELTA_REFINE_LLH;
    }

    public double getDELTA_REFINE_LLH() {
        return DELTA_REFINE_LLH;
    }

    public boolean isPLOT() {
        return PLOT;
    }

    public void setPLOT(boolean PLOT) {
        this.PLOT = PLOT;
    }

    public String getSnapshotDir() {
        return Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "snapshots" + File.separator;
    }

    public boolean isSNAPSHOTS() {
        return SNAPSHOTS;
    }

    public void setSNAPSHOTS(boolean SNAPSHOTS) {
        this.SNAPSHOTS = SNAPSHOTS;
    }

    public ExecutorService getExecutor() {
        return executor;
    }

    public void setHammingMax(int hammingMax) {
        this.hammingMax = hammingMax;
    }

    public void setFORCE_NO_RECOMB(boolean FORCE_NO_RECOMB) {
        this.FORCE_NO_RECOMB = FORCE_NO_RECOMB;
    }

    public boolean isFORCE_NO_RECOMB() {
        return FORCE_NO_RECOMB;
    }

    public synchronized double getMAX_LLH() {
        return MAX_LLH;
    }

    public synchronized void minMIN_BIC(double bic) {
        MIN_BIC = Math.min(MIN_BIC, bic);
    }

    public synchronized void maxMAX_LLH(double llh) {
        if (MAX_LLH == -1 || llh > MAX_LLH) {
            MAX_LLH = llh;
        }
    }

    public void setSTOP_QUICK(boolean STOP_QUICK) {
        this.STOP_QUICK = STOP_QUICK;
    }

    public boolean isSTOP_QUICK() {
        return STOP_QUICK;
    }

    public int getALIGNMENT_BEGIN() {
        return ALIGNMENT_BEGIN;
    }

    public int getALIGNMENT_END() {
        return ALIGNMENT_END;
    }

    public boolean isFLAT_EPSILON_PRIOR() {
        return FLAT_EPSILON_PRIOR;
    }

    public double getPCHANGE() {
        return PCHANGE;
    }

    public double getBETA_Z() {
        return BETA_Z;
    }

    public double getALPHA_Z() {
        return ALPHA_Z;
    }

    public double getALPHA_H() {
        return ALPHA_H;
    }

    public String[] getHAPLOTYPE_ARRAY_EMPIRICAL() {
        return HAPLOTYPE_ARRAY_EMPIRICAL;
    }

    public String getSAVEPATH() {
        return SAVEPATH;
    }

    public double getESTIMATION_EPSILON() {
        return ESTIMATION_EPSILON;
    }

    public double getSAMPLING_EPSILON() {
        return SAMPLING_EPSILON;
    }

    public double getDELTA_LLH() {
        return DELTA_LLH;
    }

    public ForkJoinPool getFjPool() {
        return fjPool;
    }

    public boolean isDEBUG() {
        return DEBUG;
    }

    public List<Integer> getRuntime() {
        return runtime;
    }

    public boolean isNO_RECOMB() {
        return NO_RECOMB;
    }

    public int getREPEATS() {
        return REPEATS;
    }

    public int getDESIRED_REPEATS() {
        return DESIRED_REPEATS;
    }

    public boolean isPARALLEL_RESTARTS() {
        return PARALLEL_RESTARTS;
    }

    public boolean isLOG_BIC() {
        return LOG_BIC;
    }

    public boolean isLOGGING() {
        return LOGGING;
    }

    public boolean isPRINT() {
        return PRINT;
    }

    public boolean isMODELSELECTION() {
        return MODELSELECTION;
    }

    public int getSAMPLING_NUMBER() {
        return SAMPLING_NUMBER;
    }

    public StringBuilder getLOG() {
        return LOG;
    }

    public long getStart() {
        return start;
    }

    public double getPERCENTAGE() {
        return PERCENTAGE;
    }

    public void setALIGNMENT_BEGIN(int ALIGNMENT_BEGIN) {
        this.ALIGNMENT_BEGIN = ALIGNMENT_BEGIN;
    }

    public void setALIGNMENT_END(int ALIGNMENT_END) {
        this.ALIGNMENT_END = ALIGNMENT_END;
    }

    public void setFLAT_EPSILON_PRIOR(boolean FLAT_EPSILON_PRIOR) {
        this.FLAT_EPSILON_PRIOR = FLAT_EPSILON_PRIOR;
    }

    public void setPCHANGE(double PCHANGE) {
        this.PCHANGE = PCHANGE;
    }

    public void setBETA_Z(double BETA_Z) {
        this.BETA_Z = BETA_Z;
    }

    public void setALPHA_Z(double ALPHA_Z) {
        this.ALPHA_Z = ALPHA_Z;
    }

    public void setALPHA_H(double ALPHA_H) {
        this.ALPHA_H = ALPHA_H;
    }

    public void setHAPLOTYPE_ARRAY_EMPIRICAL(String[] HAPLOTYPE_ARRAY_EMPIRICAL) {
        this.HAPLOTYPE_ARRAY_EMPIRICAL = HAPLOTYPE_ARRAY_EMPIRICAL;
    }

    public void setSAVEPATH(String SAVEPATH) {
        this.SAVEPATH = SAVEPATH;
    }

    public void setESTIMATION_EPSILON(double ESTIMATION_EPSILON) {
        this.ESTIMATION_EPSILON = ESTIMATION_EPSILON;
    }

    public void setSAMPLING_EPSILON(double SAMPLING_EPSILON) {
        this.SAMPLING_EPSILON = SAMPLING_EPSILON;
    }

    public void setDELTA_LLH(double DELTA_LLH) {
        this.DELTA_LLH = DELTA_LLH;
    }

    public void setDEBUG(boolean DEBUG) {
        this.DEBUG = DEBUG;
    }

    public void setRuntime(List<Integer> runtime) {
        this.runtime = runtime;
    }

    public void setNO_RECOMB(boolean NO_RECOMB) {
        this.NO_RECOMB = NO_RECOMB;
    }

    public void setREPEATS(int REPEATS) {
        this.REPEATS = REPEATS;
    }

    public void setDESIRED_REPEATS(int DESIRED_REPEATS) {
        this.DESIRED_REPEATS = DESIRED_REPEATS;
    }

    public void setPARALLEL_RESTARTS(boolean PARALLEL_RESTARTS) {
        this.PARALLEL_RESTARTS = PARALLEL_RESTARTS;
    }

    public void setMAX_LLH(double MAX_LLH) {
        this.MAX_LLH = MAX_LLH;
    }

    public void setMIN_BIC(double MAX_BIC) {
        this.MIN_BIC = MAX_BIC;
    }

    public void setLOG_BIC(boolean LOG_BIC) {
        this.LOG_BIC = LOG_BIC;
    }

    public void setLOGGING(boolean LOGGING) {
        this.LOGGING = LOGGING;
    }

    public void setPRINT(boolean PRINT) {
        this.PRINT = PRINT;
    }

    public void setMODELSELECTION(boolean MODELSELECTION) {
        this.MODELSELECTION = MODELSELECTION;
    }

    public void setSAMPLING_NUMBER(int SAMPLING_NUMBER) {
        this.SAMPLING_NUMBER = SAMPLING_NUMBER;
    }

    public void setLOG(StringBuilder LOG) {
        this.LOG = LOG;
    }

    public void setPERCENTAGE(double PERCENTAGE) {
        this.PERCENTAGE = PERCENTAGE;
    }

    public boolean isPAIRED() {
        return PAIRED;
    }

    public void setPAIRED(boolean PAIRED) {
        this.PAIRED = PAIRED;
    }

    public int getCpus() {
        return cpus;
    }

    public boolean isSTORAGE() {
        return STORAGE;
    }

    public void setSTORAGE(boolean STORAGE) {
        this.STORAGE = STORAGE;
    }

    public String getOPTIMUM() {
        return OPTIMUM;
    }

    public void setOPTIMUM(String OPTIMUM) {
        this.OPTIMUM = OPTIMUM;
    }

    public int getPERTURB() {
        return PERTURB;
    }

    public void setPERTURB(int PERTURB) {
        this.PERTURB = PERTURB;
    }

    public boolean isNOSAMPLE() {
        return NOSAMPLE;
    }

    public void setNOSAMPLE(boolean NOSAMPLE) {
        this.NOSAMPLE = NOSAMPLE;
    }

    public boolean isCIRCOS() {
        return CIRCOS;
    }

    public void setCIRCOS(boolean CIRCOS) {
        this.CIRCOS = CIRCOS;
    }

    public double getMULT_RHO() {
        return MULT_RHO;
    }

    public void setMULT_RHO(double MULT_RHO) {
        this.MULT_RHO = MULT_RHO;
    }

    public double getMULT_MU() {
        return MULT_MU;
    }

    public void setMULT_MU(double MULT_MU) {
        this.MULT_MU = MULT_MU;
    }

    public String getGENOME() {
        return GENOME;
    }

    public void setGENOME(String GENOME) {
        this.GENOME = GENOME;
    }

    public boolean isPRINT_ALIGNMENT() {
        return PRINT_ALIGNMENT;
    }

    public void setPRINT_ALIGNMENT(boolean PRINT_ALIGNMENT) {
        this.PRINT_ALIGNMENT = PRINT_ALIGNMENT;
    }

    public boolean isML() {
        return ML;
    }

    public void setML(boolean ML) {
        this.ML = ML;
    }

    public boolean isSILENT() {
        return SILENT;
    }

    public void setSILENT(boolean SILENT) {
        this.SILENT = SILENT;
    }

    public boolean isBIAS_MU() {
        return BIAS_MU;
    }

    public void setBIAS_MU(boolean BIAS_MU) {
        this.BIAS_MU = BIAS_MU;
    }

    public int getSTEPS() {
        return STEPS;
    }

    public void setSTEPS(int STEPS) {
        this.STEPS = STEPS;
    }

    public double getINTERPOLATE_MU() {
        return INTERPOLATE_MU;
    }

    public void setINTERPOLATE_MU(double INTERPOLATE_MU) {
        this.INTERPOLATE_MU = INTERPOLATE_MU;
    }

    public double getINTERPOLATE_RHO() {
        return INTERPOLATE_RHO;
    }

    public void setINTERPOLATE_RHO(double INTERPOLATE_RHO) {
        this.INTERPOLATE_RHO = INTERPOLATE_RHO;
    }
}
