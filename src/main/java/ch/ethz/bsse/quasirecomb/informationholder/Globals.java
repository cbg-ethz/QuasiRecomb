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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;

/**
 * Information holder for all necessary given and inferred parameters.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Globals {

    private static final Globals INSTANCE = new Globals();

    public Globals getInstance() {
        return INSTANCE;
    }

    private Globals() {
        df.setTimeZone(TimeZone.getTimeZone("GMT"));
    }
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
    private double PCHANGE;
    private double BETA_Z;
    private double ALPHA_Z;
    private double ALPHA_H;
    private double ESTIMATION_EPSILON;
    private double SAMPLING_EPSILON;
    private double DELTA_LLH;
    private double DELTA_REFINE_LLH;
    private double MAX_LLH = Double.NEGATIVE_INFINITY;
    private double[][] MU_PRIOR;
    private int ALIGNMENT_BEGIN = Integer.MAX_VALUE;
    private int ALIGNMENT_END = Integer.MIN_VALUE;
    private int REPEATS;
    private int DESIRED_REPEATS;
    private int STEPSIZE = 100;
    private int PARALLEL_RESTARTS_UPPER_BOUND = 10;
    private int SAMPLING_NUMBER;
    private int PERTURB;
    private final int cpus = Runtime.getRuntime().availableProcessors();
    private List<Integer> runtime = new LinkedList<>();
    private long start = System.currentTimeMillis();
    private String OPTIMUM;
    private String[] HAPLOTYPE_ARRAY_EMPIRICAL;
    private String SAVEPATH;
    private StringBuilder LOG = new StringBuilder();
    private final DateFormat df = new SimpleDateFormat("HH:mm:ss");
//    private final ExecutorService executor = Executors.newFixedThreadPool(1);
    private final ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1);
    private final ForkJoinPool fjPool = new ForkJoinPool();

    public synchronized void log(Object o) {
        if (PRINT) {
            System.out.print(o);
        } else {
            if (LOGGING) {
                LOG.append(o);
            }
        }
    }
    private double PERCENTAGE = 0;

    public void incPercentage() {
        PERCENTAGE += 100d / REPEATS;
    }

    public void printBIC(int K, int bic) {
        if (MODELSELECTION) {
            System.out.print("\r" + time() + " Model selection [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[BIC: " + (int) bic + "]                 ");
        } else {
            System.out.print("\r" + time() + " Model training  [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[BIC: " + (int) bic + "]                 ");
        }
    }
    private String oldOut = "";

    public void print(String s) {
        if (!oldOut.equals(s)) {
            this.oldOut = s;
            System.out.print("\r" + time() + " " + s);
        }
    }

    public void println(String s) {
        System.out.print("\n" + time() + " " + s);
    }

    public void printPercentage(int K) {
        if (!DEBUG) {
            if (MODELSELECTION) {
                System.out.print("\r" + time() + " Model selection [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[LLH: " + MAX_LLH + "]                 ");
            } else {
                System.out.print("\r" + time() + " Model training  [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[LLH: " + MAX_LLH + "]                 ");
            }
        }
    }
    private long oldTime = 0;

    public void resetTimer() {
        this.oldTime = 0;
    }

    public void printPercentage(int K, double read) {
        if (!oldOut.equals(time())) {
            this.oldOut = time();
            if (oldTime == 0) {
                oldTime = System.currentTimeMillis();
            }
            long time = System.currentTimeMillis() - oldTime;
            System.out.print("\r" + time() + " Model selection [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[ ETA:" + df.format(new Date((long) ((1 - read) * time / read))) + "]");
        }
    }
    private double hammingCount = 0;
    private int hammingMax = 0;

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

    public synchronized void maxMAX_LLH(double llh) {
        MAX_LLH = Math.max(MAX_LLH, llh);
    }

    public static Globals getINSTANCE() {
        return INSTANCE;
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

    public int getSTEPSIZE() {
        return STEPSIZE;
    }

    public int getPARALLEL_RESTARTS_UPPER_BOUND() {
        return PARALLEL_RESTARTS_UPPER_BOUND;
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

    public void setSTEPSIZE(int STEPSIZE) {
        this.STEPSIZE = STEPSIZE;
    }

    public void setPARALLEL_RESTARTS_UPPER_BOUND(int PARALLEL_RESTARTS_UPPER_BOUND) {
        this.PARALLEL_RESTARTS_UPPER_BOUND = PARALLEL_RESTARTS_UPPER_BOUND;
    }

    public void setPARALLEL_RESTARTS(boolean PARALLEL_RESTARTS) {
        this.PARALLEL_RESTARTS = PARALLEL_RESTARTS;
    }

    public void setMAX_LLH(double MAX_LLH) {
        this.MAX_LLH = MAX_LLH;
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

    public double[][] getMU_PRIOR() {
        return MU_PRIOR;
    }

    public void setMU_PRIOR(double[][] MU_PRIOR) {
        this.MU_PRIOR = MU_PRIOR;
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
    
}
