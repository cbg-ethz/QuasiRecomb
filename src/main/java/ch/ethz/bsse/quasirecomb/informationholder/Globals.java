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
package ch.ethz.bsse.quasirecomb.informationholder;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Information holder for all necessary given and inferred parameters.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Globals {

    private static final Globals INSTANCE = new Globals();

    public static Globals getINSTANCE() {
        return INSTANCE;
    }
    private boolean GRADIENT;
    private boolean ANNEALING;
    private boolean NO_GAPS;
    private boolean MAX;
    private boolean COVERAGE;
    private boolean BOOTSTRAP;
    private boolean REFINEMENT;
    private boolean SAMPLE_PROTEINS;
    private boolean SAMPLE_READS;
    private boolean NO_QUALITY;
    private boolean WINDOW;
    private boolean UNPAIRED;
    private boolean USER_OPTIMUM;
    private boolean BIAS_MU;
    private boolean SILENT;
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
    private double MAX_DEL;
    private double MAX_OVERALL_DEL;
    private double CUTOFF;
    private double PCHANGE;
    private double MULT_RHO;
    private double MULT_RHO_MIN;
    private double MULT_MU;
    private double MULT_MU_MIN;
    private double BETA_Z;
    private double ALPHA_Z;
    private double ALPHA_H;
    private double ESTIMATION_EPSILON;
    private double DELTA_LLH;
    private double DELTA_REFINE_LLH;
    private boolean PRUNE;
    private double INTERPOLATE_MU;
    private double INTERPOLATE_RHO;
    private double CURRENT_DELTA_LLH = 0;
    private double MAX_LLH = -1;
    private double MIN_BIC = Double.MIN_VALUE;
    private int READ_MINLENGTH;
    private int WINDOW_BEGIN;
    private int WINDOW_END;
    private int ALIGNMENT_BEGIN = Integer.MAX_VALUE;
    private int ALIGNMENT_END = Integer.MIN_VALUE;
    private int STEPS;
    private int DESIRED_REPEATS;
    private int SAMPLING_NUMBER;
    private int NREAL;
    private int REPEATS;
    private final int cpus = Runtime.getRuntime().availableProcessors();
    private List<Integer> runtime = new LinkedList<>();
    private String GENOME;
    private String OPTIMUM;
    private String SAVEPATH;
    private StringBuilder LOG = new StringBuilder();
    private final ForkJoinPool fjPool = new ForkJoinPool();
    private final AtomicInteger MERGED_COUNT = new AtomicInteger(0);
    private final AtomicInteger PAIRED_COUNT = new AtomicInteger(0);
    private int hammingMax = 0;
    private TauOmega TAU_OMEGA;

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

    public TauOmega getTAU_OMEGA() {
        return TAU_OMEGA;
    }

    public void setTAU_OMEGA(Read[] reads, int L) {
        this.TAU_OMEGA = new TauOmega(reads, L);
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

    public int getPAIRED_COUNT() {
        return PAIRED_COUNT.get();
    }

    public void incPAIRED() {
        PAIRED_COUNT.getAndIncrement();
    }

    public int getMERGED() {
        return MERGED_COUNT.get();
    }

    public void incMERGED() {
        MERGED_COUNT.getAndIncrement();
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

    public String getSAVEPATH() {
        return SAVEPATH;
    }

    public double getESTIMATION_EPSILON() {
        return ESTIMATION_EPSILON;
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

    public void setSAVEPATH(String SAVEPATH) {
        this.SAVEPATH = SAVEPATH;
    }

    public void setESTIMATION_EPSILON(double ESTIMATION_EPSILON) {
        this.ESTIMATION_EPSILON = ESTIMATION_EPSILON;
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

    public void setLOG(StringBuilder LOG) {
        this.LOG = LOG;
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

    public void setUSER_OPTIMUM(boolean USER_OPTIMUM) {
        this.USER_OPTIMUM = USER_OPTIMUM;
    }

    public boolean isUSER_OPTIMUM() {
        return USER_OPTIMUM;
    }

    public void setUNPAIRED(boolean UNPAIRED) {
        this.UNPAIRED = UNPAIRED;
    }

    public boolean isUNPAIRED() {
        return UNPAIRED;
    }

    public double getCUTOFF() {
        return CUTOFF;
    }

    public void setCUTOFF(double CUTOFF) {
        this.CUTOFF = CUTOFF;
    }

    public boolean isWINDOW() {
        return WINDOW;
    }

    public void setWINDOW(boolean WINDOW) {
        this.WINDOW = WINDOW;
    }

    public int getWINDOW_BEGIN() {
        return WINDOW_BEGIN;
    }

    public void setWINDOW_BEGIN(int WINDOW_BEGIN) {
        this.WINDOW_BEGIN = WINDOW_BEGIN;
    }

    public int getWINDOW_END() {
        return WINDOW_END;
    }

    public void setWINDOW_END(int WINDOW_END) {
        this.WINDOW_END = WINDOW_END;
    }

    public boolean isNO_QUALITY() {
        return NO_QUALITY;
    }

    public void setNO_QUALITY(boolean NO_QUALITY) {
        this.NO_QUALITY = NO_QUALITY;
    }

    public int getREAD_MINLENGTH() {
        return READ_MINLENGTH;
    }

    public void setREAD_MINLENGTH(int READ_MINLENGTH) {
        this.READ_MINLENGTH = READ_MINLENGTH;
    }

    public int getNREAL() {
        return NREAL;
    }

    public void setNREAL(int NREAL) {
        this.NREAL = NREAL;
    }

    public boolean isSAMPLE_READS() {
        return SAMPLE_READS;
    }

    public void setSAMPLE_READS(boolean SAMPLE_READS) {
        this.SAMPLE_READS = SAMPLE_READS;
    }

    public boolean isSAMPLE_PROTEINS() {
        return SAMPLE_PROTEINS;
    }

    public void setSAMPLE_PROTEINS(boolean SAMPLE_PROTEINS) {
        this.SAMPLE_PROTEINS = SAMPLE_PROTEINS;
    }

    public boolean isREFINEMENT() {
        return REFINEMENT;
    }

    public void setREFINEMENT(boolean REFINEMENT) {
        this.REFINEMENT = REFINEMENT;
    }

    public boolean isBOOTSTRAP() {
        return BOOTSTRAP;
    }

    public void setBOOTSTRAP(boolean BOOTSTRAP) {
        this.BOOTSTRAP = BOOTSTRAP;
    }

    public boolean isCOVERAGE() {
        return COVERAGE;
    }

    public void setCOVERAGE(boolean COVERAGE) {
        this.COVERAGE = COVERAGE;
    }

    public boolean isMAX() {
        return MAX;
    }

    public void setMAX(boolean MAX) {
        this.MAX = MAX;
    }

    public int getHammingMax() {
        return hammingMax;
    }

    public int getREPEATS() {
        return REPEATS;
    }

    public void setREPEATS(int REPEATS) {
        this.REPEATS = REPEATS;
    }

    public boolean isNO_GAPS() {
        return NO_GAPS;
    }

    public void setNO_GAPS(boolean NO_GAPS) {
        this.NO_GAPS = NO_GAPS;
    }

    public boolean isANNEALING() {
        return ANNEALING;
    }

    public void setANNEALING(boolean ANNEALING) {
        this.ANNEALING = ANNEALING;
    }

    public boolean isGRADIENT() {
        return GRADIENT;
    }

    public void setGRADIENT(boolean GRADIENT) {
        this.GRADIENT = GRADIENT;
    }

    public double getMULT_MU_MIN() {
        return MULT_MU_MIN;
    }

    public void setMULT_MU_MIN(double MULT_MU_MIN) {
        this.MULT_MU_MIN = MULT_MU_MIN;
    }

    public double getMAX_DEL() {
        return MAX_DEL;
    }

    public void setMAX_DEL(double MAX_DEL) {
        this.MAX_DEL = MAX_DEL;
    }

    public double getMAX_PERC_DEL() {
        return MAX_OVERALL_DEL;
    }

    public void setMAX_OVERALL_DEL(double MAX_OVERALL_DEL) {
        this.MAX_OVERALL_DEL = MAX_OVERALL_DEL;
    }

    public double getMULT_RHO_MIN() {
        return MULT_RHO_MIN;
    }

    public void setMULT_RHO_MIN(double MULT_RHO_MIN) {
        this.MULT_RHO_MIN = MULT_RHO_MIN;
    }
}
