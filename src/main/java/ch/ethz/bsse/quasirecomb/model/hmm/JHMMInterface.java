package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;
import org.javatuples.Triplet;

/**
 *
 * @author toepfera
 */
public interface JHMMInterface {
    
    void computeSNVPosterior();

    void free(int id);

    double[][] getSnv();

    double[] getAntieps();

    double getBeta();

    double[] getEps();

    int getK();

    int getL();

    double getLoglikelihood();

    double[][][] getMu();

    int getMuChanged();

    int getMuFlats();

    int getN();

    int getNjklFlats();

    int getNjkvFlats();

    double[] getPi();

    int getRestart();

    double[][][] getRho();

    int getRhoChanged();

    int getRhoFlats();

    TempJHMMStorage getStorage();

    double[][] getTauOmega();

    int getn();

    void incBeta();

    Triplet<Integer, Integer, Double> minKL();

    void restart();
}
