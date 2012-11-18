package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;

/**
 *
 * @author XLR
 */
public interface JHMMI {

    void free(int id);

    double[] getAntieps();

    double[] getEps();

    int getK();

    int getL();

    double getLoglikelihood();

    double[][][] getMu();

    int getMuFlats();

    int getN();

    int getNjklFlats();

    int getNjkvFlats();

    int getParametersChanged();

    double[] getPi();

    int getRestart();

    double[][][] getRho();

    int getRhoFlats();

    TempJHMMStorage getStorage();

    double[][] getTauOmega();

    int getn();

    void restart();
    
}
