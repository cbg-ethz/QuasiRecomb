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
package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;
import org.javatuples.Triplet;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
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
