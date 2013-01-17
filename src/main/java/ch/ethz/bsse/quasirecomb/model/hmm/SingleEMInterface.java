/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;

/**
 *
 * @author toepfera
 */
public interface SingleEMInterface {

    void calcBic();

    double getLoglikelihood();

    OptimalResult getOptimalResult();

    String getOptimumPath();

    void printMeanTime();
    
}
