/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.model.hmm.parallel;

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.JHMM;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMM;
import java.util.concurrent.Callable;

/**
 *
 * @author XLR
 */
public class CallableReadHMM implements Callable<ReadHMM> {

    private JHMM jhmm;
    private Read read;

    public CallableReadHMM(JHMM jhmm, Read read) {
        this.jhmm = jhmm;
        this.read = read;
    }
    @Override
    public ReadHMM call() throws Exception {
        return new ReadHMM(jhmm, read);
    }
    
}
