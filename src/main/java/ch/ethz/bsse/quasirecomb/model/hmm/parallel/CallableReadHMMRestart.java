/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.model.hmm.parallel;

import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMM;
import java.util.concurrent.Callable;

/**
 *
 * @author XLR
 */
public class CallableReadHMMRestart implements Callable<Void> {

    private ReadHMM read;

    public CallableReadHMMRestart(ReadHMM read) {
        this.read = read;
    }
    @Override
    public Void call() throws Exception {
        read.recalc();
        return null;
    }
    
}
