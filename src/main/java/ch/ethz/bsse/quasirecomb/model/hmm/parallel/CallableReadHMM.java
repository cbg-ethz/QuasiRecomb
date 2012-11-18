package ch.ethz.bsse.quasirecomb.model.hmm.parallel;

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.JHMMI;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMMStatic;
import java.util.concurrent.Callable;

/**
 *
 * @author XLR
 */
public class CallableReadHMM implements Callable<Double> {

    private JHMMI jhmm;
    private Read read;

    public CallableReadHMM(JHMMI jhmm, Read read) {
        this.jhmm = jhmm;
        this.read = read;
    }
    @Override 
    public Double call() throws Exception {
        return ReadHMMStatic.computeFB(jhmm, read);
    }
    
}
