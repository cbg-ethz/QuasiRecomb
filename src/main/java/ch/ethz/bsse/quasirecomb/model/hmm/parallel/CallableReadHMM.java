package ch.ethz.bsse.quasirecomb.model.hmm.parallel;

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.JHMM;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMMStatic;
import java.util.concurrent.Callable;

/**
 *
 * @author XLR
 */
public class CallableReadHMM implements Callable<Double> {

    private JHMM jhmm;
    private Read read;

    public CallableReadHMM(JHMM jhmm, Read read) {
        this.jhmm = jhmm;
        this.read = read;
    }
    @Override
    public Double call() throws Exception {
        return ReadHMMStatic.computeFB(jhmm, read);
    }
    
}
