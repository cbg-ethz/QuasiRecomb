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
package ch.ethz.bsse.quasirecomb.model.hmm.parallel;

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.JHMMInterface;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMMStatic;
import java.util.concurrent.Callable;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class CallableReadHMMList implements Callable<Double> {

    private JHMMInterface jhmm;
    private Read[] reads;

    public CallableReadHMMList(JHMMInterface jhmm, Read[] reads) {
        this.jhmm = jhmm;
        this.reads = reads;
    }
    @Override 
    public Double call() throws Exception {
        double d = 0;
        for (int i = 0; i < reads.length; i++) {
            d += ReadHMMStatic.computeFB(jhmm, reads[i]);
        }
        return d;
    }
    
}
