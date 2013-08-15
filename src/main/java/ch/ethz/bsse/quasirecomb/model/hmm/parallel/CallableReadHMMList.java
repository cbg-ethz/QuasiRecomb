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

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.ParallelJHMMStorage;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMMStatic;
import ch.ethz.bsse.quasirecomb.model.hmm.JHMM;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMMStatic_NR;
import java.util.concurrent.Callable;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class CallableReadHMMList implements Callable<Double> {

    private final JHMM jhmm;
    private final Read[] reads;
    private final ParallelJHMMStorage p;

    public CallableReadHMMList(JHMM jhmm, Read[] reads, ParallelJHMMStorage p) {
        this.jhmm = jhmm;
        this.reads = reads;
        this.p = p;
    }

    @Override
    public Double call() throws Exception {
        double d = 0;
        for (Read read : reads) {
            if (Globals.getINSTANCE().isNO_RECOMB()) {
                d += ReadHMMStatic_NR.computeFB(jhmm, read, p);
            } else {
                d += ReadHMMStatic.computeFB(jhmm, read, p);
            }
        }
        return d;
    }
}
