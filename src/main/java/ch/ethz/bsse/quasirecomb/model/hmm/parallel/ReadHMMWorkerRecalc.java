/**
 * Copyright (c) 2011-2012 Armin Töpfer
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
import ch.ethz.bsse.quasirecomb.model.hmm.JHMM;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMM;
import java.util.concurrent.RecursiveTask;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadHMMWorkerRecalc extends RecursiveTask<Void> {

    private static final long serialVersionUID = 1L;
    private int start;
    private int end;
    private JHMM jhmm;
    private ReadHMM[] readHMMArray;

    public ReadHMMWorkerRecalc(JHMM jhmm, ReadHMM[] readHMMArray, int start, int end) {
        this.readHMMArray = readHMMArray;
        this.start = start;
        this.end = end;
        this.jhmm = jhmm;
    }

    @Override
    protected Void compute() {
        if (end - start < Globals.getINSTANCE().getSTEPSIZE()) {
            for (int i = start; i < end; i++) {
                ReadHMM r = this.readHMMArray[i];
                r.recalc();
            }
            return null;
        } else {
            final int mid = start + (end - start) / 2;
            ReadHMMWorkerRecalc left = new ReadHMMWorkerRecalc(jhmm, readHMMArray, start, mid);
            ReadHMMWorkerRecalc right = new ReadHMMWorkerRecalc(jhmm, readHMMArray, mid, end);
            left.fork();
            right.compute();
            return null;
        }
    }
}
