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
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.model.hmm.JHMM;
import ch.ethz.bsse.quasirecomb.model.hmm.ReadHMM;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.RecursiveTask;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadHMMWorker extends RecursiveTask<List<ReadHMM>> {

    private static final long serialVersionUID = 1L;
    private JHMM jhmm;
    private Read[] reads;
    private int start;
    private int end;

    public ReadHMMWorker(JHMM jhmm, Read[] reads, int start, int end) {
        this.jhmm = jhmm;
        this.reads = reads;
        this.start = start;
        this.end = end;
    }

    @Override
    protected List<ReadHMM> compute() {
        if (end - start <= Globals.getINSTANCE().getSTEPSIZE()) {
            List<ReadHMM> list = new LinkedList<>();
            for (int i = start; i < end; i++) {
                ReadHMM r = new ReadHMM(jhmm, reads[i]);
                list.add(r);
            }
            return list;
        } else {
            final int mid = start + (end - start) / 2;
            ReadHMMWorker left = new ReadHMMWorker(jhmm, reads, start, mid);
            ReadHMMWorker right = new ReadHMMWorker(jhmm, reads, mid, end);
            left.fork();
            List<ReadHMM> list = right.compute();
            list.addAll(left.join());
            return list;
        }
    }
}
