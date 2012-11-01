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
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ReadHMMWorker extends RecursiveTask<Pair<List<ReadHMM>, EInfo>> {

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
    protected Pair<List<ReadHMM>, EInfo> compute() {
        if (end - start < Globals.getINSTANCE().getSTEPSIZE()) {
            List<ReadHMM> list = new LinkedList<>();
            EInfo einfo = new EInfo(jhmm.getK(), jhmm.getL(), jhmm.getn());
            for (int i = start; i < end; i++) {
                ReadHMM r = new ReadHMM(jhmm, reads[i]);
                list.add(r);
                int offset = r.getBegin();
                int times = r.getCount();
                //CONTINUE
                for (int j = 0; j < r.getLength(); j++) {
                    if (r.getRead().isHit(j)) {
                        int jGlobal = offset + j;
                        for (int k = 0; k < jhmm.getK(); k++) {
                            einfo.nJK[jGlobal][k] += r.gamma(j, k) * times;
                            if (j > 0) {
                                for (int l = 0; l < jhmm.getK(); l++) {
                                    einfo.nJKL[jGlobal][k][l] += r.xi(j, k, l) * times;
                                    if (k == l) {
                                        einfo.nJeq[jGlobal] += r.xi(j, k, l) * times;
                                    } else {
                                        einfo.nJneq[jGlobal] += r.xi(j, k, l) * times;
                                    }
                                }
                            }
                            for (int v = 0; v < jhmm.getn(); v++) {
                                einfo.nJKV[jGlobal][k][v] += r.gamma(j, k, v) * times;
                            }
                        }

                        for (int v = 0; v < jhmm.getn(); v++) {
                            byte b = r.getSequence(j);
                            for (int k = 0; k < jhmm.getK(); k++) {
                                if (v != b) {
                                    einfo.nneqPos[jGlobal] += r.gamma(j, k, v) * times;
                                } else {
                                    einfo.neqPos[jGlobal] += r.gamma(j, k, v) * times;
                                }
                            }
                        }
                    }
                }
            }
            return Pair.with(list, einfo);
        } else {
            final int mid = start + (end - start) / 2;
            ReadHMMWorker left = new ReadHMMWorker(jhmm, reads, start, mid);
            ReadHMMWorker right = new ReadHMMWorker(jhmm, reads, mid, end);
            left.fork();
            List<ReadHMM> list = new LinkedList<>();
            Pair<List<ReadHMM>, EInfo> compute = right.compute();
            list.addAll(compute.getValue0());
            Pair<List<ReadHMM>, EInfo> join = left.join();
            list.addAll(join.getValue0());
            EInfo e = compute.getValue1();
            e.add(join.getValue1());
            return Pair.with(list, e);
        }
    }
}
