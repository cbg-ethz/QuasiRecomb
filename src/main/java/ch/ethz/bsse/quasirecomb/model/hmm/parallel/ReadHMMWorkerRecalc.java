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
public class ReadHMMWorkerRecalc extends RecursiveTask<EInfo> {

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
    protected EInfo compute() {
        if (end - start < Globals.getINSTANCE().getSTEPSIZE()) {
            EInfo einfo = new EInfo(jhmm.getK(), jhmm.getL(), jhmm.getn());
            for (int i = start; i < end; i++) {
                ReadHMM r = this.readHMMArray[i];
                r.recalc();
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
            return einfo;
        } else {
            final int mid = start + (end - start) / 2;
            ReadHMMWorkerRecalc left = new ReadHMMWorkerRecalc(jhmm, readHMMArray, start, mid);
            ReadHMMWorkerRecalc right = new ReadHMMWorkerRecalc(jhmm, readHMMArray, mid, end);
            left.fork();
            EInfo einfo = right.compute();
            einfo.add(left.join());
            return einfo;
        }
    }
}
