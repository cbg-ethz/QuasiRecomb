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
import ch.ethz.bsse.quasirecomb.utils.DistanceUtils;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.RecursiveTask;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class HammerInnerWorker extends RecursiveTask<Map<String, Integer>> {

    private static final long serialVersionUID = 1L;
    private String[] fasta;
    private int start;
    private int end;
    private int index;

    public HammerInnerWorker(String[] fasta, int index, int start, int end) {
        this.start = start;
        this.end = end;
        this.fasta = fasta;
        this.index = index;
    }

    @Override
    protected Map<String, Integer> compute() {
        if (end - start <= 1) {
            Map<String, Integer> mapInner = new HashMap<>();
            for (int i = start; i < end; i++) {
                mapInner.put(fasta[i], DistanceUtils.calcHamming(fasta[index], fasta[i]));
                Globals.getINSTANCE().printHamming(1);
            }
            return mapInner;
        } else {
            final int mid = start + (end - start) / 2;
            HammerInnerWorker left = new HammerInnerWorker(fasta, index, start, mid);
            HammerInnerWorker right = new HammerInnerWorker(fasta, index, mid, end);
            left.fork();
            Map<String, Integer> compute = right.compute();
            compute.putAll(left.join());
            return compute;
        }
    }
}
