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

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.model.hmm.SingleEM;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.RecursiveTask;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class RestartWorker extends RecursiveTask<List<OptimalResult>> {

    private int N;
    private int K;
    private int L;
    private int n;
    private Map<byte[], Integer> reads;
    private byte[][] haplotypesArray;
    private OptimalResult or;
    private double delta;
    private int start;
    private int end;

    public RestartWorker(int N, int K, int L, int n, Map<byte[], Integer> reads, byte[][] haplotypesArray, double delta, int start, int end) {
        this.N = N;
        this.K = K;
        this.L = L;
        this.n = n;
        this.reads = reads;
        this.haplotypesArray = haplotypesArray;
        this.delta = delta;
        this.start = start;
        this.end = end;
    }

    
    @Override
    protected List<OptimalResult> compute() {
        if (end - start <= Globals.PARALLEL_RESTARTS_UPPER_BOUND || !Globals.PARALLEL_RESTARTS) {
            final List<OptimalResult> list = new ArrayList<>();
            Globals.log("+"+start+":"+end);
            for (int i = start; i < end; i++) {
                Globals.log("\tx:"+i);
                final SingleEM singleEm = new SingleEM(N, K, L, n, reads, haplotypesArray, delta);
                list.add(singleEm.getOptimalResult());
            }
            Globals.log("-"+start+":"+end);
            return list;
        } else {
            final int mid = start + (end - start) / 2;
            final RestartWorker left = new RestartWorker(N, K, L, n, reads, haplotypesArray,  delta, start, mid);
            final RestartWorker right = new RestartWorker(N, K, L, n, reads, haplotypesArray,  delta, mid, end);
            left.fork();
            final List<OptimalResult> list = new LinkedList<>();
            list.addAll(right.compute());
            list.addAll(left.join());
            return list;
        }
    }
}
