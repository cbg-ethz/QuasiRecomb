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
package ch.ethz.bsse.quasirecomb.model;

import ch.ethz.bsse.quasirecomb.distance.DistanceUtils;
import ch.ethz.bsse.quasirecomb.model.hmm.ModelSelection;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.simulation.Sampling;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.util.*;
import org.javatuples.Pair;

/**
 * Forwards parameters depending if input is from an experimental dataset or
 * artifical haplotypes from which has to be sampled.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class ArtificialExperimentalForwarder {

    /**
     * Entry point. Forwards invokes of the specified workflow.
     *
     * @param L length of the reads
     * @param exp is it an experimental dataset
     * @param input path to the fasta file
     * @param Kmin minimal amount of generators
     * @param Kmax minimal amount of generators
     * @param n size of the alphabet
     * @param f distribution of the haplotypes if sampling has to be done
     * @param N amount of reads in case exp is false
     */
    public static void forward(String input, int Kmin, int Kmax, int N) {
            byte[][] haplotypesArray = Utils.parseInput(input);

        if (Globals.BOOTSTRAP) {
            Map<byte[], Double> hapMap = new HashMap<>();
            int hapL = haplotypesArray.length;

            for (int i = 0; i < hapL; i++) {
                hapMap.put(haplotypesArray[i], 1d / hapL);
            }

            Frequency<byte[]> f = new Frequency<>(hapMap);

            byte[][] haplotypesBootstrap = new byte[hapL][];
            for (int i = 0; i < hapL; i++) {
                haplotypesBootstrap[i] = f.roll();
            }
            haplotypesArray = haplotypesBootstrap;
        }

        Map<byte[], Integer> clusterReads = Utils.clusterReads(haplotypesArray);

        int n = countChars(haplotypesArray);
        ModelSelection ms = new ModelSelection(clusterReads, Kmin, Kmax, N, haplotypesArray[0].length, n, haplotypesArray);
    }


    private static int countChars(byte[][] bbb) {
        Map<Byte, Boolean> map = new HashMap<>();
        for (byte[] bb : bbb) {
            for (byte b : bb) {
                map.put(b, Boolean.TRUE);
            }
        }
        return map.keySet().size();
    }
}
