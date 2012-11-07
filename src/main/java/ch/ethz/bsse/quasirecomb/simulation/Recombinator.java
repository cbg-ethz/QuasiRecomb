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
package ch.ethz.bsse.quasirecomb.simulation;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import com.google.common.collect.Sets;
import java.util.*;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Recombinator {

    public static void recombine(String path, int[] spots, String output) {
        Map<String, String> haplotypes = Utils.parseHaplotypeFile(path);
        if (Globals.getINSTANCE().isDEBUG()) {
            System.out.println(Arrays.toString(spots));
            for (String s : haplotypes.keySet()) {
                System.out.println(s);
            }
        }

        List<Integer> spotsList = new LinkedList<>();
        spotsList.add(0);

        for (int i : spots) {
            spotsList.add(i);
        }

        spotsList.add(haplotypes.keySet().iterator().next().length());

        List<Set<String>> sets = new ArrayList<>();
        for (int i = 1; i < spotsList.size(); i++) {
            Set<String> tmp_set = new HashSet<>();
            for (String s : haplotypes.keySet()) {
                tmp_set.add(s.substring(spotsList.get(i - 1), spotsList.get(i)));
            }
            sets.add(tmp_set);
        }

        Set<List<String>> cartesianProduct = Sets.cartesianProduct(sets);
        List<String> recombinants = new ArrayList<>();
        for (List<String> ll : cartesianProduct) {
            StringBuilder sb = new StringBuilder();
            for (String s : ll) {
                sb.append(s);
            }
            recombinants.add(sb.toString());
        }
        StringBuilder sb = new StringBuilder();
        int i = 0;
        for (String s : haplotypes.keySet()) {
            recombinants.remove(s);
            sb.append(">generator-").append(i++).append("\n").append(s).append("\n");
        }
        i = 0;
        for (String s : recombinants) {
//            System.out.println(s);
            sb.append(">recombinant-").append(i++).append("\n").append(s).append("\n");
        }
//        if (output.endsWith(File.separator)) {
//            output+="recombinants.fasta";
//        }
        System.out.println(sb.toString());
//        Utils.saveFile(output, sb.toString());
    }
}
