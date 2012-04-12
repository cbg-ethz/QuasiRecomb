/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.simulation;

import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import com.google.common.collect.Sets;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author XLR
 */
public class Recombinator {

    public static void main(String[] args) {
        recombine("/Users/XLR/Dropbox/QuasiRecomb/src/main/resources/haplotypes/generators_3.fasta", new int[]{186}, "");
    }

    public static void recombine(String path, int[] spots, String output) {
        Map<String, String> haplotypes = Utils.parseHaplotypeFile(path);
        if (Globals.DEBUG) {
            System.out.println(Arrays.toString(spots));
            for (String s : haplotypes.keySet()) {
                System.out.println(s);
            }
        }
//        haplotypes = new HashMap<>();
//        haplotypes.put("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA","a");
//        haplotypes.put("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG","g");
//        haplotypes.put("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC","c");
//        haplotypes.put("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT","t");

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
        if (output.endsWith(File.separator)) {
            output+="recombinants.fasta";
        }
        System.out.println(sb.toString());
//        Utils.saveFile(output, sb.toString());
    }
}
