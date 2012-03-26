/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.distance;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import org.javatuples.Pair;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class DistanceUtils {

    public static enum PhiCalcMethod {

        META_RECOMBINANT,
        ONLY_GIVEN_HAPLOTYPES
    }

    public static Pair[] calculatePhi(Map<String, String> PH, Map<String, Integer> PHtest) {
        return calculatePhi(PH, PHtest, PhiCalcMethod.ONLY_GIVEN_HAPLOTYPES);
    }

    public static Pair[] calculatePhi(Map<String, String> original, Map<String, Integer> ph, PhiCalcMethod method) {
        int max = 0;
        for (Integer x : ph.values()) {
            max += x;
        }
        int N = 20;
        Pair[] pairs = new Pair[N];
        double[] dist = new double[N];
        double[] recombs = new double[N];
        for (int DELTA = 0; DELTA < N; DELTA++) {

            Map<String, Double> distMap = new HashMap<>();
            for (String s : original.keySet()) {
                distMap.put(s, 0d);
            }
            double recombinant = 0d;
            Map<String, List<Pair<String, Double>>> PHtoPHtest = new HashMap<>();
            for (String ptest : ph.keySet()) {
                List<Pair<String, Double>> list = new LinkedList<>();
                for (String ptrue : original.keySet()) {
                    double error = 0d;
                    error = calcHamming(ptrue, ptest);
                    if (error <= DELTA) {
                        list.add(Pair.with(ptrue, error));
                    }
                }
                if (list.isEmpty()) {
                    if (method == PhiCalcMethod.META_RECOMBINANT) {
                        recombs[DELTA]++;
                    }
                }
                PHtoPHtest.put(ptest, list);
                double minError = Integer.MAX_VALUE;
                List<Pair<String, Integer>> best = new LinkedList<>();
                for (Pair p : list) {
                    if ((double) p.getValue1() < minError) {
                        minError = (double) p.getValue1();
                        best.clear();
                        best.add(p);
                    } else if ((double) p.getValue1() == minError) {
                        best.add(p);
                    }
                }
                for (Pair p : best) {
                    String ptrue = (String) p.getValue0();
                    distMap.put(ptrue, distMap.get(ptrue) + ph.get(ptest) / (double) best.size());
                }
            }
            double sum = 0d;
            for (String hap : original.keySet()) {
                sum += distMap.get(hap);
            }
            dist[DELTA] = sum / (double) max;
            pairs[DELTA] = Pair.with(dist[DELTA], recombs[DELTA]/(double)max);
        }
        return pairs;
    }

    public static int calcHamming(String ptrue, String ptest) {
        int error = 0;
        char[] p1 = ptrue.toCharArray();
        char[] p2 = ptest.toCharArray();
        for (int i = 0; i < ptest.length(); i++) {
            if (p1[i] != p2[i]) {
                error++;
            }
        }
        return error;
    }

    public static Double calculateKLD(Map<String, Integer> PH, Map<String, Integer> PhatH) {

        for (String s : PH.keySet()) {
            if (!PhatH.containsKey(s)) {
                PhatH.put(s, 1);
            }
        }
        for (String s : PhatH.keySet()) {
            if (!PH.containsKey(s)) {
                PH.put(s, 1);
            }
        }

        int PHsize = 0;
        for (Integer i : PH.values()) {
            PHsize += i;
        }
        int PhatHsize = 0;
        for (Integer i : PhatH.values()) {
            PhatHsize += i;
        }

        Double result = 0.0;
        Double frequency2 = 0.0;
        StringBuilder sb = new StringBuilder();
        for (String sequence : PH.keySet()) {
            sb.append(sequence).append("\t").append(PH.get(sequence)).append("\t").append(PhatH.get(sequence)).append("\n");
            Double frequency1 = (double) PH.get(sequence) / PHsize;
            frequency2 = (double) PhatH.get(sequence) / PhatHsize;
            result += frequency1 * (Math.log(frequency1 / frequency2) / Math.log(2));
        }
        System.out.println(sb);
        return result;
    }
    
    public static Double calculateKLD2(Map<String, Integer> P, Map<String, Integer> Q) {

        double Psize = (double) P.values().size();
        double Qsize = (double) Q.values().size();

        double result = 0.0;
        for (String sequence : P.keySet()) {
            double fP = P.get(sequence) / Psize;
            double fQ = Q.containsKey(sequence) ? Q.get(sequence) / Qsize : 0d;
            fQ += fP;
            fQ /= 2;
            result += fP * (Math.log(fP / fQ) / Math.log(2));
        }
        return result;
    }
}
