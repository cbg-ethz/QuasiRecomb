/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.distance;

import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.utils.FastaParser;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

/**
 *
 * @author xlr
 */
public class IntersectQuasispecies {

    public static StringBuilder compute(String path1, String path2) {
        Map<String, Double> q1 = FastaParser.parseQuasispeciesFile(path1);
        Map<String, Double> q2 = FastaParser.parseQuasispeciesFile(path2);
        Map<String, Double> intersect = compute(q1, q2);

        StringBuilder sb = new StringBuilder();
        int i = 0;
        for (Object o : ModelSampling.sortMapByValue(intersect).entrySet()) {
            Entry<String, Double> e = (Entry) o;
            sb.append(">INTERSECT").append(i++).append("_").append(shorten(e.getValue())).append("\n");
            sb.append(e.getKey()).append("\n");
        }
        return sb;
    }

    public static Map<String, Double> compute(Map<String, Double> q1, Map<String, Double> q2) {
        Map<String, Double> intersect = new HashMap<>();
        double sum = 0d;
        for (String s1 : q1.keySet()) {
            if (q2.containsKey(s1)) {
                intersect.put(s1, q1.get(s1) + q2.get(s1));
                sum += q1.get(s1) + q2.get(s1);
            }
        }
        StringBuilder sb = new StringBuilder();
        int i = 0;
        for (Object o : ModelSampling.sortMapByValue(intersect).entrySet()) {
            Entry<String, Double> e = (Entry) o;
            intersect.put(e.getKey(), e.getValue() / sum);
        }
        return intersect;
    }

    public static String shorten(double value) {
        String s;
        if (value < 1e-20) {
            s = "0      ";
        } else if (value == 1.0) {
            s = "1      ";
        } else {
            String t = "" + value;
            String r;
            if (t.length() > 7) {
                r = t.substring(0, 7);
                if (t.contains("E")) {
                    r = r.substring(0, 4);
                    r += "E" + t.split("E")[1];
                }
                s = r;
            } else {
                s = String.valueOf(value);
            }
        }
        return s;
    }
}
