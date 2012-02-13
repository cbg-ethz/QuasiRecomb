package ch.ethz.bsse.quasirecomb.utils;

import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

/**
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class Frequency<T> {

    private final NavigableMap<Double, T> table = new TreeMap<>();
    private final double max;

    public Frequency(Map<T, Double> frequency) {
        double total = 0;
        
        for (Map.Entry<T, Double> e : frequency.entrySet()) {
            if (e.getValue() > 1E-5) {
                total += e.getValue();
                table.put(total, e.getKey());
            }
        }
        max = total;
    }

    /** 
     * Choose a random symbol. The choices are weighted by frequency.
     */
    public T roll() {
        Double key = Math.random() * max;
//        Double key = Math.random() * max;
//        Double key = generator.nextDouble() * max;
        return table.higherEntry(key).getValue();
    }
}
