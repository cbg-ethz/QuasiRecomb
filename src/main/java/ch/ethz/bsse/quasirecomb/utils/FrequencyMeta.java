/**
 * Copyright (c) 2011-2013 Armin Töpfer
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
package ch.ethz.bsse.quasirecomb.utils;

import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class FrequencyMeta {

    private final NavigableMap<Double, Integer> table = new TreeMap<>();
    private final double max;

    public FrequencyMeta(Map<Integer, Double> frequency) {
        double total = 0;

        for (Map.Entry<Integer, Double> e : frequency.entrySet()) {
            if (e.getValue() > 1E-6) {
                total += e.getValue();
                table.put(total, e.getKey());
            }
        }
        table.put(1d, -1);
        max = 1;
    }

    /**
     * Choose a random symbol. The choices are weighted by frequency.
     */
    public Integer roll() {
        try {
            Double key = Math.random() * max;
            return table.higherEntry(key).getValue();
        } catch (Exception e) {
            System.out.println("error rolling");
            return null;
        }
    }
}
