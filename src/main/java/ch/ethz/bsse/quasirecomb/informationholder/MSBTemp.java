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
package ch.ethz.bsse.quasirecomb.informationholder;

import com.google.common.collect.Multimap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math.stat.descriptive.rank.Median;


/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class MSBTemp {

    private Map<Integer, SelectionResultBootstrap> srMap = new HashMap<>();
    private int bestK;

    public MSBTemp(Multimap<Integer, Double> map) {
        for (Map.Entry<Integer, Collection<Double>> e : map.asMap().entrySet()) {
            this.add(e.getValue(), e.getKey());
        }
    }

    public void add(Collection<Double> bics, int K) {
        srMap.put(K, new SelectionResultBootstrap(bics));
    }

    private void select() {
        double maxValue = Double.NEGATIVE_INFINITY;
        int maxK = 0;
        for (Map.Entry<Integer, SelectionResultBootstrap> entry : srMap.entrySet()) {
            SelectionResultBootstrap srTmp = entry.getValue();
            if (maxValue < srTmp.median) {
                maxValue = srTmp.median;
                maxK = entry.getKey();
            }
        }
        if (maxK == 1) {
            this.bestK = maxK;
        } else {
            while (maxK - 1 > 0) {
                if (srMap.containsKey(maxK - 1)) {
                    SelectionResultBootstrap previous = srMap.get(maxK - 1);
                    if (srMap.get(maxK).median < previous.lowerBound) {
                        maxK--;
                    } else {
                        this.bestK = maxK;
                        break;
                    }
                }
            }
        }
    }

    public int getBestK() {
        this.select();
        return this.bestK;
    }

    public Map<Integer, SelectionResultBootstrap> getSrMap() {
        return srMap;
    }
    
}
class SelectionResultBootstrap {

    double median;
    double lowerBound;

    public SelectionResultBootstrap(Collection<Double> bics) {
        List<Double> list = new ArrayList<>(bics);
        double[] bicsTmp = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            bicsTmp[i] = list.get(i);
        }
        this.median = new Median().evaluate(bicsTmp);
        this.lowerBound = median - new StandardDeviation().evaluate(bicsTmp) * Math.sqrt(1 + 1d / bicsTmp.length);
    }
}