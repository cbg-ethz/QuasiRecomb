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

import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class TauOmega implements Serializable {

    private Map<Integer, Map<Integer, Double>> tauWatsonMap = new HashMap<>();
    private Map<Integer, Map<Integer, Double>> omegaWatsonMap = new HashMap<>();
    private Map<Integer, Map<Integer, Double>> tauCrickMap = new HashMap<>();
//    private Map<Integer, Map<Integer, Double>> deletionMap = new HashMap<>();
    private int[] coverage;

    public TauOmega(Read[] reads, int L) {
        this.init(reads, L);
    }

    private void init(Read[] reads, int L) {
        double N = Globals.getINSTANCE().getNREAL();
        this.coverage = new int[L];
        for (Read r : reads) {
            for (int i = r.getWatsonBegin(); i < r.getWatsonEnd(); i++) {
                this.coverage[i] += r.getCount();
            }
            setEntry(tauWatsonMap, r.getWatsonBegin(), r.getWatsonLength(), r.getCount() / N);
            if (r.isPaired()) {
                setEntry(omegaWatsonMap, r.getWatsonEnd(), r.getInsertSize(), r.getCount() / N);
                for (int i = r.getCrickBegin(); i < r.getCrickEnd(); i++) {
                    this.coverage[i] += r.getCount();
                }
                setEntry(tauCrickMap, r.getCrickBegin(), r.getCrickLength(), r.getCount() / N);
            }
        }

        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < L; j++) {
            sb.append(getPostitionProb(this.tauWatsonMap, j)).append("\t");
            sb.append(getPostitionProb(this.omegaWatsonMap, j)).append("\t");
            sb.append(getPostitionProb(this.tauCrickMap, j)).append("\n");
        }

        if (Globals.getINSTANCE().isDEBUG()) {
            Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "support" + File.separator + "twtw", sb.toString());
        }
    }

    private double getPostitionProb(Map<Integer, Map<Integer, Double>> map, int j) {
        if (map.containsKey(j)) {
            double sum = 0;
            for (double d : map.get(j).values()) {
                sum += d;
            }
            return sum;
        } else {
            return 0d;
        }
    }

    private void setEntry(Map<Integer, Map<Integer, Double>> map, int s, int l, double f) {
        if (!map.containsKey(s)) {
            map.put(s, new HashMap<Integer, Double>());
        }
        map.get(s).put(l, f);
    }

    public double getTauWatsonProbability(int s, int l) {
        return this.tauWatsonMap.get(s).get(l);
    }

    public double getTauCrickProbability(int s, int l) {
        return this.tauCrickMap.get(s).get(l);
    }

    public double getOmegaWatsonProbability(int s, int l) {
        return this.omegaWatsonMap.get(s).get(l);
    }

    public Integer[] getTauWatsonLengths(int s) {
        return this.tauWatsonMap.get(s).keySet().toArray(new Integer[this.tauWatsonMap.get(s).keySet().size()]);
    }

    public Integer[] getTauCrickLengths(int s) {
        return this.tauCrickMap.get(s).keySet().toArray(new Integer[this.tauCrickMap.get(s).keySet().size()]);
    }

    public Integer[] getOmegaWatsonLengths(int s) {
        return this.omegaWatsonMap.get(s).keySet().toArray(new Integer[this.omegaWatsonMap.get(s).keySet().size()]);
    }

    public Double[] getTauWatsonProbabilities(int s) {
        return this.tauWatsonMap.get(s).values().toArray(new Double[this.tauWatsonMap.get(s).values().size()]);
    }

    public Double[] getTauCrickProbabilities(int s) {
        return this.tauCrickMap.get(s).values().toArray(new Double[this.tauCrickMap.get(s).values().size()]);
    }

    public Double[] getOmegaWatsonProbabilities(int s) {
        return this.omegaWatsonMap.get(s).values().toArray(new Double[this.omegaWatsonMap.get(s).values().size()]);
    }

    public Map<Integer, Map<Integer, Double>> getTauWatsonMap() {
        return tauWatsonMap;
    }

    public Map<Integer, Map<Integer, Double>> getOmegaWatsonMap() {
        return omegaWatsonMap;
    }

    public Map<Integer, Map<Integer, Double>> getTauCrickMap() {
        return tauCrickMap;
    }

    public int[] getCoverage() {
        return coverage;
    }

    public boolean isPaired() {
        return !tauCrickMap.isEmpty();
    }
}
