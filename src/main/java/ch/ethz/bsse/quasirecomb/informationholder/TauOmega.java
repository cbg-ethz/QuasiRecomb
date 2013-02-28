///**
// * Copyright (c) 2011-2013 Armin Töpfer
// *
// * This file is part of QuasiRecomb.
// *
// * QuasiRecomb is free software: you can redistribute it and/or modify it under
// * the terms of the GNU General Public License as published by the Free Software
// * Foundation, either version 3 of the License, or any later version.
// *
// * QuasiRecomb is distributed in the hope that it will be useful, but WITHOUT
// * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// * details.
// *
// * You should have received a copy of the GNU General Public License along with
// * QuasiRecomb. If not, see <http://www.gnu.org/licenses/>.
// */
//package ch.ethz.bsse.quasirecomb.informationholder;
//
//import ch.ethz.bsse.quasirecomb.utils.Utils;
//import java.io.File;
//import java.util.HashMap;
//import java.util.Map;
//
///**
// * @author Armin Töpfer (armin.toepfer [at] gmail.com)
// */
//public class TauOmega {
////    private Map<Integer, Map<Integer,Double>> tauWatsonMap = new HashMap<>();
////    private Map<Integer, Map<Integer,Double>> omegaWatsonMap = new HashMap<>();
////    private Map<Integer, Map<Integer,Double>> tauCrickMap = new HashMap<>();
////    private Map<Integer, Map<Integer,Double>> omegaCrickMap = new HashMap<>();
////    private Map<Integer, Map<Integer,Double>> deletionMap = new HashMap<>();
//
//    private Map<Integer, Double> tauWatsonMap = new HashMap<>();
//    private Map<Integer, Double> omegaWatsonMap = new HashMap<>();
//    private Map<Integer, Double> tauCrickMap = new HashMap<>();
//    private Map<Integer, Double> omegaCrickMap = new HashMap<>();
//    private Map<Integer, Map<Integer, Double>> deletionMap = new HashMap<>();
//    private int[] coverage;
//
//    private void init(Read[] reads, int L) {
//        double N = 0;
//        for (Read r : reads) {
//            N += r.getCount();
//        }
//        this.coverage = new int[L];
//        for (Read r : reads) {
//            for (int i = r.getWatsonBegin(); i < r.getWatsonEnd(); i++) {
//                this.coverage[i] += r.getCount();
//            }
//            setEntry(tauWatsonMap,r.getWatsonBegin(), r.getCount()/N);
//            setEntry(omegaWatsonMap,r.getWatsonEnd(), r.getCount()/N);
//            if (r.isPaired()) {
//                for (int i = r.getCrickBegin(); i < r.getCrickEnd(); i++) {
//                    this.coverage[i] += r.getCount();
//                }
//                tau2[r.getCrickBegin()] += r.getCount();
//                omega2[r.getCrickEnd()] += r.getCount();
//            }
//        }
//    }
//
//    private void setEntry(Map<Integer, Double> map, int j, double f) {
//        map.put(j, map.containsKey(j) ? map.get(j) + f : f);
//    }
////    private void setEntry(Map<Integer, Map<Integer, Double>> map, int j, int x, double f) {
////        if (!map.containsKey(j)) {
////            map.put(j, new HashMap<Integer, Double>());
////        }
////        map.get(j).put(f / N)
////    }
//}
