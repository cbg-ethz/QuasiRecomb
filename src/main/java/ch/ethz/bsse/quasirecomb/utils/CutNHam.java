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

import ch.ethz.bsse.quasirecomb.distance.HammerWorker;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class CutNHam {

    private Map<String, String> seqs = new HashMap<>();

    public CutNHam(String input, int from, int to, int windowSize) {
        StringBuilder sb = new StringBuilder();
        this.read(input);
        if (windowSize == 0) {
//            Triplet<Integer, String, String> min = cutandham(from, to);
//            System.out.println("");
//            System.out.println(min);
        } else {
            for (int i = from; i < to; i += 1) {
                int b = i + windowSize;
                if (b > to) {
                    b = to;
                    break;
                }
                Pair<Integer, List<Pair<String, String>>> cutandham = cutandham(i, b);

                sb.append(i).append("\t").append(cutandham.getValue0()).append("\t");
                List<String> strains = new LinkedList<>();
                for (Pair<String, String> p : cutandham.getValue1()) {
                    if (!strains.contains(p.getValue0())) {
                        strains.add(p.getValue0());
                    }
                    if (!strains.contains(p.getValue1())) {
                        strains.add(p.getValue1());
                    }
                }
                for (String s : new String[]{"HXB2", "NL4-3", "YU-2", "89.6", "JRCSF"}) {
                    if (strains.contains(s)) {
                        sb.append(s);
                    } else {
                        sb.append(" ");
                    }
                    sb.append("\t");

                }
                sb.append("\n");
            }
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "distdist_" + windowSize + ".txt", sb.toString());
    }

    private Map<String, String> cut(int from, int to) {
        Map<String, String> cuts = new HashMap<>();
        for (Map.Entry<String, String> e : seqs.entrySet()) {
            cuts.put(e.getKey(), e.getValue().substring(from - 1, to));
        }
        return cuts;
    }

    private void read(String input) {
        try {
            FileInputStream fstream = new FileInputStream(input);
            StringBuilder sb;
            String head = null;
            int i = 0;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();

                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (Globals.getINSTANCE().isDEBUG()) {
                            i++;
                            if (i % 100 == 0) {
                                System.out.println(i);
                            }
                        }
                        if (sb.length() > 0) {
                            seqs.put(head, sb.toString());
                            sb.setLength(0);
                        }
                        head = strLine;
                    } else {
                        sb.append(strLine);
                    }
                }
                seqs.put(head, sb.toString());
            }
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
    }

    private Pair<Integer, List<Pair<String, String>>> cutandham(int from, int to) {
        Map<String, String> cut = this.cut(from, to);
        String[] sequences = new String[cut.size()];
        Map<String, String> identify = new HashMap<>();
        int i = 0;
        for (Map.Entry<String, String> e : cut.entrySet()) {
            sequences[i++] = e.getValue();
            identify.put(e.getValue(), e.getKey().replaceAll(">", ""));
        }
        Globals.getINSTANCE().setHammingMax((int) Math.pow(sequences.length, 2));
        Map<String, Map<String, Integer>> invoke = Globals.getINSTANCE().getFjPool().invoke(new HammerWorker(sequences, 0, sequences.length));
        int min = Integer.MAX_VALUE;
        List<Pair<String, String>> mins = new LinkedList<>();
        if (invoke.size() == 1) {
            mins.add(Pair.with("all", "all"));
            return Pair.with(0, mins);
        }
        if (from == 3001) {
            System.out.println("");
        }
        for (Map.Entry<String, Map<String, Integer>> entry : invoke.entrySet()) {
            for (Map.Entry<String, Integer> inner : entry.getValue().entrySet()) {
                if (inner.getValue() < min) {
                    min = inner.getValue();
                    mins = new LinkedList<>();
                    String a = identify.get(entry.getKey());
                    String b = identify.get(inner.getKey());
                    mins.add(Pair.with(a, b));
                } else if (inner.getValue() == min) {
                    String a = identify.get(entry.getKey());
                    String b = identify.get(inner.getKey());
                    if (!mins.contains(Pair.with(b, a))) {
                        mins.add(Pair.with(a, b));
                    }
                }
            }
        }
        return Pair.with(min, mins);
    }
}
