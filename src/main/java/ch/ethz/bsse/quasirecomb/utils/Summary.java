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

import ch.ethz.bsse.quasirecomb.distance.KullbackLeibler;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import org.javatuples.Pair;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Summary extends Utils {

    public String snvs(OptimalResult or) {
        StringBuilder sb = new StringBuilder();
        char[] alphabet = new char[]{'A', 'C', 'G', 'T', '-'};

        for (int v = 0; v < or.getn(); v++) {
            sb.append("\t").append(alphabet[v]);
        }
        sb.append("\n");
        for (int j = 0; j < or.getL(); j++) {
            sb.append(j);
            double sum = 0d;
            for (int v = 0; v < or.getn(); v++) {
                sum += or.getSnv()[j][v]>1e-5?or.getSnv()[j][v]:0;
            }
            for (int v = 0; v < or.getn(); v++) {
                sb.append("\t").append(shortenSmall(or.getSnv()[j][v]/sum));
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    
    public static String shortenSmall(double value) {
        String s;
        if (value < 1e-5) {
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

    private String opacity(double value) {
        String s;
        if (value < 1e-3) {
            s = "0";
        } else if (value == 1.0) {
            s = "1";
        } else {
            s = "" + value;
        }
        return s;
    }

    public void printAlignment(Read[] reads) {
        Map<Integer, List<Read>> readMap = new HashMap<>();
        for (Read r : reads) {
            if (!readMap.containsKey(r.getBegin())) {
                readMap.put(r.getBegin(), new ArrayList<Read>());
            }
            readMap.get(r.getBegin()).add(r);
        }
        Integer[] order = readMap.keySet().toArray(new Integer[readMap.size()]);
        Arrays.sort(order);
        int min = order[0];
        StringBuilder sb = new StringBuilder();
        for (int item : order) {
            List<Read> currentReads = readMap.get(item);
            for (Read currentRead : currentReads) {
                for (int i = min; i < item; i++) {
                    sb.append(" ");
                }
                String s = BitMagic.toString(currentRead.getSequence());
                for (char c : s.toCharArray()) {
                    sb.append(reverse(Integer.parseInt("" + c)));
                }
                sb.append("\n");
            }
        }
        Utils.saveFile(Globals.getINSTANCE().getSAVEPATH() + "alignment.txt", sb.toString());
    }

    public void circos(int L, int[][] alignment) {
        String[] genome = FastaParser.parseFarFile(Globals.getINSTANCE().getGENOME());
        //FastaParser.parseFarFile("/Users/XLR/Dropbox/Projects/karin/PacBio/run1/wglobal/ref_5343-8525.fasta");
        int[][][] alignmentG = new int[genome.length][L][5];
        double[][] alignmentE = new double[L][5];
        for (int i = 0; i < genome.length; i++) {
            byte[] r = splitReadIntoByteArray(genome[i]);
            for (int j = 0; j < genome[0].length(); j++) {
                alignmentG[i][j][r[j]]++;
                alignmentE[j][r[j]] += 1d / genome.length;
            }
        }
        int i = 0;
//        for (int i = 0; i < L; i += 500) 
        {
            int from = i;
            int to = L;
//            int to = from + 500;
//            if (to >= L) {
//                to = L;
//            }
            StringBuilder dataSB = new StringBuilder();
            StringBuilder rawEntropyFill = new StringBuilder();
            StringBuilder[] genomeEntropyFill = new StringBuilder[genome.length];
            for (int x = 0; x < genome.length; x++) {
                genomeEntropyFill[x] = new StringBuilder();
            }
            StringBuilder ticks = new StringBuilder();
            ticks.append("show_ticks          = yes\n"
                    + "show_tick_labels    = yes\n"
                    + "\n"
                    + "<ticks>\n"
                    + "offset = - 50p\n"
                    + "radius         = dims(ideogram,radius_outer)\n"
                    + "multiplier     = 1/1u\n"
                    + "color          = black\n"
                    + "force_display  = yes\n"
                    + "thickness      = 2p\n"
                    + "show_label     = yes\n"
                    + "format         = %d\n");
            StringBuilder coverage = new StringBuilder();
            int h = 0;
            int g = 0;
            boolean nextTick = false;
            for (int j = from; j < to; j++) {
                double max = 0d;
                for (int v = 0; v < 5; v++) {
                    max = Math.max(max, alignmentE[j][v]);
                }
                if ((5433 + j) % 250 == 0) {
                    nextTick = true;
                }
                if (max < 1d) {
                    if (nextTick) {
                        nextTick = false;
                        ticks.append("<tick>\n"
                                + "show           = yes\n"
                                + "position       = " + g + "u\n"
                                + "size           = 15p\n"
                                + "label_size     = 36p\n"
                                + "label          = " + (j + 5433) + "\n"
                                + "</tick>\n");
                    }
                    double sum = 0d;
                    double sumG[] = new double[genome.length];
                    for (int v = 0; v < alignment[j].length; v++) {
                        sum += alignment[j][v];
                        for (int x = 0; x < genome.length; x++) {
                            sumG[x] += alignmentG[x][j][v];
                        }
                    }
                    rawEntropyFill.append("h0 ").append(g).append(" ").append(g + 1).append(" ").append(alignment[j][0] / sum);
                    for (int x = 0; x < genome.length; x++) {
                        genomeEntropyFill[x].append("h0 ").append(g).append(" ").append(g + 1).append(" ").append(alignmentG[x][j][0] / sumG[x]);
                    }
                    for (int v = 1; v < alignment[j].length; v++) {
                        rawEntropyFill.append(",").append(alignment[j][v] / sum);
                        for (int x = 0; x < genome.length; x++) {
                            genomeEntropyFill[x].append(",").append(alignmentG[x][j][v] / sumG[x]);
                        }
                    }
                    rawEntropyFill.append("\n");
                    for (int x = 0; x < genome.length; x++) {
                        genomeEntropyFill[x].append("\n");
                    }
                    coverage.append("h0 ").append(g).append(" ").append(g + 1).append(" ").append(sum).append("\n");
                    g++;
                }
            }
            ticks.append("</ticks>\n");
            dataSB.append("chr - h0 1 ").append(from).append(" ").append(g - 1).append(" blue\n");
            new File("/Users/XLR/Dropbox/basicPlot/" + i + "/support/").mkdirs();
            Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/ticks.conf", ticks.toString());
            Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/data.txt", dataSB.toString());
            Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/rawEntropy.txt", rawEntropyFill.toString());
            for (int x = 0; x < genome.length; x++) {
                Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/genomeEntropy" + x + ".txt", genomeEntropyFill[x].toString());
            }
//            Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/rawEntropy1.txt", rawEntropy1.toString());
//            Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/rawEntropy2.txt", rawEntropy2.toString());
//            Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/rawEntropy3.txt", rawEntropy3.toString());
//            Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/rawEntropy4.txt", rawEntropy4.toString());
            Utils.saveFile("/Users/XLR/Dropbox/basicPlot/" + i + "/support/coverage.txt", coverage.toString());
        }

//        Utils.saveFile(path + "rho_" + j, sb.toString());
    }

    public String kl(OptimalResult or) {
        StringBuilder sb = new StringBuilder();
        double[][][] mu = or.getMu();
        Set<Pair<Integer, Integer>> a = new HashSet<>();
        for (int k = 0; k < or.getK(); k++) {
            for (int l = 0; l < or.getK(); l++) {
                if (k != l) {
                    if (!a.contains(Pair.with(k, l)) && !a.contains(Pair.with(l, k))) {
                        a.add(Pair.with(k, l));
                        sb.append(k).append(" <-> ").append(l).append(" = ").append(KullbackLeibler.symmetric(mu, k, l)).append("\n");
                    }
                }
            }
        }
        return sb.toString();
    }

    public String minimal(OptimalResult or) {
        StringBuilder sb = new StringBuilder();
        sb.setLength(0);
        sb.append("#loglikelihood:").append(or.getLlh()).append("\n");
        sb.append("#BIC:").append(or.getBIC()).append("\n");
        double[][][] mu = or.getMu();
        double[][][] rho = or.getRho();
        sb.append("\nPI:\t\t");
        for (int k = 0; k < or.getK(); k++) {
            sb.append(shorten(or.getPi()[k]));
            sb.append("\t\t\t\t\t");
        }
        sb.append("\n\n");
        StringBuilder sb2 = new StringBuilder();
        sb2.append("j\teps\t");
        for (int k = 0; k < or.getK(); k++) {
            sb2.append("Generator ").append(k).append("\t\t\t\t");
        }
        sb.append("---");
        for (int i = 0; i < sb2.toString().length(); i++) {
            sb.append("-");
        }
        for (int k = 0; k < or.getK() - 1; k++) {
            sb.append("-----------------------------");
        }
        sb.append("\n");
        sb.append(sb2);
        sb.append("\n");
        sb.append("---");
        for (int i = 0; i < sb2.length(); i++) {
            sb.append("-");
        }
        for (int k = 0; k < or.getK() - 1; k++) {
            sb.append("-----------------------------");
        }
        sb.append("\n");
        for (int j = 0; j < or.getL(); j++) {


            //mu
            boolean flatMuK = false;
            for (int k = 0; k < or.getK(); k++) {
                for (int v = 0; v < or.getMu()[0][0].length; v++) {
                    if (or.getMu()[j][k][v] > 1e-20 && or.getMu()[j][k][v] != 1d) {
                        flatMuK = true;
                        break;
                    }
                }
            }
            boolean flatRho = false;
            if (j < or.getL() - 1) {
                for (int k = 0; k < or.getK(); k++) {
                    for (int l = 0; l < or.getK(); l++) {
                        if (or.getRho()[j][k][l] > 1e-20 && or.getRho()[j][k][l] != 1d) {
                            flatRho = true;
                            break;
                        }
                        if (or.getRho()[j][k][l] == 1d && k != l) {
                            flatRho = true;
                            break;
                        }
                    }
                    if (flatRho) {
                        break;
                    }
                }
            }
            if (flatRho || flatMuK) {
                sb.append(j + 1).append("\t");
                sb.append(shorten(or.getEps()[j])).append("\t");
            }
            if (flatMuK) {
                for (int k = 0; k < or.getK(); k++) {
                    boolean flatMu = false;
                    for (int v = 0; v < or.getMu()[0][0].length; v++) {
                        if (or.getMu()[j][k][v] > 1e-20 && or.getMu()[j][k][v] != 1d) {
                            flatMu = true;
                            break;
                        }
                    }
                    if (flatMu) {
                        sb.append("[");
                        for (int v = 0; v < or.getn(); v++) {
                            sb.append(shorten(or.getMu()[j][k][v]));
                            if (v + 1 < or.getn()) {
                                sb.append(", ");
                            }
                        }
                        sb.append("]\t");
                    } else {
                        double max = Double.MIN_VALUE;
                        Map<Double, Integer> m = new ConcurrentHashMap<>();
                        for (int v = 0; v < or.getMu()[0][0].length; v++) {
                            max = Math.max(max, or.getMu()[j][k][v]);
                            m.put(or.getMu()[j][k][v], v);
                        }
                        if (max == Double.MIN_VALUE) {
                            sb.append(0).append("\t\t\t\t\t");
                        } else {
                            sb.append(reverse(m.get(max))).append("\t\t\t\t\t");
                        }
                    }
                }
            }
            if (flatRho || flatMuK) {
                sb.append("\n");
            }
            //recombination
            if (j < or.getL() - 1) {
//                boolean flatRho = false;
//                for (int k = 0; k < or.getK(); k++) {
//                    for (int l = 0; l < or.getK(); l++) {
//                        if (or.getRho()[j][k][l] > 1e-20 && or.getRho()[j][k][l] != 1d) {
//                            flatRho = true;
//                            break;
//                        }
//                        if (or.getRho()[j][k][l] == 1d && k != l) {
//                            flatRho = true;
//                            break;
//                        }
//                    }
//                    if (flatRho) {
//                        break;
//                    }
//                }
                if (flatRho) {
                    sb.append("RECOMBINATION").append("\t");
                    for (int k = 0; k < or.getK(); k++) {
                        sb.append("[");
                        for (int l = 0; l < or.getK(); l++) {
                            sb.append(shorten(or.getRho()[j][k][l]));
                            if (l + 1 < or.getK()) {
                                sb.append(", ");
                            }

                        }
                        sb.append("]\t\t\t");
                    }
                    sb.append("\n");
                }
            }

        }
        return sb.toString();
    }

    public String print(OptimalResult or) {
        StringBuilder sb = new StringBuilder();
        sb.setLength(0);
        sb.append("#loglikelihood:").append(or.getLlh()).append("\n");
        sb.append("#BIC:").append(or.getBIC()).append("\n");
        double[][][] mu = or.getMu();
//        double mue = 0d;
//        for (int j = 0; j < or.getL(); j++) {
//            for (int k = 0; k < or.getK(); k++) {
//                for (int v = 0; v < or.getn(); v++) {
//                    mue -= mu[j][k][v] * Math.log(mu[j][k][v]) / Math.log(or.getn());
//                }
//            }
//        }
//        sb.append("#MUE:").append(mue / (or.getK() * or.getL())).append("\n");
//        double rhoe = 0d;
        double[][][] rho = or.getRho();
//        for (int j = 0; j < or.getL() - 1; j++) {
//            for (int k = 0; k < or.getK(); k++) {
//                for (int l = 0; l < or.getK(); l++) {
//                    rhoe -= rho[j][k][l] * Math.log(rho[j][k][l]) / Math.log(or.getK());
//                }
//            }
//        }
//        sb.append("#RHOE:").append(mue / (or.getK() * or.getL())).append("\n");
//        sb.append("#PE:").append((mue - rhoe) / (or.getK() * or.getL())).append("\n");
        sb.append("\nPI:\t\t");
        for (int k = 0; k < or.getK(); k++) {
            sb.append(shorten(or.getPi()[k]));
            sb.append("\t\t\t\t\t");
        }
        sb.append("\n\n");
        StringBuilder sb2 = new StringBuilder();
        sb2.append("j\teps\t");
        for (int k = 0; k < or.getK(); k++) {
            sb2.append("Generator ").append(k).append("\t\t\t\t");
        }
        sb.append("---");
        for (int i = 0; i < sb2.toString().length(); i++) {
            sb.append("-");
        }
        for (int k = 0; k < or.getK() - 1; k++) {
            sb.append("-----------------------------");
        }
        sb.append("\n");
        sb.append(sb2);
        sb.append("\n");
        sb.append("---");
        for (int i = 0; i < sb2.length(); i++) {
            sb.append("-");
        }
        for (int k = 0; k < or.getK() - 1; k++) {
            sb.append("-----------------------------");
        }
        sb.append("\n");
        for (int j = 0; j < or.getL(); j++) {

            sb.append(j + 1).append("\t");
            sb.append(shorten(or.getEps()[j])).append("\t");
            //mu
            for (int k = 0; k < or.getK(); k++) {
                boolean flatMu = false;
                for (int v = 0; v < or.getMu()[0][0].length; v++) {
                    if (or.getMu()[j][k][v] > 1e-20 && or.getMu()[j][k][v] != 1d) {
                        flatMu = true;
                        break;
                    }
                }
                if (flatMu) {
                    sb.append("[");
                    for (int v = 0; v < or.getn(); v++) {
//                        if (or.getMu()[j][k][v] > 0) {
//                            entropy -= or.getMu()[j][k][v] * Math.log(or.getMu()[j][k][v]);
//                        }
                        sb.append(shorten(or.getMu()[j][k][v]));
                        if (v + 1 < or.getn()) {
                            sb.append(", ");
                        }
                    }
                    sb.append("]\t");
                } else {
                    double max = Double.MIN_VALUE;
                    Map<Double, Integer> m = new ConcurrentHashMap<>();
                    for (int v = 0; v < or.getMu()[0][0].length; v++) {
                        max = Math.max(max, or.getMu()[j][k][v]);
                        m.put(or.getMu()[j][k][v], v);
                    }
                    if (max == Double.MIN_VALUE) {
//                        System.out.println(Arrays.toString(m.keySet().toArray(new Double[m.size()])));
                        sb.append(0).append("\t\t\t\t\t");
                    } else {
                        sb.append(reverse(m.get(max))).append("\t\t\t\t\t");
                    }
                }
            }
            sb.append("\n");
            //recombination
            if (j < or.getL() - 1) {
                boolean flat = false;
                for (int k = 0; k < or.getK(); k++) {
                    for (int l = 0; l < or.getK(); l++) {
                        if (or.getRho()[j][k][l] > 1e-20 && or.getRho()[j][k][l] != 1d) {
                            flat = true;
                            break;
                        }
                        if (or.getRho()[j][k][l] == 1d && k != l) {
                            flat = true;
                            break;
                        }
                    }
                    if (flat) {
                        break;
                    }
                }
                if (flat) {
                    sb.append("RECOMBINATION").append("\t");
                    for (int k = 0; k < or.getK(); k++) {
                        sb.append("[");
                        for (int l = 0; l < or.getK(); l++) {
                            sb.append(shorten(or.getRho()[j][k][l]));
                            if (l + 1 < or.getK()) {
                                sb.append(", ");
                            }

                        }
                        sb.append("]\t\t\t");
                    }
                    sb.append("\n");
                }
            }

        }
        return sb.toString();
    }

    public String html(OptimalResult or) {
        List<Integer> breaks = new ArrayList<>();
        for (int j = 0; j < or.getL(); j++) {
            for (int k = 0; k < or.getK(); k++) {
                boolean flat = false;
                if (j < or.getL() - 1) {
                    for (int l = 0; l < or.getK(); l++) {
//                        if (or.getRho()[j][k][l] > 1e-20 && or.getRho()[j][k][l] != 1d) {
                        if (k == l && or.getRho()[j][k][l] != 1d) {
                            flat = true;
                            break;
                        }
                    }
                }
                if (flat) {
                    if (!breaks.contains(j)) {
                        breaks.add(j);
                    }
                    break;
                }
            }
        }


        StringBuilder sb = new StringBuilder();
        String A = "#00A1CB";
        String C = "#61AE24";
        String G = "#E54028";
        String T = "#D0D102";
        String gap = "#616161";
        sb.append("<html>");
        sb.append("<style>"
                + "* {font-size: 10pt}\n"
                + "table {border-top: 1px solid #2d2d2d}\n"
                + "td {min-width: 50px; text-align:center; border-left: 1px solid #2d2d2d; border-bottom:1px solid #2d2d2d; padding: 2px 5px}"
                + "</style>");
        sb.append("<body>");
        sb.append("<table cellpadding=\"0\" cellspacing=\"0\">");
        sb.append("<tr><td colspan=\"2\">Position</td>");
        for (int j = 0; j < or.getL(); j++) {
            sb.append("<td rowspan=\"2\">").append(j + 1).append("</td>");
            if (breaks.contains(j)) {
                sb.append("<td rowspan=\"2\" style=\"background:#D70060\">").append("R").append("</td>");
            }
        }
        sb.append("</tr>");
        sb.append("<tr><td>Generator</td><td>Pi</td></tr>");
        for (int k = 0; k < or.getK(); k++) {
            sb.append("<tr>");
            sb.append("<td>G").append(k).append("</td>");
            sb.append("<td>");
            sb.append(shorten(or.getPi()[k]));
            sb.append("</td>");

            for (int j = 0; j < or.getL(); j++) {

                //mu
                sb.append("<td");
                boolean flatMu = false;
                for (int v = 0; v < or.getMu()[0][0].length; v++) {
                    if (or.getMu()[j][k][v] > 1e-20 && or.getMu()[j][k][v] != 1d) {
                        flatMu = true;
                        break;
                    }
                }
                if (flatMu) {
                    sb.append(" style=\"padding:0\">");
                    sb.append("<table cellpadding=\"0\" cellspacing=\"0\" style=\"border:0\">");
                    for (int v = 0; v < or.getn(); v++) {
                        if (v + 1 < or.getn()) {
                            sb.append("<tr><td style=\"border:none;border-bottom:1px solid #2d2d2d;");
                        } else {
                            sb.append("<tr><td style=\"border:none;");
                        }
                        switch (v) {
                            case 0:
                                sb.append(" background-color: rgba(0, 161, 203,");
                                break;
                            case 1:
                                sb.append(" background-color: rgba(97,174,36,");
                                break;
                            case 2:
                                sb.append(" background-color: rgba(229,64,40,");
                                break;
                            case 3:
                                sb.append(" background-color: rgba(241,141,5,");
                                break;
                            case 4:
                                sb.append(" background-color: rgba(97,97,97,");
                                break;
                            default:
                                break;
                        }
                        sb.append(opacity(or.getMu()[j][k][v])).append(")\">");
                        sb.append(shorten(or.getMu()[j][k][v]));
                        sb.append("</tr></td>");
                    }
                    sb.append("</table>");
                } else {
                    double max = Double.MIN_VALUE;
                    Map<Double, Integer> m = new ConcurrentHashMap<>();
                    for (int v = 0; v < or.getMu()[0][0].length; v++) {
                        max = Math.max(max, or.getMu()[j][k][v]);
                        m.put(or.getMu()[j][k][v], v);
                    }
                    if (max != Double.MIN_VALUE) {
                        switch (m.get(max)) {
                            case 0:
                                sb.append(" style=\"background-color:").append(A).append("\"");
                                break;
                            case 1:
                                sb.append(" style=\"background-color:").append(C).append("\"");
                                break;
                            case 2:
                                sb.append(" style=\"background-color:").append(G).append("\"");
                                break;
                            case 3:
                                sb.append(" style=\"background-color:").append(T).append("\"");
                                break;
                            case 4:
                                sb.append(" style=\"background-color:").append(gap).append("\"");
                                break;
                            default:
                                break;
                        }
                    }
                    sb.append(">");
                    if (max != Double.MIN_VALUE) {
                        sb.append(reverse(m.get(max)));
                    } else {
                        sb.append(0);
                    }
                }
                sb.append("</td>");
                if (breaks.contains(j)) {
                    sb.append("<td style=\"padding:0\">");
                    sb.append("<table cellpadding=\"0\" cellspacing=\"0\" style=\"border:0; width: 100%\">");
                    int r = 0;
                    for (int l = 0; l < or.getK(); l++) {
                        if (k != l) {
                            if (or.getRho()[j][k][l] > 1e-10) {
                                r++;
                            }
                        }

                    }
                    if (r != 0) {
                        for (int l = 0; l < or.getK(); l++) {
                            if (k != l) {
                                if (or.getRho()[j][k][l] > 1e-10) {
                                    r--;
                                    if (r > 0) {
                                        sb.append("<tr><td style=\"border:none;border-bottom:1px solid #2d2d2d;\">");
                                    } else {
                                        sb.append("<tr><td style=\"border:none;\">");
                                    }
                                    sb.append("G").append(l).append(":").append(shorten(or.getRho()[j][k][l]));
                                    sb.append("</td></tr>");
                                }
                            }
                        }
                    }
                    sb.append("</table>");
                    sb.append("</td>");
                }

            }

            sb.append("</tr>");
        }
        sb.append("<tr><td colspan=\"2\">Epsilon</td>");
        for (int j = 0; j < or.getL(); j++) {
            sb.append("<td>").append(shorten(or.getEps()[j])).append("</td>");
            if (breaks.contains(j)) {
                sb.append("<td style=\"background:#D70060\"></td>");
            }
        }
        sb.append("</tr>");
        sb.append("</table></body></html>");

        return sb.toString();
    }
}
