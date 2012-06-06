/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.modelsampling;

import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResult;
import ch.ethz.bsse.quasirecomb.informatioholder.OptimalResultsList;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

/**
 *
 * @author XLR
 */
public class OptimaModelSampling {

    private String savePath;
    private OptimalResultsList ors;

    public OptimaModelSampling(String input, String output) {
        try {
            FileInputStream fis = new FileInputStream(input);
            try (ObjectInputStream in = new ObjectInputStream(fis)) {
                ors = (OptimalResultsList) in.readObject();
            }
            this.sample();
        } catch (IOException | ClassNotFoundException ex) {
            System.err.println(ex);
        }
        this.savePath = output;
    }

    private void sample() {
//        System.out.println("Sampling " + input);

        List<Map<byte[], Integer>> reads = new ArrayList<>();
        Map<byte[], Integer> uniqueReads = new HashMap<>();
        double[] llhs = new double[ors.getOrs().size()];
        double max = 0;
        int i = 0;
        OptimalResult orMax = null;
        for (OptimalResult or : ors.getOrs()) {
            llhs[i++] = or.getLlh();
            if (max == 0 || or.getLlh() > max) {
                max = or.getLlh();
                orMax = or;
            }
        }
        double sd = new StandardDeviation().evaluate(llhs);
        double ose = sd / Math.sqrt(llhs.length);
        i = 0;
        System.out.println("sd\n"+sd);
        System.out.println("ose\n"+ose);
        uniqueReads.putAll(new ModelSampling(orMax).getReads());
        for (OptimalResult or : ors.getOrs()) {
//            System.out.print(Math.round(i * 100d /llhs.length) + "\r");
            System.out.println(or.getLlh());
            if (Math.abs(max - or.getLlh()) < ose && max != or.getLlh()) {
                ModelSampling simulation = new ModelSampling(or);
                reads.add(simulation.getReads());
                Map<byte[], Integer> reads_tmp = simulation.getReads();
                for (byte[] b : reads_tmp.keySet()) {
                    if (uniqueReads.containsKey(b)) {
                        uniqueReads.put(b, uniqueReads.get(b) + reads_tmp.get(b));
                    } else {
                        uniqueReads.remove(b);
                    }
                }
                i++;
            }
        }
        System.out.print("                  \r");
        System.out.println(i);
        int amount = 0;
        for (Integer x : uniqueReads.values()) {
            amount += x;
        }
        StringBuilder sb = new StringBuilder();
        i = 0;
        Map sortMapByValue = ModelSampling.sortMapByValue(uniqueReads);
        for (Object o : sortMapByValue.keySet()) {
            byte[] read = (byte[]) o;
            sb.append(">read").append(i++).append("_").append(((double) uniqueReads.get(read)) / amount).append("\n");
            for (int r : read) {
                sb.append(Utils.reverse(r));
            }
            sb.append("\n");
        }
        Utils.saveFile(savePath + "quasispecies_core.fasta", sb.toString());
    }
}
