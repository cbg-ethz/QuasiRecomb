package ch.ethz.bsse.quasirecomb.diversity;

import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.*;
import java.util.HashMap;
import java.util.Map;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

/**
 *
 * @author toepfera
 */
public class Diversity {

    private String input;
    private String output;
    private String description;
    private Map<String, byte[]> content;
    private int N;
    private int L;
    private boolean plot;

    public Diversity(String input, String output, String description, boolean plot) {
        this.plot = plot;
        this.input = input;
        this.description = description;
        this.output = output;
        this.content = this.parseFarFile(input);
        this.output += File.separator + description + File.separator;
//        String tmp = "";
//        if (input.contains(File.separator)) {
//            tmp = input.split(File.separator)[input.split(File.separator).length - 1];
//        } else {
//            tmp = input;
//        }
//        if (tmp.contains(".")) {
//            System.out.println(Arrays.toString(tmp.split("\\.")));
//            tmp = tmp.split("\\.")[0];
//        }
//        this.output += tmp;
        new File(this.output).mkdirs();
        this.N = content.size();
        this.L = content.values().iterator().next().length;
//        normalizedShannonEntropy();
    }

//    public void normalizedShannonEntropy() {
//        double sum = 0d;
//        double[] frequencies = new double[N];
//        int i = 0;
//        for (String head : content.keySet()) {
//            try {
//                frequencies[i] = Double.parseDouble(head);
//            } catch (NumberFormatException e) {
//                frequencies[i] = Double.parseDouble(head.split("_")[1]);
//            }
//            sum += frequencies[i];
//            i++;
//        }
//        if (sum > 1.00001 || sum < 9.99999) {
//            for (int j = 0; j < N; j++) {
//                frequencies[j] /= sum;
//            }
//        }
//        double shannonIndex = 0;
//        for (int j = 0; j < N; j++) {
//            shannonIndex -= frequencies[j] * Math.log(frequencies[j]) / Math.log(N);
//        }
//        Utils.saveFile(output + description + "_normalizedShannonIndex", description + "\t" + shannonIndex + "\n");
//    }
//    public void shannonIndex() {
//        double sum = 0d;
//        double[] frequencies = new double[N];
//        int i = 0;
//        for (String head : content.keySet()) {
//            try {
//                frequencies[i] = Double.parseDouble(head);
//            } catch (NumberFormatException e) {
//                frequencies[i] = Double.parseDouble(head.split("_")[1]);
//            }
//            sum += frequencies[i];
//            i++;
//        }
//        if (sum > 1.00001 || sum < 9.99999) {
//            for (int j = 0; j < N; j++) {
//                frequencies[j] /= sum;
//            }
//        }
//        double shannonIndex = 0;
//        for (int j = 0; j < N; j++) {
//            shannonIndex -= frequencies[j] * Math.log(frequencies[j]);
//        }
//        Utils.saveFile(output + description + "_shannonIndex", description + "\t" + shannonIndex + "\n");
//    }
    public void entropy() {
        double[] positionEntropy = new double[L];
        int[] coverage = new int[L];
        double[] weightedPositionEntropy = this.calcEntropy();

//        Map<String, Double> freqMap = new HashMap<>();
//        double sum = 0d;
//        for (String head : content.keySet()) {
//            try {
//                freqMap.put(head, Double.parseDouble(head));
//            } catch (NumberFormatException e) {
//                freqMap.put(head, Double.parseDouble(head.split("_")[1]));
//            }
//            sum += freqMap.get(head);
//        }
//        if (sum > 1.00001 || sum < 9.99999) {
//            for (String head : freqMap.keySet()) {
//                freqMap.put(head, freqMap.get(head) / sum);
//            }
//        }
//
//        for (int j = 0; j < L; j++) {
//            int l = 0;
//
//            Map<Byte, Double> map = new HashMap<>();
//
//            for (String head : content.keySet()) {
//                byte c = content.get(head)[j];
//                l++;
//                if (!map.containsKey(c)) {
//                    map.put(c, 0d);
//                }
//                map.put(c, map.get(c) + 1);
//            }
//            Double result = 0.0;
//            for (byte base : map.keySet()) {
//                Double frequency = (double) map.get(base) / l;
//                result -= frequency * Math.log(frequency) / Math.log(5);
//            }
//            positionEntropy[j] = result;
//            coverage[j] = l;
//        }
        StringBuilder wsb = new StringBuilder();
        double weightedOverallEntropy = 0d;
        for (double p : weightedPositionEntropy) {
            wsb.append(p).append("\n");
            weightedOverallEntropy += p;
        }
        weightedOverallEntropy /= L;
        double weightedSdev = 0d;
        for (double p : weightedPositionEntropy) {
            weightedSdev += Math.pow(p - weightedOverallEntropy, 2);
        }
        weightedSdev /= L;
//        double overallEntropy = 0d;
//        for (double p : positionEntropy) {
//            sb.append(p).append("\n");
//            overallEntropy += p;
//        }
//        overallEntropy /= L;
//        StringBuilder cb = new StringBuilder();
//        for (int i : coverage) {
//            cb.append(i).append("\n");
//        }
        String tmp = output + "tmp" + File.separator;
        new File(tmp).mkdirs();
//        Utils.saveFile(tmp + "Data", sb.toString());
        Utils.saveFile(output + description + "_siteEntropy", wsb.toString());
        if (plot) {
            Utils.saveFile(tmp + "DataWeighted", wsb.toString());
//        Utils.saveFile(tmp + "Coverage", cb.toString());
//        sb.append("png(\"" + output + description + "_siteEntropy.png\")\n");
//        sb.append("data <- read.delim(\"" + tmp + "Data\", header=F)\n");
//        sb.append("plot(seq(1,length(data$V1)),data$V1,type='l',ylab=\"Entropy\",xlab=\"site\",main=\"" + description + " Entropy\")\n");
//        sb.append("dev.off()");
//        Utils.saveFile(tmp + "_histogram.R", sb.toString());
//        sb.setLength(0);
//        execute("R CMD BATCH " + tmp + "_histogram.R");

            StringBuilder sb = new StringBuilder();
            sb.append("png(\"" + output + description + "_siteWeightedEntropy.png\")\n");
            sb.append("data <- read.delim(\"" + tmp + "DataWeighted\", header=F)\n");
            sb.append("plot(seq(1,length(data$V1)),data$V1,type='l',ylab=\"Entropy\",xlab=\"site\",main=\"" + description + " Weighted Entropy\")\n");
            sb.append("dev.off()");
            Utils.saveFile(tmp + "_weightedhistogram.R", sb.toString());
            execute("R CMD BATCH " + tmp + "_weightedhistogram.R");
        }
//        sb.setLength(0);
//        sb.append("png(\"" + output + description + "_coverage.png\")\n");
//        sb.append("data <- read.delim(\"" + tmp + "Coverage\", header=F)\n");
//        sb.append("plot(seq(1,length(data$V1)),data$V1,type='l',ylab=\"Coverage\",xlab=\"site\",main=\"" + description + " Coverage\")\n");
//        sb.append("dev.off()");
//        Utils.saveFile(tmp + "_coverage.R", sb.toString());
//        execute("R CMD BATCH " + tmp + "_coverage.R");
        Utils.saveFile(output + description + "_weightedOverallEntropy", description + "\t" + weightedOverallEntropy + "\t" + weightedSdev + "\n");
//        Utils.saveFile(output + description + "_overallEntropy", description + "\t" + overallEntropy + "\n");


//        MannWhitneyUTest mwut = new MannWhitneyUTest();
//        mwut.


    }

    public double[] calcEntropy() {
        double[] siteEntropy = new double[L];

        Map<String, Double> freqMap = new HashMap<>();
        double sum = 0d;
        for (String head : content.keySet()) {
            try {
                freqMap.put(head, Double.parseDouble(head));
            } catch (NumberFormatException e) {
                freqMap.put(head, Double.parseDouble(head.split("_")[1]));
            }
            sum += freqMap.get(head);
        }
        if (sum > 1.00001 || sum < 9.99999) {
            for (String head : freqMap.keySet()) {
                freqMap.put(head, freqMap.get(head) / sum);
            }
        }

        for (int j = 0; j < L; j++) {
            int l = 0;
            double lsum = 0d;
            Map<Byte, Double> weightedMap = new HashMap<>();

            for (String head : content.keySet()) {
                byte c = content.get(head)[j];
                l++;
                if (!weightedMap.containsKey(c)) {
                    weightedMap.put(c, 0d);
                }
                weightedMap.put(c, weightedMap.get(c) + freqMap.get(head));
                lsum += freqMap.get(head);
            }


            // calculate the entropy
            Double weightedResult = 0.0;
            for (byte base : weightedMap.keySet()) {
                Double frequency = (double) weightedMap.get(base) / lsum;
                weightedResult -= frequency * Math.log(frequency) / Math.log(5);
            }
            siteEntropy[j] = weightedResult;
        }
        return siteEntropy;
    }

    private Map<String, byte[]> parseFarFile(String location) {
        Map<String, byte[]> map = new HashMap<>();
        try {
            FileInputStream fstream = new FileInputStream(location);
            StringBuilder sb;
            String head = null;
            try (DataInputStream in = new DataInputStream(fstream)) {
                BufferedReader br = new BufferedReader(new InputStreamReader(in));
                String strLine;
                sb = new StringBuilder();
                while ((strLine = br.readLine()) != null) {
                    if (strLine.startsWith(">")) {
                        if (sb.length() > 0) {
                            map.put(head, Utils.splitReadIntoByteArray(sb.toString()));
                        }
//                        String xy[] = strLine.split("###")[1].split("-");
//                        Pair p = Pair.with(Integer.parseInt(xy[0]), Integer.parseInt(xy[1]));
                        sb.setLength(0);
                        head = strLine;
                    } else {
                        sb.append(strLine);
                    }
                }
            }
            map.put(head, Utils.splitReadIntoByteArray(sb.toString()));
        } catch (Exception e) {
            System.err.println("Error Far: " + e.getMessage());
        }
        return map;
    }

    private void execute(String cmd) {
        try {
            String line;
            Process p = Runtime.getRuntime().exec(cmd);
            BufferedReader bre;
            try (BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
                bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                while ((line = bri.readLine()) != null) {
                    System.out.println(line);
                }
            }
            while ((line = bre.readLine()) != null) {
                System.out.println(line);
            }
            bre.close();
            p.waitFor();
//            System.out.println("Plotting done.");
        } catch (IOException | InterruptedException err) {
            System.err.println(err);
        }
    }
}
