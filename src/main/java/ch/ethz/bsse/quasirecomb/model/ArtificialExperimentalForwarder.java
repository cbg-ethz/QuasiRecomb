package ch.ethz.bsse.quasirecomb.model;

import ch.ethz.bsse.quasirecomb.distance.DistanceUtils;
import ch.ethz.bsse.quasirecomb.model.hmm.ModelSelection;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.simulation.Sampling;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.javatuples.Pair;

/**
 * Forwards parameters depending if input is from an experimental dataset or
 * artifical haplotypes from which has to be sampled.
 *
 * @author Armin Töpfer (armin.toepfer@bsse.ethz.ch)
 */
public class ArtificialExperimentalForwarder {

    /**
     * Entry point. Forwards invokes of the specified workflow.
     *
     * @param L length of the reads
     * @param exp is it an experimental dataset
     * @param input path to the fasta file
     * @param Kmin minimal amount of generators
     * @param Kmax minimal amount of generators
     * @param n size of the alphabet
     * @param f distribution of the haplotypes if sampling has to be done
     * @param N amount of reads in case exp is false
     */
    public static void forward(boolean exp, String input, int Kmin, int Kmax, double[] f, int N) {
        if (exp) {
            expDataset(input, Kmin, Kmax);
        } else {
            singleDataset(input, f, Kmin, Kmax, N);
        }
    }

    private static ModelSampling calc(Map<byte[], Integer> clusterReads, int Kmin, int Kmax, int N, int L, int n, byte[][] haplotypesArray) {
        Utils.saveClusteredReads(clusterReads);
        ModelSelection ms = new ModelSelection(clusterReads, Kmin, Kmax, N, L, n, haplotypesArray);
        return new ModelSampling(L, n, ms.getBestK(), ms.getRho(), ms.getPi(), ms.getMu(), Globals.savePath);
    }

    private static int countChars(byte[][] bbb) {
        Map<Byte, Boolean> map = new HashMap<>();
        for (byte[] bb : bbb) {
            for (byte b : bb) {
                map.put(b, Boolean.TRUE);
            }
        }
        return map.keySet().size();
    }

    private static void expDataset(String path, int Kmin, int Kmax) {
        byte[][] haplotypesArray = Utils.parseInput(path);

        if (Globals.BOOTSTRAP) {
            Map<byte[], Double> hapMap = new HashMap<>();
            int hapL = haplotypesArray.length;

            for (int i = 0; i < hapL; i++) {
                hapMap.put(haplotypesArray[i], 1d / hapL);
            }

            Frequency<byte[]> f = new Frequency<>(hapMap);

            byte[][] haplotypesBootstrap = new byte[hapL][];
            for (int i = 0; i < hapL; i++) {
                haplotypesBootstrap[i] = f.roll();
            }
            haplotypesArray = haplotypesBootstrap;
        }

        Map<byte[], Integer> clusterReads = Utils.clusterReads(haplotypesArray);

        expSub(path, clusterReads, Kmin, Kmax, haplotypesArray[0].length, countChars(haplotypesArray), haplotypesArray);
    }

    private static void expSub(String path, Map<byte[], Integer> clusterReads, int Kmin, int Kmax, int L, int n, byte[][] haplotypesArray) {
        // count reads
        int N = 0;
        for (Integer i : clusterReads.values()) {
            N += i;
        }
        if (Globals.CROSSVALIDATION) {
            N -= N / 10;
            String[] haps = Utils.parseFarFile(path);

            double[][] dists = new double[20][10];
            double[] KLdists = new double[10];
//            double[] KLdistsR = new double[10];
            double[][] distsR = new double[20][10];
            int Nstart = 0;
            for (int i = 0; i < 10; i++) {
                System.out.println("+++++++");
                System.out.println("+++++" + i);
                System.out.println("+++++++");
                String[] hapEmp = new String[N / 10];
                List<String> hapTestList = new LinkedList();
                for (int j = 0; j < N; j++) {
                    if (j >= Nstart && j < Nstart + N / 10) {
                        hapEmp[j - Nstart] = haps[j];
                    } else {
                        hapTestList.add(haps[j]);
                    }
                }
                String[] hapTest = hapTestList.toArray(new String[hapTestList.size()]);

                ModelSampling ms = calc(Utils.clusterReads(Utils.splitReadsIntoByteArrays(hapTest)), Kmin, Kmax, N, L, n, haplotypesArray);
                Map<String, Integer> simulatedReads = ms.getReadsReversed();
                Map<String, String> simulatedReadsString = new HashMap<>();
                for (String s : simulatedReads.keySet()) {
                    simulatedReadsString.put(s, "" + simulatedReads.get(s));
                }
                Map<String, String> clusteredEmpirical = new HashMap<>();
                Map<String, Integer> clusteredEmpiricalTest = new HashMap<>();
                for (String s : hapEmp) {
                    if (!clusteredEmpirical.containsKey(s)) {
                        clusteredEmpirical.put(s, "0");
                        clusteredEmpiricalTest.put(s, 0);
                    }
                    clusteredEmpirical.put(s, String.valueOf(Integer.parseInt(clusteredEmpirical.get(s)) + 1));
                    clusteredEmpiricalTest.put(s, clusteredEmpiricalTest.get(s) + 1);
                }
                Pair[] ds = DistanceUtils.calculatePhi(clusteredEmpirical, simulatedReads);
                int j = 0;
                for (Pair d : ds) {
                    dists[j++][i] = (double) d.getValue0();
                }
                ds = DistanceUtils.calculatePhi(simulatedReadsString, clusteredEmpiricalTest);
                j = 0;
                for (Pair d : ds) {
                    distsR[j++][i] = (double) d.getValue0();
                }
                Nstart += N / 10;
                KLdists[i] = Math.sqrt(DistanceUtils.calculateKLD2(clusteredEmpiricalTest, simulatedReads) + DistanceUtils.calculateKLD2(simulatedReads, clusteredEmpiricalTest));
            }
            StringBuilder sb = new StringBuilder();
            StringBuilder sbR = new StringBuilder();
            sb.append("mean").append("\t").append("sd").append("\n");
            sbR.append("mean").append("\t").append("sd").append("\n");
            for (int j = 0; j < 20; j++) {
                Mean m = new Mean();
                StandardDeviation sd = new StandardDeviation();
                sb.append(m.evaluate(dists[j])).append("\t").append(sd.evaluate(dists[j])).append("\n");
                sbR.append(m.evaluate(distsR[j])).append("\t").append(sd.evaluate(distsR[j])).append("\n");
            }

            Utils.saveFile(Globals.savePath + File.separator + "KL-" + Kmin + ".txt", "" + new Mean().evaluate(KLdists) + "\t" + new StandardDeviation().evaluate(KLdists) + "\n");
            Utils.saveFile(Globals.savePath + File.separator + "crossvalidation-" + Kmin + ".txt", sb.toString());
            Utils.saveFile(Globals.savePath + File.separator + "crossvalidationR-" + Kmin + ".txt", sbR.toString());
        } else {
            Map<byte[], Integer> simulatedReads = calc(clusterReads, Kmin, Kmax, N, L, n, haplotypesArray).getReads();
            if (Globals.DISTANCE && Globals.SIMULATION) {
                Map<String, Integer> PHtrue = Utils.reverse(clusterReads);
                Map<String, String> PHoriginal = new HashMap<>();
                for (String key : PHtrue.keySet()) {
                    PHoriginal.put(key, PHtrue.get(key).toString());
                }

                Map<String, Integer> PHtest = Utils.reverse(simulatedReads);
                Pair[] dists = DistanceUtils.calculatePhi(PHoriginal, PHtest);
                StringBuilder sb = new StringBuilder();
                for (Pair d : dists) {
                    sb.append(d).append("\n");
                }
                if (Globals.rho0force) {
                    Utils.saveFile(Globals.savePath + File.separator + "distance-rho0.txt", "rho0\n" + sb.toString());
                } else {
                    Utils.saveFile(Globals.savePath + File.separator + "distance-phat.txt", "PHat\n" + sb.toString());
                }
            }
        }
    }

    private static void singleDataset(String path, double[] f, int Kmin, int Kmax, int N) {
        byte[][] haplotypesArray = Utils.parseInput(path);
        singleSub(path, N, haplotypesArray[0].length, f, countChars(haplotypesArray), Kmin, Kmax, haplotypesArray);
    }

    private static void singleSub(String path, int N, int L, double[] f, int n, int Kmin, int Kmax, byte[][] haplotypesArray) {
        Map<String, Integer> haplotypes;

        haplotypes = Sampling.fromHaplotypes(path, N, L, Globals.SAMPLING_EPSILON, f, n, Globals.savePath);

        List<String> reads = new LinkedList<>();
        List<String> hapL = new LinkedList<>();
        for (String read : haplotypes.keySet()) {
            hapL.add(read);
            for (int i = 0; i < haplotypes.get(read); i++) {
                reads.add(read);
            }
        }
        System.out.println("Unique: " + hapL.size());
//        Utils.save(haplotypes, Globals.savePath + "reads.txt");
        Map<byte[], Integer> clusterReads = Utils.clusterReads(Utils.splitReadsIntoByteArrays(reads.toArray(new String[N])));

        if (Globals.CROSSVALIDATION) {
            N -= N / 10;
            String[] haps = Sampling.fromHaplotypesCross(path, N, L, Globals.SAMPLING_EPSILON, f, n, Globals.savePath);

            double[][] dists = new double[20][10];
            double[] KLdists = new double[10];
//            double[] KLdistsR = new double[10];
            double[][] distsR = new double[20][10];
            int Nstart = 0;
            for (int i = 0; i < 10; i++) {
                System.out.println("+++++++");
                System.out.println("+++++" + i);
                System.out.println("+++++++");
                String[] hapEmp = new String[N / 10];
                List<String> hapTestList = new LinkedList();
                for (int j = 0; j < N; j++) {
                    if (j >= Nstart && j < Nstart + N / 10) {
                        hapEmp[j - Nstart] = haps[j];
                    } else {
                        hapTestList.add(haps[j]);
                    }
                }
                String[] hapTest = hapTestList.toArray(new String[hapTestList.size()]);

                ModelSampling ms = calc(Utils.clusterReads(Utils.splitReadsIntoByteArrays(hapTest)), Kmin, Kmax, N, L, n, haplotypesArray);
                Map<String, Integer> simulatedReads = ms.getReadsReversed();
                Map<String, String> simulatedReadsString = new HashMap<>();
                for (String s : simulatedReads.keySet()) {
                    simulatedReadsString.put(s, "" + simulatedReads.get(s));
                }
                Map<String, String> clusteredEmpirical = new HashMap<>();
                Map<String, Integer> clusteredEmpiricalTest = new HashMap<>();
                for (String s : hapEmp) {
                    if (!clusteredEmpirical.containsKey(s)) {
                        clusteredEmpirical.put(s, "0");
                        clusteredEmpiricalTest.put(s, 0);
                    }
                    clusteredEmpirical.put(s, String.valueOf(Integer.parseInt(clusteredEmpirical.get(s)) + 1));
                    clusteredEmpiricalTest.put(s, clusteredEmpiricalTest.get(s) + 1);
                }
                Pair[] ds = DistanceUtils.calculatePhi(clusteredEmpirical, simulatedReads);
                int j = 0;
                for (Pair d : ds) {
                    dists[j++][i] = (double) d.getValue0();
                }
                ds = DistanceUtils.calculatePhi(simulatedReadsString, clusteredEmpiricalTest);
                j = 0;
                for (Pair d : ds) {
                    distsR[j++][i] = (double) d.getValue0();
                }
                Nstart += N / 10;
                KLdists[i] = Math.sqrt(DistanceUtils.calculateKLD2(clusteredEmpiricalTest, simulatedReads) + DistanceUtils.calculateKLD2(simulatedReads, clusteredEmpiricalTest));
            }
            StringBuilder sb = new StringBuilder();
            StringBuilder sbR = new StringBuilder();
            sb.append("mean").append("\t").append("sd").append("\n");
            sbR.append("mean").append("\t").append("sd").append("\n");
            for (int j = 0; j < 20; j++) {
                Mean m = new Mean();
                StandardDeviation sd = new StandardDeviation();
                sb.append(m.evaluate(dists[j])).append("\t").append(sd.evaluate(dists[j])).append("\n");
                sbR.append(m.evaluate(distsR[j])).append("\t").append(sd.evaluate(distsR[j])).append("\n");
            }

            Utils.saveFile(Globals.savePath + File.separator + "KL-" + Kmin + ".txt", "" + new Mean().evaluate(KLdists) + "\t" + new StandardDeviation().evaluate(KLdists) + "\n");
            Utils.saveFile(Globals.savePath + File.separator + "crossvalidation-" + Kmin + ".txt", sb.toString());
            Utils.saveFile(Globals.savePath + File.separator + "crossvalidationR-" + Kmin + ".txt", sbR.toString());
        } else {
            Map<byte[], Integer> simulatedReads = calc(clusterReads, Kmin, Kmax, N, L, n, haplotypesArray).getReads();
            if (Globals.DISTANCE && Globals.SIMULATION) {
                Map<String, String> PHoriginal = new HashMap<>();
                for (String key : haplotypes.keySet()) {
                    PHoriginal.put(key, haplotypes.get(key).toString());
                }
                Map<String, Integer> PHtest = Utils.reverse(simulatedReads);
                Pair[] dists = DistanceUtils.calculatePhi(PHoriginal, PHtest);
                StringBuilder sb = new StringBuilder();
                for (Pair d : dists) {
                    sb.append(d).append("\n");
                }
                if (Globals.rho0force) {
                    Utils.saveFile(Globals.savePath + File.separator + "distance-rho0.txt", "rho0\n" + sb.toString());
                } else {
                    Utils.saveFile(Globals.savePath + File.separator + "distance-phat.txt", "PHat\n" + sb.toString());
                }
            }
        }
    }
}
