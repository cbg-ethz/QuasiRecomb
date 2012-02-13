package ch.ethz.bsse.quasirecomb.model;

import ch.ethz.bsse.quasirecomb.distance.DistanceUtils;
import ch.ethz.bsse.quasirecomb.model.hmm.ModelSelection;
import ch.ethz.bsse.quasirecomb.modelsampling.ModelSampling;
import ch.ethz.bsse.quasirecomb.simulation.Sampling;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
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

    private static Map<byte[], Integer> calc(Map<byte[], Integer> clusterReads, int Kmin, int Kmax, int N, int L, int n, byte[][] haplotypesArray) {
        Utils.saveClusteredReads(clusterReads);
        ModelSelection ms = new ModelSelection(clusterReads, Kmin, Kmax, N, L, n, haplotypesArray);
        if (Globals.SIMULATION) {
            ModelSampling s = new ModelSampling(L, n, ms.getBestK(), ms.getRho(), ms.getPi(), ms.getMu(), Globals.savePath);
            return s.getReads();
        }
        return null;
    }

    private static int countChars(byte[][] bbb) {
        Map<Byte,Boolean> map = new HashMap<>();
        for (byte[] bb : bbb) {
            for (byte b : bb) {
                map.put(b, Boolean.TRUE);
            }
        }
        return map.keySet().size();
    }
    
    private static void expDataset(String path, int Kmin, int Kmax) {
        byte[][] haplotypesArray = Utils.parseInput(path);
        Map<byte[], Integer> clusterReads = Utils.clusterReads(
                Utils.splitReadsIntoByteArrays(Utils.parseFarFile(path)));
        expSub(clusterReads, Kmin, Kmax, haplotypesArray[0].length, countChars(haplotypesArray), haplotypesArray);
    }

    private static void expSub(Map<byte[], Integer> clusterReads, int Kmin, int Kmax, int L, int n, byte[][] haplotypesArray) {
        // count reads
        int N = 0;
        for (Integer i : clusterReads.values()) {
            N += i;
        }
        if (Globals.MASK_RHO_THRESHOLD == -1) {
            Globals.MASK_RHO_THRESHOLD = 1d / (3 * N);
        }
        Map<byte[], Integer> simulatedReads = calc(clusterReads, Kmin, Kmax, N, L, n, haplotypesArray);
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

    private static void singleDataset(String path, double[] f, int Kmin, int Kmax, int N) {
        byte[][] haplotypesArray = Utils.parseInput(path);
        singleSub(path, N, haplotypesArray[0].length, f, countChars(haplotypesArray), Kmin, Kmax, haplotypesArray);
    }

    private static void singleSub(String path, int N, int L, double[] f, int n, int Kmin, int Kmax, byte[][] haplotypesArray) {
        Map<String, Integer> haplotypes = Sampling.fromHaplotypes(path, N, L, Globals.SAMPLING_EPSILON, f, n, Globals.savePath);
        List<String> reads = new LinkedList<>();
        List<String> hapL = new LinkedList<>();
        for (String read : haplotypes.keySet()) {
            hapL.add(read);
            for (int i = 0; i < haplotypes.get(read); i++) {
                reads.add(read);
            }
        }
        Utils.save(haplotypes, Globals.savePath + "reads.txt");
        Map<byte[], Integer> clusterReads = Utils.clusterReads(Utils.splitReadsIntoByteArrays(reads.toArray(new String[N])));

        Map<byte[], Integer> simulatedReads = calc(clusterReads, Kmin, Kmax, N, L, n, haplotypesArray);
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
