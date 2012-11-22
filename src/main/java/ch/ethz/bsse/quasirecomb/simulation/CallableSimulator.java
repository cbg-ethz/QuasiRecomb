package ch.ethz.bsse.quasirecomb.simulation;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.utils.BitMagic;
import ch.ethz.bsse.quasirecomb.utils.Frequency;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import org.javatuples.Pair;

/**
 *
 * @author XLR
 */
public class CallableSimulator implements Callable<Pair<Read, Read>> {

    private int length;
    private double epsilon;
    private int n;
    private String[] haplotypes;
    private int hap;
    private int start;
    private int start2;

    public CallableSimulator(int length, double epsilon, int n, String[] haplotypes, int hap, int start, int start2) {
        this.length = length;
        this.epsilon = epsilon;
        this.n = n;
        this.haplotypes = haplotypes;
        this.hap = hap;
        this.start = start;
        this.start2 = start2;
    }

    @Override
    public Pair<Read, Read> call() throws Exception {
        char[] readArray = new char[length];
        for (int j = 0; j < length; j++) {
            //error
            if (epsilon > 0d) {
                Map<Character, Double> baseMap = new ConcurrentHashMap<>();
                for (int v = 0; v < n; v++) {
                    char x = reverse(v);
                    if (haplotypes[hap].charAt(j + start) == x) {
                        baseMap.put(x, 1.0 - (n - 1.0) * epsilon);
                    } else {
                        baseMap.put(x, epsilon);
                    }
                }
                Frequency<Character> errorF = new Frequency<>(baseMap);
                readArray[j] = errorF.roll();
            } else {
                readArray[j] = haplotypes[hap].charAt(j + start);
            }
        }
        StringBuilder sb = new StringBuilder(length);
        for (int j = 0; j < length; j++) {
            sb.append(readArray[j]);
        }
        Read r1 = new Read(BitMagic.splitReadIntoBytes(sb.toString()), start, start + length);

        readArray = new char[length];
        for (int j = 0; j < length; j++) {
            final int pos = j + start2;
            //error
            if (epsilon > 0d) {
                Map<Character, Double> baseMap = new ConcurrentHashMap<>();
                for (int v = 0; v < n; v++) {
                    char x = reverse(v);
                    if (haplotypes[hap].charAt(pos) == x) {
                        baseMap.put(x, 1.0 - (n - 1.0) * epsilon);
                    } else {
                        baseMap.put(x, epsilon);
                    }
                }
                Frequency<Character> errorF = new Frequency<>(baseMap);
                readArray[j] = errorF.roll();
            } else {
                readArray[j] = haplotypes[hap].charAt(pos);
            }
        }
        StringBuilder sb2 = new StringBuilder(length);
        for (int j = 0; j < length; j++) {
            sb2.append(readArray[j]);
        }
        Read r2 = new Read(BitMagic.splitReadIntoBytes(sb2.toString()), start2, start2 + length);
        return Pair.with(r1, r2);
    }

    private static char reverse(int v) {
        switch ((short) v) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            case 4:
                return '-';
            default:
                throw new IllegalStateException("cannot reverse " + v);
        }
    }
}
