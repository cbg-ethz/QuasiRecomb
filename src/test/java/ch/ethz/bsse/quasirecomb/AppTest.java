package ch.ethz.bsse.quasirecomb;

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.utils.FastaParser;
import ch.ethz.bsse.quasirecomb.utils.Utils;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for simple App.
 */
public class AppTest
        extends TestCase {

    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest(String testName) {
        super(testName);
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite() {
        return new TestSuite(AppTest.class);
    }

    /**
     * Rigourous Test :-)
     */
    public void testApp() {
//        Read[] reads = FastaParser.parseFastaPairedEnd("C:/Users/XLR/Dropbox/simulationStudy/a1.fasta");
//        String genome = FastaParser.parseFarFile("C:/Users/XLR/Dropbox/simulationStudy/haplotypes/haploytpes_1_1.fasta")[0];
//        for (Read r : reads) {
//            if (genome.contains(Utils.reverse(r.getSequence()))) {
//                System.out.println("a");
//            } else {
//                System.out.println("'''''''");
//            }
//        }
        //        double eArray[][][] = new double[5000][500][500];
        //        long minus = 0;
        //        long plus = 0;
        //        for (int x = 0; x < 10; x++) {
        //
        //            long time = System.currentTimeMillis();
        //            for (int j = 4999; j >= 0; j--) {
        //                for (int k = 499; k >= 0; k--) {
        //                    for (int i = 499; i >= 0; i--) {
        //                        eArray[j][k][i] = 1d / 5;
        //                    }
        //                }
        //            }
        //            minus += (System.currentTimeMillis() - time);
        //            time = System.currentTimeMillis();
        //            for (int j = 0; j < 5000; j++) {
        //                for (int k = 0; k < 500; k++) {
        //                    for (int i = 0; i < 500; i++) {
        //                        eArray[j][k][i] = 1d / 5;
        //                    }
        //                }
        //            }
        //            plus += (System.currentTimeMillis() - time);
        //        }
        //        System.out.println(minus/10d);
        //        System.out.println(plus/10d);
    }
}
