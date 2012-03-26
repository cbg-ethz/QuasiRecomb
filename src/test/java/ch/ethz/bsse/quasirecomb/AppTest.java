package ch.ethz.bsse.quasirecomb;

import ch.ethz.bsse.quasirecomb.distance.DistanceUtils;
import java.util.HashMap;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import java.util.Map;

/**
 * Unit test for simple App.
 */
public class AppTest 
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    /**
     * Rigourous Test :-)
     */
    public void testApp()
    {
        Map<String,Integer> q = new HashMap<>();
        q.put("A", 100);
        q.put("B", 100);
        q.put("C", 100);
        q.put("D", 100);
        for (int i = 1; i < 50; i++) {
            q.put(""+i, 1);
        }
        Map<String,Integer> p = new HashMap<>();
        p.put("A", 100);
        p.put("B", 100);
        p.put("C", 100);
        p.put("D", 100);
        for (int i = 1; i < 100; i++) {
            p.put(""+i, 1);
        }
        
        System.out.println(Math.sqrt(DistanceUtils.calculateKLD2(p, q)+DistanceUtils.calculateKLD2(q, p)));
        System.out.println();
//        Sampling.fromHaplotypesGlobal("/Users/toepfera/Dropbox/QuasiRecomb/src/main/resources/haplotypes/dataset_1.fasta", 
//                5000, 300, Globals.SAMPLING_EPSILON, new double[]{.8,.1,.05,.05}, 4, "/Users/toepfera/Desktop/DATASET1");
    }
}
