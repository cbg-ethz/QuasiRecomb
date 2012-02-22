package ch.ethz.bsse.quasirecomb;

import ch.ethz.bsse.quasirecomb.model.Globals;
import ch.ethz.bsse.quasirecomb.simulation.Sampling;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

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
//        Sampling.fromHaplotypesGlobal("/Users/toepfera/Dropbox/QuasiRecomb/src/main/resources/haplotypes/dataset_1.fasta", 
//                5000, 300, Globals.SAMPLING_EPSILON, new double[]{.8,.1,.05,.05}, 4, "/Users/toepfera/Desktop/DATASET1");
    }
}
