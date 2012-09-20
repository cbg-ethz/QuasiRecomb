/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.informationholder;

import ch.ethz.bsse.quasirecomb.utils.FastaParser;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author toepfera
 */
public class ReadTest {

    private Read read;

    public ReadTest() {
        
//        this.read = FastaParser.parseFastq("/Users/XLR/Dropbox/QuasiAsterisk/QuasiRecomb/src/main/resources/haplotypes/singlePairedEnd.fastq")[0];
//        System.out.println("");
//        Globals.getINSTANCE().setALIGNMENT_BEGIN(Math.min(read.getBegin(), Globals.getINSTANCE().getALIGNMENT_BEGIN()));
//        Globals.getINSTANCE().setALIGNMENT_END(Math.max(read.getEnd(), Globals.getINSTANCE().getALIGNMENT_END()));
    }

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of setCount method, of class Read.
     */
    @Test
    public void testInsertSize() {
//        Read read1 = new Read(new byte[]{1,2,3}, 0, 12);
//        assertTrue(read1.equals(new Read(new byte[]{1,2,3}, 0, 12)));
//        assertFalse(read1.equals(new Read(new byte[]{1,1,3}, 0, 12)));
//        assertFalse(read1.equals(new Read(new byte[]{1,2,3}, 1, 12)));
//        assertFalse(read1.equals(new Read(new byte[]{1,2,3}, 0, 2)));
        //        int hits = 0;
        //        for (int i = this.read.getBegin(); i < this.read.getEnd(); i++) {
        //            boolean hit = this.read.isHit(i-Globals.getINSTANCE().getALIGNMENT_BEGIN());
        //            if (hit) {hits++;}
        //            System.out.println(i + "\t" + hit+"\t"+(hit?this.read.getBase(i-Globals.getINSTANCE().getALIGNMENT_BEGIN()):""));
        //        }
        //        assertEquals(300, hits);
        //        assertEquals(500, this.read.getLength());
        //        assertEquals(200, this.read.getInsertSize());
    }
}
