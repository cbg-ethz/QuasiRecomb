/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.informationholder;

import ch.ethz.bsse.quasirecomb.utils.BitMagic;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author toepfera
 */
public class ReadTest {

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    private Read read;

    public ReadTest() {
//        Read read1 = new Read(BitMagic.splitReadIntoBytes("ACGT"), 0, 4, BitMagic.splitReadIntoBytes("ACGT"), 10, 14);
//        for (int i = read1.getBegin(); i <= read1.getEnd(); i++) {
//            System.out.println(i + ":" + read1.getPosition(i) + " " + read1.getBase(i));
//        }
////        Read[] reads = FastaParser.parseFastaPairedEnd("/Users/XLR/Dropbox/simulationStudy/reads.fasta");
//        System.out.println("");
        //        this.read = FastaParser.parseFastq("/Users/XLR/Dropbox/QuasiAsterisk/QuasiRecomb/src/main/resources/haplotypes/singlePairedEnd.fastq")[0];
        //        System.out.println("");
        //        Globals.getINSTANCE().setALIGNMENT_BEGIN(Math.min(read.getBegin(), Globals.getINSTANCE().getALIGNMENT_BEGIN()));
        //        Globals.getINSTANCE().setALIGNMENT_END(Math.max(read.getEnd(), Globals.getINSTANCE().getALIGNMENT_END()));
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
//        Globals.getINSTANCE().setDEBUG(true);
//        Read read1 = new Read(BitMagic.splitReadIntoBytes("ACGT"), 0, 4, BitMagic.splitReadIntoBytes("GTAC"), 2, 6);
//        System.out.println("");
//        for (int i = 0; i < read1.getLength(); i++) {
//            System.out.print(read1.getBase(i));
//        }
//        System.out.println("");
//        assertTrue(read1.equals(new Read(new byte[]{1,2,3}, 0, 12,new byte[]{4,5,6,1}, 50, 100)));
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
