/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.utils;

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
public class BitMagicTest {

    public BitMagicTest() {
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
     * Test of splitReadIntoBytes method, of class BitMagic.
     */
    @Test
    public void testSplitReadIntoBytes() {
        System.out.print("splitReadIntoBytes:\t");
        int l = (int) (Math.random() * 100);
        String s = "";
        for (int i = 0; i < l; i++) {
            double r = Math.random();
            if (r < .2) {
                s += "A";
            } else if (r < .4) {
                s += "C";
            } else if (r < .6) {
                s += "G";
            } else if (r < .8) {
                s += "T";
            } else {
                s += "-";
            }
        }
        byte[] expResult = null;
        byte[] packed = BitMagic.splitReadIntoBytes(s);
        String x = "";
        for (int j = 0; j < s.length(); j++) {
            switch (BitMagic.getPosition(packed, j)) {
                case 0:
                    x += "A";
                    break;
                case 1:
                    x += "C";
                    break;
                case 2:
                    x += "G";
                    break;
                case 3:
                    x += "T";
                    break;
                case 4:
                    x += "-";
                    break;
            }
        }
        assertEquals(s, x);
        System.out.println("done");
        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
    }

    /**
     * Test of getLength method, of class BitMagic.
     */
    @Test
    public void testGetLength() {
        System.out.print("getLength:");
        String s = "ACGT-";
        byte[] packed = BitMagic.splitReadIntoBytes(s);
        int expResult = s.length();
        int result = BitMagic.getLength(packed);
        assertEquals(expResult, result);
        System.out.println("\t\tdone");
        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
    }
}