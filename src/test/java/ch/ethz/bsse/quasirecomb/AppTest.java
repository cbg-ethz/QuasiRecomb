package ch.ethz.bsse.quasirecomb;

import cc.mallet.types.Dirichlet;
import ch.ethz.bsse.quasirecomb.utils.BitMagic;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apfloat.Apfloat;

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

    private static int getBit(byte[] data, int pos) {
        int posByte = pos / 8;
        int posBit = pos % 8;
        byte valByte = data[posByte];
        int valInt = valByte >> (8 - (posBit + 1)) & 0x0001;
        return valInt;
    }

    private static void setBit(byte[] data, int pos, int val) {
        int posByte = pos / 8;
        int posBit = pos % 8;
        byte oldByte = data[posByte];
        oldByte = (byte) (((0xFF7F >> posBit) & oldByte) & 0x00FF);
        byte newByte = (byte) ((val << (8 - (posBit + 1))) | oldByte);
        data[posByte] = newByte;
    }

    private static String byteToBits(byte b) {
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < 8; i++) {
            buf.append((int) (b >> (8 - (i + 1)) & 0x0001));
        }
        return buf.toString();
    }

    public static byte[] charToBytes(String s) {
        int l = s.length() * 3;
        int byteCount = l / 8 + (l % 8 != 0 ? 1 : 0);
        byte[] packed = new byte[byteCount];
        int pos = 0;
        for (char c : s.toCharArray()) {
            switch ((short) c) {
                case 65:
                    break;
                case 67:
                    setBit(packed, pos * 3 + 2, 1);
                    break;
                case 71:
                    setBit(packed, pos * 3 + 1, 1);
                    break;
                case 84:
                    setBit(packed, pos * 3 + 1, 1);
                    setBit(packed, pos * 3 + 2, 1);
                    break;
                case 45:
                    setBit(packed, pos * 3 + 0, 1);
                    break;
                default:
                    break;
            }
            pos++;
//            for (int i = 0; i < packed.length; i++) {
//                System.out.print(byteToBits(packed[i]) + " ");
//            }
//            System.out.println("");
        }
        return packed;
    }

    public static char getPosition(byte[] packed, int i) {
        if (getBit(packed, i * 3 + 0) == 1) {
            return '-';
        } else if (getBit(packed, i * 3 + 1) == 1) {
            if (getBit(packed, i * 3 + 2) == 1) {
                return 'T';
            } else {
                return 'G';
            }
        } else if (getBit(packed, i * 3 + 2) == 1) {
            return 'C';
        } else {
            return 'A';
        }
    }

    /**
     * Rigourous Test :-)
     */
    public void testApp() throws MathException {
        double a = 0.9819770210743201;
        double b = Math.log(a);
        double c = b * 1;
        System.out.println(c);
//        for (int i = 0; i < 10; i++) {
//        String s = "ACG";
//        byte[] packed = BitMagic.splitReadIntoBytes(s);
//        assertEquals(s.length(), BitMagic.getLength(packed));
        //            long time = System.currentTimeMillis();
        //            byte[] packed = charToBytes(s);
        //
        //            String x = "";
        //            for (int j = 0; j < s.length(); j++) {
        //                x += getPosition(packed, j);
        ////            System.out.println(getBit(packed, i * 3) + "" + getBit(packed, i * 3 + 1) + "" + getBit(packed, i * 3 + 2));
        //            }
        //            System.out.println(System.currentTimeMillis() - time);
        //            assertEquals(s, x);
        //        }









    }
}
