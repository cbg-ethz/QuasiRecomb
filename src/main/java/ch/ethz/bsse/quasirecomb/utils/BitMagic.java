/**
 * Copyright (c) 2011-2012 Armin Töpfer
 *
 * This file is part of QuasiRecomb.
 *
 * QuasiRecomb is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * QuasiRecomb is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * QuasiRecomb. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.quasirecomb.utils;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class BitMagic {

    private static int getBit(byte[] data, int pos) {
        return data[pos / 8] >> (8 - ((pos % 8) + 1)) & 0x0001;
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

    public static byte[] pack(byte[] s) {
        int l = s.length * 3;
        int hangover = (l + 3) % 8;
        int byteCount = (3 + l) / 8 + (hangover != 0 ? 1 : 0);

        byte[] packed = new byte[byteCount];
        
        switch (8 - hangover) {
            case 8:
                break;
            case 1:
                setBit(packed, 2, 1);
                break;
            case 2:
                setBit(packed, 1, 1);
                break;
            case 3:
                setBit(packed, 1, 1);
                setBit(packed, 2, 1);
                break;
            case 4:
                setBit(packed, 0, 1);
                break;
            case 5:
                setBit(packed, 0, 1);
                setBit(packed, 2, 1);
                break;
            case 6:
                setBit(packed, 0, 1);
                setBit(packed, 1, 1);
                break;
            case 7:
                setBit(packed, 0, 1);
                setBit(packed, 1, 1);
                setBit(packed, 2, 1);
                break;
            default:
                throw new IllegalStateException("BitMagic length problem: " + (8 - hangover));
        }
        int pos = 1;
        for (byte c : s) {
            switch ((short) c) {
                case 0:
                    break;
                case 1:
                    setBit(packed, pos * 3 + 2, 1);
                    break;
                case 2:
                    setBit(packed, pos * 3 + 1, 1);
                    break;
                case 3:
                    setBit(packed, pos * 3 + 1, 1);
                    setBit(packed, pos * 3 + 2, 1);
                    break;
                case 4:
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
    public static byte[] splitReadIntoBytes(String s) {
        int l = s.length() * 3;
        int hangover = (l + 3) % 8;
        int byteCount = (3 + l) / 8 + (hangover != 0 ? 1 : 0);

        byte[] packed = new byte[byteCount];
        
        switch (8 - hangover) {
            case 8:
                break;
            case 1:
                setBit(packed, 2, 1);
                break;
            case 2:
                setBit(packed, 1, 1);
                break;
            case 3:
                setBit(packed, 1, 1);
                setBit(packed, 2, 1);
                break;
            case 4:
                setBit(packed, 0, 1);
                break;
            case 5:
                setBit(packed, 0, 1);
                setBit(packed, 2, 1);
                break;
            case 6:
                setBit(packed, 0, 1);
                setBit(packed, 1, 1);
                break;
            case 7:
                setBit(packed, 0, 1);
                setBit(packed, 1, 1);
                setBit(packed, 2, 1);
                break;
            default:
                throw new IllegalStateException("BitMagic length problem: " + (8 - hangover));
        }
        int pos = 1;
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

    public static int getLength(byte[] packed) {
        StringBuilder buf = new StringBuilder();
        buf.append((int) (packed[0] >> (8 - (1)) & 0x0001));
        buf.append((int) (packed[0] >> (8 - (2)) & 0x0001));
        buf.append((int) (packed[0] >> (8 - (3)) & 0x0001));
        int hangover = 0;
        switch (buf.toString()) {
            case "000":
                hangover = 0;
                break;
            case "001":
                hangover = 1;
                break;
            case "010":
                hangover = 2;
                break;
            case "011":
                hangover = 3;
                break;
            case "100":
                hangover = 4;
                break;
            case "101":
                hangover = 5;
                break;
            case "110":
                hangover = 6;
                break;
            case "111":
                hangover = 7;
                break;
        }
        return (packed.length * 8 - 3 - hangover) / 3;
    }
    
    public static String toString(byte[] packed) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < getLength(packed); i++) {
            sb.append(getPosition(packed, i));
        }
        return sb.toString();
    }

    public static byte getPosition(byte[] packed, int i) {
        i += 1;
        if (getBit(packed, i * 3 + 0) == 1) {
            return 4;
        } else if (getBit(packed, i * 3 + 1) == 1) {
            if (getBit(packed, i * 3 + 2) == 1) {
                return 3;
            } else {
                return 2;
            }
        } else if (getBit(packed, i * 3 + 2) == 1) {
            return 1;
        } else {
            return 0;
        }
    }
}
