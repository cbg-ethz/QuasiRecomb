/**
 * Copyright (c) 2011-2013 Armin Töpfer
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

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.ReadTMP;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import net.sf.samtools.CigarElement;
import static net.sf.samtools.CigarOperator.D;
import static net.sf.samtools.CigarOperator.EQ;
import static net.sf.samtools.CigarOperator.H;
import static net.sf.samtools.CigarOperator.I;
import static net.sf.samtools.CigarOperator.M;
import static net.sf.samtools.CigarOperator.N;
import static net.sf.samtools.CigarOperator.P;
import static net.sf.samtools.CigarOperator.S;
import static net.sf.samtools.CigarOperator.X;
import net.sf.samtools.SAMRecord;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SFRComputing implements Callable<List<ReadTMP>> {

    final List<SAMRecord> samRecordList;

    public SFRComputing(final List<SAMRecord> samRecord) {
        this.samRecordList = samRecord;
    }

    @Override
    public List<ReadTMP> call() {
        List<ReadTMP> results = new LinkedList<>();
        for (SAMRecord s : samRecordList) {
            ReadTMP r = single(s);
            if (r != null) {
                results.add(r);
            }
        }
        return results;
    }

    private ReadTMP single(SAMRecord samRecord) {
        try {
            if (samRecord.getAlignmentBlocks().isEmpty()) {
                return null;
            }
            int refStart = samRecord.getAlignmentStart() - 1;
            int readStart = 0;
            int readLength_ = samRecord.getCigar().getReadLength();
            List<Byte> buildRead = new ArrayList<>(readLength_);
            List<Double> buildQuality = new ArrayList<>(readLength_);
            boolean hasQuality;
            if (Globals.getINSTANCE().isNO_QUALITY()) {
                hasQuality = false;
            } else {
                hasQuality = samRecord.getBaseQualities().length > 1;
            }
            List<Boolean> buildCigar = new ArrayList<>(readLength_);
            byte[] readBases_pure = samRecord.getReadBases();
            for (CigarElement c : samRecord.getCigar().getCigarElements()) {
                switch (c.getOperator()) {
                    case X:
                    case EQ:
                    case M:
                        if ((readStart + c.getLength()) > readBases_pure.length) {
                            System.out.println("\nInput alignment is corrupt.\nCIGAR is longer than actual read-length.");
                            System.exit(9);
                        }
                        for (int i = 0; i < c.getLength(); i++) {
                            byte b = readBases_pure[readStart];
                            buildRead.add(b);
                            if (hasQuality) {
                                double q = 1 - Math.pow(10, -(samRecord.getBaseQualities()[readStart]) / 10d);
                                if (q == 0) {
//                                    q = 0.79432823472;
                                    q = 0.01;
                                }
                                buildQuality.add(q);
                            } else {
                                buildQuality.add(1d);
                            }
                            buildCigar.add(true);
                            readStart++;
                        }
                        break;
                    case I:
                        readStart += c.getLength();
                        break;
                    case D:
                        for (int i = 0; i < c.getLength(); i++) {
                            buildRead.add((byte) "-".charAt(0));
                            double q;

                            if (c.getLength() % 3 == 0) {
                                q = 0.79432823472;
                            } else {
                                if (Globals.getINSTANCE().isNO_GAPS()) {
                                    q = 0.0;
                                } else {
                                    q = 0.01;
                                }
                            }
                            buildQuality.add(q);
                            buildCigar.add(false);
                        }
                        break;
                    case S:
                        readStart += c.getLength();
                        break;
                    case H:
                        break;
                    case P:
                        System.out.println("P");
                        System.exit(9);
                        break;
                    case N:
                        System.out.println("N");
                        System.exit(9);
                        break;
                    default:
                        break;
                }
            }

            double[] quality = null;
            boolean[] cigar = null;
            //---
            //cut read
            byte[] readBases = Utils.convertRead(buildRead.toArray(new Byte[buildRead.size()]));
            int readEnd = refStart + readBases.length;

            double delCountSum = 0;
            double delCountMax = 0;
            double delCount = 0;
            if (Globals.getINSTANCE().isWINDOW()) {

                int from = Globals.getINSTANCE().getWINDOW_BEGIN();
                int to = Globals.getINSTANCE().getWINDOW_END();
                int length = readBases.length;
                try {
                    if (refStart > to || readEnd < from) {
                        return null;
                    }
                    //leftover
                    if (refStart < from && readEnd <= to) {
//                        readBases = Arrays.copyOfRange(readBases, from - refStart, length);
//                        cigar = Arrays.copyOfRange(cigar, from - refStart, length);
//                        quality = Arrays.copyOfRange(quality, from - refStart, length);
                        int l = length - (from - refStart);
                        byte[] readBases_ = new byte[l];
                        cigar = new boolean[l];
                        quality = new double[l];
                        for (int i = from - refStart, j = 0; i < length; i++, j++) {
                            readBases_[j] = readBases[i];
                            cigar[j] = buildCigar.get(i);
                            quality[j] = buildQuality.get(i);

                            if (!cigar[j]) {
                                delCount++;
                                delCountSum++;
                            } else {
                                if (delCount > delCountMax) {
                                    delCountMax = delCount;
                                }
                                delCount = 0;
                            }
                        }
                        readBases = readBases_;
                        refStart = from;
                    } else if (refStart >= from && readEnd > to) {
                        //rightover
//                        readBases = Arrays.copyOfRange(readBases, 0, to - refStart);
//                        cigar = Arrays.copyOfRange(cigar, 0, to - refStart);
//                        quality = Arrays.copyOfRange(quality, 0, to - refStart);

                        int l = to - refStart;
                        byte[] readBases_ = new byte[l];
                        cigar = new boolean[l];
                        quality = new double[l];
                        for (int i = 0; i < l; i++) {
                            readBases_[i] = readBases[i];
                            cigar[i] = buildCigar.get(i);
                            quality[i] = buildQuality.get(i);

                            if (!cigar[i]) {
                                delCount++;
                                delCountSum++;
                            } else {
                                if (delCount > delCountMax) {
                                    delCountMax = delCount;
                                }
                                delCount = 0;
                            }
                        }
                        readBases = readBases_;
                    } else if (refStart >= from && readEnd <= to) {
                        //inner
                        int l = buildCigar.size();
                        cigar = new boolean[l];
                        quality = new double[l];
                        for (int i = 0; i < l; i++) {
                            cigar[i] = buildCigar.get(i);
                            quality[i] = buildQuality.get(i);

                            if (!cigar[i]) {
                                delCount++;
                                delCountSum++;
                            } else {
                                if (delCount > delCountMax) {
                                    delCountMax = delCount;
                                }
                                delCount = 0;
                            }
                        }
                    } else if (refStart < from && readEnd > to) {
                        //outer
//                        readBases = Arrays.copyOfRange(readBases, from - refStart, to - refStart);
//                        cigar = Arrays.copyOfRange(cigar, from - refStart, to - refStart);
//                        quality = Arrays.copyOfRange(quality, from - refStart, to - refStart);
                        int l = (to - refStart) - (from - refStart);
                        byte[] readBases_ = new byte[l];
                        cigar = new boolean[l];
                        quality = new double[l];
                        for (int i = from - refStart, j = 0; i < to - refStart; i++, j++) {
                            readBases_[j] = readBases[i];
                            cigar[j] = buildCigar.get(i);
                            quality[j] = buildQuality.get(i);

                            if (!cigar[j]) {
                                delCount++;
                                delCountSum++;
                            } else {
                                if (delCount > delCountMax) {
                                    delCountMax = delCount;
                                }
                                delCount = 0;
                            }
                        }
                        readBases = readBases_;
                        refStart = from;
                    } else {
                        System.err.println("");
                        System.err.println("start: " + refStart + "\t+end:" + readEnd);
                        System.err.println("w00t");
                    }
                } catch (Exception e) {
                    System.err.println("");
                    System.err.println("start: " + refStart + "\t+end:" + readEnd);
                    System.err.println("w00t");
                }
            } else {
                int l = buildCigar.size();
                cigar = new boolean[l];
                quality = new double[l];
                for (int i = 0; i < l; i++) {
                    cigar[i] = buildCigar.get(i);
                    quality[i] = buildQuality.get(i);

                    if (!cigar[i]) {
                        delCount++;
                        delCountSum++;
                    } else {
                        if (delCount > delCountMax) {
                            delCountMax = delCount;
                        }
                        delCount = 0;
                    }
                }
            }

            if (readBases.length < Globals.getINSTANCE().getREAD_MINLENGTH()) {
                return null;
            }

            if (delCount > delCountMax) {
                delCountMax = delCount;
            }
            if (Globals.getINSTANCE().getMAX_DEL() != -1) {
                if (delCountMax > Globals.getINSTANCE().getMAX_DEL()) {
                    return null;
                }
            }
            if (cigar == null || delCountSum / cigar.length > Globals.getINSTANCE().getMAX_PERC_DEL()) {
                return null;
            }

            String name = samRecord.getReadName();

            return new ReadTMP(name, quality, readBases, refStart, true, cigar);

        } catch (ArrayIndexOutOfBoundsException e) {
            System.err.println();
            System.err.println(e);
            System.err.println();
        } catch (Exception e) {
            System.err.println("WOOT:" + e);
            // Sometimes CIGAR is not correct. In that case we simply ignore it/
        }
        return null;
    }
}
