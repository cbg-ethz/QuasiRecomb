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
package ch.ethz.bsse.quasirecomb.align;

//import org.biojava3.alignment.Alignments;
//import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
//import org.biojava3.alignment.SimpleGapPenalty;
//import org.biojava3.alignment.SubstitutionMatrixHelper;
//import org.biojava3.alignment.template.SequencePair;
//import org.biojava3.alignment.template.SubstitutionMatrix;
//import org.biojava3.core.sequence.DNASequence;
//import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
//import org.biojava3.core.sequence.compound.NucleotideCompound;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class AlignToReference {

//    public static void main(String[] args) {
//        String targetSeq = "GCGGGTGCGGTTGCTGGAAAGATGCATCTATAAC";
//
//        String querySeq = "ACGAGTGCGTGTTTTCCCGCCTGGTCCCCAGGCCCCCTTTCCGTCCTCAGGAA"
//                + "GACAGAGGAGGAGCCCCTCGGGCTGCAGGTGGTGGGCGTTGCGGCGGCGGCCGGTTAAGGT"
//                + "TCCCAGTGCCCGCACCCGGCCCACGGGAGCCCCGGACTGGCGGCGTCACTGTCAGTGTCTT"
//                + "CTCAGGAGGCCGCCTGTGTGACTGGATCGTTCGTGTCCCCACAGCACGTTTCTTGGAGTAC"
//                + "TCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTCCTGGACAGATACT"
//                + "TCCATAACCAGGAGGAGAACGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGAC"
//                + "GGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACATCCTGGAAGACGAG"
//                + "CGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGTGAGAGCTTCACCGTGCA"
//                + "GCGGCGAGACGCACTCGT";
//        DNASequence target = new DNASequence(targetSeq,
//                AmbiguityDNACompoundSet.getDNACompoundSet());
//        DNASequence query = new DNASequence(querySeq,
//                AmbiguityDNACompoundSet.getDNACompoundSet());
//
//        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
//
//        SimpleGapPenalty gapP = new SimpleGapPenalty();
//        gapP.setOpenPenalty((short) 5);
//        gapP.setExtensionPenalty((short) 2);
//
//        SequencePair<DNASequence, NucleotideCompound> psa =
//                Alignments.getPairwiseAlignment(query, target,
//                PairwiseSequenceAlignerType.LOCAL, gapP, matrix);
//
//        System.out.println(psa.getAlignedSequences().get(0).getLocationInAlignment());
//    }
}
