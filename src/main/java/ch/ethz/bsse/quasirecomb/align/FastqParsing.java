/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.align;

//import java.io.File;
//import java.io.IOException;
//import java.util.Iterator;
//import java.util.logging.Level;
//import java.util.logging.Logger;
//import org.biojava3.sequencing.io.fastq.Fastq;
//import org.biojava3.sequencing.io.fastq.FastqReader;
//import org.biojava3.sequencing.io.fastq.SangerFastqReader;

/**
 *
 * @author toepfera
 */
public class FastqParsing {

    public static void main(String[] args) {

//        FastqReader reader = new SangerFastqReader();
//        try {
//            Iterable<Fastq> read = reader.read(new File("/Users/XLR/Dropbox/SRR069887.fastq"));
//            Iterator<Fastq> iterator = read.iterator();
//            while (iterator.hasNext()) {
//                Fastq next = iterator.next();
//                for (int i = 0; i < next.getSequence().length(); i++) {
//                    System.out.println(i
//                            +"\t"+next.getVariant().errorProbability(next.getQuality().charAt(i))
//                            +"\t"+next.getVariant().qualityScore(next.getQuality().charAt(i)));
//                }
//                System.out.println("");
//            }
//        } catch (IOException ex) {
//            Logger.getLogger(FastqParsing.class.getName()).log(Level.SEVERE, null, ex);
//        }



//        Map<String, Double> hapMap = new HashMap<>();
//        try {
//            FileInputStream fstream = new FileInputStream("/Users/XLR/Dropbox/SRR069887.fastq");
//
//            String head = null;
//
//            try (DataInputStream in = new DataInputStream(fstream)) {
//                BufferedReader br = new BufferedReader(new InputStreamReader(in));
//                String strLine;
//                int i = 0;
//                Fastq fastq = null;
//                List<Fastq> fastqList = new ArrayList<>();
//                while ((strLine = br.readLine()) != null) {
//                    if (strLine.startsWith("@")) {
//                        if (fastq != null) {
//                            fastqList.add(fastq);
//                        }
//                        fastq = new Fastq();
//                        fastq.setIdentifier(strLine);
//                        i = 0;
//                    } else {
//                        if (i == 1) {
//                            fastq.setRead(Utils.splitReadIntoByteArray(strLine));
//                        } else if (i == 2) {
//                            fastq.setDescription(strLine);
//                        } else if (i == 3) {
//                            fastq.setQuality(Utils.splitQualityIntoByteArray(strLine));
//                        }
//                    }
//                    i++;
//                }
//                fastqList.add(fastq);
//            }
//        } catch (IOException | NumberFormatException e) {
//            System.err.println("Error fastq: " + e.getMessage());
//        }
    }
}
