# QuasiRecomb
RNA viruses are present in a single host as a population of different
but related strains. This population, shaped by the combination of
genetic change and selection, is called quasispecies. Genetic change
is due to both point mutations and recombination events. We present a
jumping hidden Markov model that describes the generation of the viral
quasispecies and a method to infer its parameters by analysing next
generation sequencing data. We offer an implementation of the
EM algorithm to find maximum a posteriori estimates of the model
parameters and a method to estimate the distribution of viral strains
in the quasispecies. The model is validated on simulated data, showing
the advantage of explicitly taking the recombination process into
account, and validated on experimental HIV samples.

### CONTENT:
This java command line application is a toolbox, combining all necessary
steps to infer a viral quasispecies from Next Generation Sequencing (NGS) data.

### CITATION:
If you use QuasiRecomb, please cite the paper <i>Töpfer et al.</i> in <a href="http://online.liebertpub.com/doi/abs/10.1089/cmb.2012.0232">Journal of Computational Biology</a>

### DOWNLOAD:
Please get the latest binary at https://sourceforge.net/projects/quasirecomb/

### FEATURES:
 - First algorithm that models the recombination process
 - Allows position-wise mutation events
 - Infers a parametric probability distribution from the underlying viral population
 - Error correction by estimating position-wise sequencing error-rates
 - Local, gene- and genome-wide reconstruction
 - Reports SNV (single nucleotide variant) posteriors
 - Incorporates paired-end information
 - Uses PHRED scores to weight each base of each read
 - Input may contain paired-end and single reads
 - Supports reads of all current sequencing technologies (454/Roche, Illumina and PacBio)
 - Suitable for amplicon and shotgun sequencing projects
 - Reports reconstructed haplotypes and their relative frequencies
 - Reports translated proteins in all three reading frames with their relative frequencies
 - Input data can be in BAM, SAM or FASTA format

- - -

#### PREREQUISITES TO RUN:
 - JDK 7 (http://jdk7.java.net/)

## RUN:
### Local reconstruction
 `java -jar QuasiRecomb.jar -i alignedReads.fasta`
 Reads need to be aligned, therefore it is only useful for local reconstruction.

### Global reconstruction
 `java -jar QuasiRecomb.jar -i alignment.bam -global`
  In this case, all insertions will be omitted, but deletions are preserved.

### Use fixed number of generators
 `java -jar QuasiRecomb.jar -i alignedReads.fasta -K 2`

### Reconstruct specific region with respect to reference genome numbering
 `java -jar QuasiRecomb.jar -i alignedReads.bam -global -r 790-2292`

### Output
 The reconstructed DNA haplotype distribution quasispecies.fasta will be saved in the working directory.
 An amino acid translation of the quasispecies in all three reading frame is saved as quasispecies_p(0|1|2).fasta.

### Help:
 Further help can be showed by running without additional parameters:
  `java -jar QuasiRecomb.jar`

## PREREQUISITES COMPILE (only for dev):
 - Maven 3 (http://maven.apache.org/)

## INSTALL (only for dev):
    cd QuasiRecomb
    mvn -DartifactId=samtools -DgroupId=net.sf -Dversion=1.8.4 -Dpackaging=jar -Dfile=src/main/resources/jars/sam-1.84.jar -DgeneratePom=false install:install-file
    mvn clean package
    java -jar QuasiRecomb/target/QuasiRecomb.jar

# CONTACT:
    Armin Töpfer
    armin.toepfer (at) gmail.com
    http://www.bsse.ethz.ch/cbg/people/toepfera

# LICENSE:
 GNU GPLv3 http://www.gnu.org/licenses/gpl-3.0
