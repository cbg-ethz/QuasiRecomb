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
 - Input data can be in BAM or SAM format

- - -

#### PREREQUISITES TO RUN:
 - JDK 7 (http://jdk7.java.net/)

## HOW-TO:
If you are new to QuasiRecomb, please read the **[Beginners' guide to viral population inference](https://github.com/armintoepfer/QuasiRecomb/blob/master/HOWTO.md)**

## RUN:
### Local / Global reconstruction
 `java -jar QuasiRecomb.jar -i alignment.bam`
 Reads need to be properly aligned.

### Conservative reconstruction
 `-conservative`
  In this case, only major haplotypes will be reconstructed.

### Disregard deletions
 `-noGaps`
  If deletions are not of interest, not expected, or only due to technical noise, all deletions will be ignored.

### Use fixed number of generators
 `-K 2`

### Reconstruct specific region with respect to reference genome numbering
 `-r 790-2292`

### Robust model selection
 `-bootstrap`

### Incorporate PHRED quality scores (slower)
 `-quality`

### Disable recombination process
 `-noRecomb`

### Output
 The reconstructed DNA haplotype distribution quasispecies.fasta will be saved in the working directory.
 An amino acid translation of the quasispecies in all three reading frame is saved as support/quasispecies_protein_(0|1|2).fasta, if `-protein` is used.
 
### Plots
 Summary statistics can be produced with R:
```
R CMD BATCH support/coverage.R
R CMD BATCH support/modelselection.R
```

## Technical details
##### Memory consumption
To minimize the memory consumption and the number of full garbage collector executions, use:
`java -XX:NewRatio=9 -jar QuasiRecomb.jar`

If your dataset is very large and you run out of memory, increase the heapspace with:
`java -XX:NewRatio=9 -Xms2G -Xmx10G -jar QuasiRecomb.jar`

### Help:
 Further help can be showed by running without additional parameters:
  `java -jar QuasiRecomb.jar`

## PREREQUISITES COMPILE (only for dev):
 - Maven 3 (http://maven.apache.org/)

## INSTALL (only for dev):
    cd QuasiRecomb
    mvn -DartifactId=samtools -DgroupId=net.sf -Dversion=1.8.9 -Dpackaging=jar -Dfile=src/main/resources/jars/sam-1.89.jar -DgeneratePom=false install:install-file
    mvn clean package
    java -jar QuasiRecomb/target/QuasiRecomb.jar

# CONTACT:
    Armin Töpfer
    armin.toepfer (at) gmail.com
    http://www.bsse.ethz.ch/cbg/people/toepfera

# LICENSE:
 GNU GPLv3 http://www.gnu.org/licenses/gpl-3.0
