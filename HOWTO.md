#Beginners' guide to viral population inference

## Input data
### 454/Roche
454 data comes in different formats, depending on your sequencing service provider.
##### SFF
Convert your SFF to fastq with `sff2fastq input.sff -o input.fastq`  
<b>sff2fastq</b> can be installed with:
```
git clone git://github.com/indraniel/sff2fastq.git;
cd sff2fastq;
make;
```

##### FNA+QUAL
Sometimes the sequences and quality scores are in different files, as input.fna and input.qual. These can be merged into one single fastq file.  
For this, a Python script is available:
```
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
import sys
name = sys.argv[1]
handle = open(name+".fastq", "w") #w=write
records = PairedFastaQualIterator(open(name+".fna"), open(name+".qual"))
count = SeqIO.write(records, handle, "fastq")
handle.close()
print "Converted %i records" % count
```
Save this code as <b>454tofastq.py</b> and execute it with `python 454tofastq.py input` and a <i>input.fastq</i> will be created.

##### FASTQ or FASTA
These formats can be used directly.

### Illumina
Illumina data is provided in the FASTQ format. If you used a MiSeq sequencing platform,  
you might get paired-end data with two files `input_R1_001.fastq` and `input_R2_001.fastq`.

### PacBio
The most reliable PacBio data are CCS (circular consesus sequences). These files are often called `input.ccs.fastq`. It is recommended to use the FASTQ.

##Alignment
The key to a good results is a high quality alignment of the input data. For this case InDelFixer will be used.  
InDelFixer can be downloaded with:
```
wget http://sourceforge.net/projects/indelfixer/files/latest/download
```
An alignment can be computed with:
```
java -jar InDelFixer.jar -i input.fastq -g referenceGenomes.fasta -illumina
```
or with paired-end data:
```
java -jar InDelFixer.jar -i input_R1_001.fastq -ir input_R2_001.fastq -g referenceGenomes.fasta -illumina
```
Depending on the sequencing platform, please use different options `-illumina`, `-454` or `-pacbio`.

For further information on InDelFixer see https://github.com/armintoepfer/InDelFixer

###Quality control
ALWAYS look at the alignment with your own eyes, to check the quality, for example with [Tablet](http://bioinf.scri.ac.uk/tablet/)

##Viral population inference
##### Coverage
The coverage is critical for inference and a MINIMAL coverage of 1,000x is needed to distinguish between sequencing errors and SNPs.  
To find reliable low-frequency variants, a coverage of >10,000x is needed.

QuasiRecomb is capable of providing regions with a minimum coverage of 1,000x and 10,000x.
```
java -jar QuasiRecomb.jar -i reads.sam -coverage

[...]
00:00:00:563 Compute coverage
00:00:00:566 To create a coverage plot, please execute: R CMD BATCH support/coverage.R
00:00:00:567 A coverage >1000x is in region 2929-7047
00:00:00:567 There is no region with a sufficient coverage of >10000x
```

One of these regions should be used:
```
java -jar QuasiRecomb.jar -i reads.sam -r 2929-7047
```

##### Model-selection
Usually model-selection is done automatically in the range of 1-5 generators, but in benchmark situations or if the underlying population is too diverse, model-selection for a larger range of generators can be activated with:
```
java -jar QuasiRecomb.jar -i reads.sam -r 2929-7047 -K 1:10
```

##### Global reconstruction
When the region of interest is larger than the average read-length, or the focus is only on the dominant haplotypes, please use `-global`.

##### Reduce number of false-positive haplotypes
If the distribution of haplotypes in the `quasispecies.fasta` file is too flat, the number of false-positives can be reduced with executing the same command-line call as before, but this time with an additional `-refine`. Of course, reducing the number of false-positives is a tradeoff with introducing false-negatives.
