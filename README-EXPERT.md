#Readme for experts
This is the perfect place if you want to get the full power out of QuasiRecomb. By now you should have used QuasiRecomb to analyse some of the simulation datasets.

This document will give deep insights how QuasiRecomb works and can be tuned. By default, all parameters are set to empirical values tested on simulated and clinical data.

The structure how parameters with input is describes, will be `-parameter TYPE [default]`, for example, `-e 0.1 [0.008]`.

The execution of `java -jar QuasiRecomb.jar` will be abbreviated with `qr`. You can use it as a function in your .bashrc or .profile 

````
qr () { 
    java -jar ~/PATH/TO/QuasiRecomb.jar $*
}
````

#Input
QuasiRecomb accepts **FASTA** and **BAM** files. All insertions in the BAM will be ignored! If you want to study small insertions, please re-align your reads again the consensus sequence of the MSA.

#Standard parameters
A typical execution for a local reconstruction looks like:

````
qr -i ../simulationStudy/dataset_1.far 

00:00:00:678 Parsing    100%
00:00:00:679 Unique reads  191
00:00:02:581 Model selection [K 1]:     100%    [BIC: -23306]
00:00:07:063 Model selection [K 2]:     100%    [BIC: -13331]
00:00:14:452 Model selection [K 3]:     100%    [BIC: -16511]
00:01:00:202 Model training  [K 2]:     100%    [BIC: -13317]
00:01:01:324 Model training [K 2]:      100%
Quasispecies saved: quasispecies.fasta            
````
In the pre-processing step, the input file is parsed and the unique reads get hashed.

###Model selection
In the next step model selection is performed. The default range is `-K 1:5`, where *K* is the number of generators and `1:5` is the interval for model selection. Model selection can be skipped by fixing the number of generators, e.g., `-K 2`.  
For each *K* in the interval, 5 EM restarts are performed. The *K* with the largest BIC is selected. The number of EM restarts can be changed with  
`-t 10 [5]`

###Training
The model gets trained with the given number of generators, otherwise *K* has been be chosen by model selection. The number of EM restarts can be changed with   
`-m 30 [50]`
####Convergence
During training of the 50 EM restarts, the convergence criterion is a defined delta loglikelihood  
`-d 1e-4 [1e-6]`

The model with the best likehood gets additional restarts to converge  
`-dd 1e-6 [1e-8]`

###Output
After the model has been trained, 10,000 haplotypes are sampled at the MAP estimate. The distribution of haplotypes is saved in the file `quasispecies.fasta`.  
In addition, a HTML `K?-summary.html` file  visualizes the position wise mutation vectors and recombination rates of the generators. Furthermore, the initial probability vector pi and the estimated error-rate is shown. The file `K?-result.txt` is a text representation of the HTML visualization and `K?-minimal.txt` only shows 