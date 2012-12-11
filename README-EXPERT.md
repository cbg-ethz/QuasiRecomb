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
In addition, a HTML `K?-summary.html` file  visualizes the position wise mutation vectors and recombination rates of the generators. Furthermore, the initial probability vector pi and the estimated error-rate is shown. The file `K?-result.txt` is a text representation of the HTML visualization and `K?-minimal.txt` only shows the position, at which mutation or recombination distributions are flat.

The file `support/best.optimum` is the binary java version of the results. 

####Haplotype sampling
The binary java *optimum* files can be used to sample the haplotype distribution with  
`qr --sample -i support/best.optimum`

####Training from intermediate
Furthermore, the binary java *optimum* files can be used as initialization for training  
`qr -i reads.sam -optimum support/best.optimum`
#####Convergence
In combination with `-optimum`, the convergence criterium can be switched to number of parameters changed `-pdelta`. The EM will converge, if no parameter has changed in the last EM step. Otherwise the delta likelihood `-dd 1e-8` can be tweaked.

#Intermediate files
One interesting overview is `support/hit_dist.txt`, which shows the position-wise allel frequencies of the input alignment. This gives a comprehensive overview of the dataset.

During training, intermediate results, as long as they are better then their predecessors, are saved in `support/snapshots/training/`. Files have the naming pattern `R(restart)_K(current generator)_(number of restarts until converged)` and the suffix *.txt* or *.optimum*. The *txt* file has the same structure as the *K?-result.txt*. The *.optimum* file is like the *support/best.optimum*, from which can be sampled or used as initialization.

#DEBUG
Debugging is important, but the normal command line output gives you some insights during execution.  
````
00:00:01:639 Model training [K 2]:          0% [ETA:00:03:01:420][cK 2]
|----------| |------------| |---|  |---------| |----------------| |---|
  Elapsed     Training or   fix K    Elapsed     Remaining time    current K
   time       Modelselect.           restarts    for current EM (only import for pruning)
                                                      step
`````

Additional information can be printed after each EM step with `-verbose -print`. The normal line gets extended by:
````
    5            62     0.04311297852050573    m(11|17)    r(82|98)     232     -11851.217734419906
|--------| |-------|   |-------------------|  |-------|   |--------|  |-------| |------------------|  
   EM      Millisecs   delta log-likelihood   flat mu's   flat rho's   #param    log-likelihood
iteration    taken                                                     changed
````
The fields flat mu's and rho's are of following structure `?(R|E)`, where E is the number of flat distributions of the expected counts and R is the number of flat distributions after regularization. Here, the incremental regularization can be observed.
#ADVANCED
###Backward pruning
Model initialization is a crucial point and one main aspect, why the EM may not find the global optimum. Therefore, backward pruning has been implemented and can be activated by `-prune`. With pruning, each restart for given *K*, starts with *K**2 and then one EM round is done. After that, the most similar generators are merged (Kullback Leibler divergence) or the generator with the highest entropy. This is done until *K* is reached. This provides more stable solutions.