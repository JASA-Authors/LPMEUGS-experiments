# BCD-experiments 

Code to reproduce simulation results, figures, and real data analysis results from the paper "Large Precision Matrix Estimation with Unknown Group Structure".

## Organization

### simulation-code  

The codes in this file requires R version 4.4.1 and the following packages:

1. ```MASS 7.3.61```
2. ```pdfCluster 1.04```
3. ```Matrix 1.7.0```
4. ```glmnet 4.1.8```
5. ```glasso 1.11```
6. ```expm 0.999.9```
7. ```flare 1.7.0.1```
8. ```POET 2.0```
9. ```foreach 1.5.2```
10. ```doParallel 1.0.17```

Besides, these codes use parallel computing in order to reduce simulation time and require at least 24 cores and 500GB of storage. The following is the guidance for each document.

1. Table1.R will generate the result for the Table 1.
2. 

### figures  

Produces Figures 1-2 in the paper. They are generated by plots.R in the real data file.

1. 1.png is Figure 1 in the paper.
2. 2.png is Figure 2 in the paper.

### data

This file contains dataset for the real data analysis. They are public data from https://bioinformatics.mdanderson.org/public-datasets/ and have been modified for data analysis.

1. label.csv is the true label.
2. data_breast.csv is the dataset.
3. old.zip is the original format of the above two documents.

### real-data-code  

The code to run the real data analysis for Figure 1-2 and Table 8. 
1. plots.R is to generate plots.
2. real_data.R is to generate Table 8.
