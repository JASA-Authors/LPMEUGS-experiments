# Replication Guide 

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
2. Table2.R will generate the result for the Table 2. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
3. Model 1.R will generate the result for the Table 3. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
4. Model 2.R will generate the result for the Table 4. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
5. Model 3.R will generate the result for the Table 5. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
6. Model 4.R will generate the result for the Table 6. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
7. Model 5.R will generate the result for the Table 7. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.

### figures  

Produces Figures 1-2 in the paper. They are generated by plots.R in the real data file.

1. 1.png is Figure 1 in the paper.
2. 2.png is Figure 2 in the paper.

### data

This file contains dataset for the real data analysis. They are public data from https://bioinformatics.mdanderson.org/public-datasets/ and have been modified for data analysis.

1. label.csv is the true label.
2. data_breast.csv is the dataset.
3. old.zip is the original format of the above two documents.

### simulation-results  

Contains the results from running the code in the simulation-code folder as described above. 


### real-data-code  

The code to run the real data analysis for Figure 1-2 and Table 8. 
1. plots.R is to generate plots.
2. real_data.R is to generate Table 8. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
