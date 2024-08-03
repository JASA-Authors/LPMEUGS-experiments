# Reproducibility Materials

This GitHub repository contains codes and data to reproduce simulation results, figures, and real data analysis results from the paper "Large Precision Matrix Estimation with Unknown Group Structure".



## Package Environments

The codes in this repository requires R version 4.3.1 and the following R packages:

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


### simulation-demo 

This folder contains the six demo codes to run the simulation examples introduced in Section 4 and Appendix B in the paper. The demo codes will run every model and setting for one replication (you can also change the replication times by changing the variable ```replication_time```) for illustration purpose. The expected running time for each document on a regular desktop PC should be within 5 minutes. 

### simulation-full  

This folder contains the codes to reproduce the full simulation results presented in Section 4 and Appendix B in the paper. As the full simulation results invovles multiple competing methods, various settings, and 200 replications, the codes utilized  parallel computing in order to reduce the running time. We would recommend to reproduce the full simulation results with a computing system of at least 24 cores and 500GB of storage. The following list provides the guidance for each file in the folder.

1. Table1.R will generate the result for the Table 1.
2. Table2.R will generate the result for the Table 2. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
3. Model 1.R will generate the result for the Table 3. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
4. Model 2.R will generate the result for the Table 4. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
5. Model 3.R will generate the result for the Table 5. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
6. Model 4.R will generate the result for the Table 6. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
7. Model 5.R will generate the result for the Table 7. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.
8. Model 6.R will generate the result for the Table 1 in the Appendix B. The codes will generate the error for each method in each time simulation. Remember to use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std.


### simulation-results  

This folder contains the results from running the code in the simulation-code folder as described above. 


### data

This folder contains the real dataset dataset for the real data analysis. The dataset is publicly avaliable at  https://bioinformatics.mdanderson.org/public-datasets/ . Please download the data file named "MDA133: Clinical Data and dChip MBEI value Files" from this link. We have modified the dataset for data analysis.

1. "label.csv" is the true label.
2. "data_breast.csv" is the dataset.
3. old.zip is the original format of the above two documents. In this file, there are two excel documents named "MDA133CompleteInfo20070319.xls" and "MDA133PredictorTrainAndValidation.xls". Here is how we modified these two datasets:
    *  We used the column "pCR" in "MDA133CompleteInfo20070319.xls" as our label which is the "label.csv".
    *  We deleted the header and genes name in "MDA133PredictorTrainAndValidation.xls" and built our "data_breast.csv".


### real-data-code  

This folder contrains the codes to reproduce the real data analysis results presented in Section 5 in the paper. 
1. plots.R will generate two plots. 1.png is Figure 1 in the paper. 2.png is Figure 2 in the paper.
2. real_data.R will generate the analysis results presented in Table 8. The codes will generate the error for each method in each time simulation. Please use ```apply(x,2,mean)``` and ```apply(x,2,sd)``` to calculate error mean and std (include it inthe code?).

### figures  

This folder contains Figures 1-2 in the paper. They are generated by plots.R in the real data file. The output contains the following two files.

1. 1.png is Figure 1 in the paper.
2. 2.png is Figure 2 in the paper.
