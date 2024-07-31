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

Besides, these codes use parallel computing in order to reduce simulation time and require at least 24 cores. 

Table1.R will generate the result for the Table 1.

