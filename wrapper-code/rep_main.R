library(foreach)
library(doParallel)
set.seed(1234)

### Indicators for running specific tables or figures
### Change the indicator to FALSE to skip an experiment.
run_table1=TRUE
run_table2=TRUE
run_table3=TRUE
run_table4=TRUE
run_table5=TRUE
run_table6=TRUE
run_table7=TRUE
run_figure1=TRUE
run_figure2=TRUE
run_table8=TRUE


### Number of simulation replications for each experiment
### The default numbers allow the user to replicate all the results on a personal computer within a managable time
repli_table1=400
repli_table2=2
repli_table3=2
repli_table4=2
repli_table5=2
repli_table6=2
repli_table7=2
repli_table8=2

### Other parameters
verbose=F # whether to print results
save_result=T # whether to save results



#################
#### Table 1 ####
#################
if(run_table1==TRUE){
  source("Source_Table1.R")
  dir.create("./Table1/") # create a new folder to store results
  
  Method_list=c("Model1","Model2","Model3","Model4","Model5")
  
  ## Function: Source_Table1(Method,repli=400,verbose=T,save_result=T)
  ## Input:
  ## Method: Model1, Model2, Model3, Model4, Model5
  ## Definition of models is from Section 4.1
  ## repli: simulation replication times
  ## verbose: whether to print the result table for the selected model
  ## default is True.
  ## save_result: whether to save result table to be a CSV file
  ## default is True.
  ##
  ## Output:
  ## result: result table for the selected model
  for(method in Method_list){
    Source_Table1(method,repli_table1,verbose,save_result)
  }
}

#################
#### Table 2 ####
#################
if(run_table2==TRUE){
  clnum<-22
  cl <- makeCluster(getOption("cl.cores", clnum))
  registerDoParallel(cl)
  
  source("Source_Table2.R")
  dir.create("./Table2/") # create a new folder to store results
  
  n_p_list=list(c(200,120),c(200,160),c(120,200),c(160,200))
  
  ## Function: Source_Table2(repli=5,n=200,p=160,verbose=T,save_result=T)
  ## Input:
  ## repli: simulation replication times
  ## n: sample size
  ## p: number of features
  ## verbose: whether to print the result table for the selected model
  ## default is True.
  ## save_result: whether to save result table to be a CSV file
  ## default is True.
  ##
  ## Output:
  ## result: result table for the selected model
  for (nums in n_p_list) {
    Source_Table2(repli_table2,n=nums[1],p=nums[2],verbose,save_result)
  }
  
  stopCluster(cl)
}

#################
#### Table 3 ####
#################
if(run_table3==TRUE){
  clnum<-24
  cl <- makeCluster(getOption("cl.cores", clnum))
  registerDoParallel(cl)
  
  source("Source_Table3.R")
  dir.create("./Table3/") # create a new folder to store results
  
  ## Function: Source_Table3(repli=5,n=200,p=160,verbose=T,save_result=T)
  ## Input:
  ## repli: simulation replication times
  ## n: sample size
  ## p: number of features
  ## verbose: whether to print the result table for the selected model
  ## default is True.
  ## save_result: whether to save result table to be a CSV file
  ## default is True.
  ##
  ## Output:
  ## result: result table for the selected model
  n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
  for (nums in n_p_list) {
    Source_Table3(repli_table3,n=nums[1],p=nums[2],verbose,save_result)
  }
  
  stopCluster(cl)
}


#################
#### Table 4 ####
#################
if(run_table4==TRUE){
  clnum<-24
  cl <- makeCluster(getOption("cl.cores", clnum))
  registerDoParallel(cl)
  
  source("Source_Table4.R")
  dir.create("./Table4/") # create a new folder to store results
  
  ## Function: Source_Table4(repli=5,n=200,p=160,verbose=T,save_result=T)
  ## Input:
  ## repli: simulation replication times
  ## n: sample size
  ## p: number of features
  ## verbose: whether to print the result table for the selected model
  ## default is True.
  ## save_result: whether to save result table to be a CSV file
  ## default is True.
  ##
  ## Output:
  ## result: result table for the selected model
  n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
  for (nums in n_p_list) {
    Source_Table4(repli_table4,n=nums[1],p=nums[2],verbose,save_result)
  }
  
  stopCluster(cl)
}


#################
#### Table 5 ####
#################
if(run_table5==TRUE){
  clnum<-24
  cl <- makeCluster(getOption("cl.cores", clnum))
  registerDoParallel(cl)
  
  source("Source_Table5.R")
  dir.create("./Table5/") # create a new folder to store results
  
  ## Function: Source_Table5(repli=5,n=200,p=160,verbose=T,save_result=T)
  ## Input:
  ## repli: simulation replication times
  ## n: sample size
  ## p: number of features
  ## verbose: whether to print the result table for the selected model
  ## default is True.
  ## save_result: whether to save result table to be a CSV file
  ## default is True.
  ##
  ## Output:
  ## result: result table for the selected model
  n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
  for (nums in n_p_list) {
    Source_Table5(repli_table5,n=nums[1],p=nums[2],verbose,save_result)
  }
  
  stopCluster(cl)
}


#################
#### Table 6 ####
#################
if(run_table6==TRUE){
  clnum<-24
  cl <- makeCluster(getOption("cl.cores", clnum))
  registerDoParallel(cl)
  
  source("Source_Table6.R")
  dir.create("./Table6/") # create a new folder to store results
  
  ## Function: Source_Table6(repli=5,n=200,p=160,verbose=T,save_result=T)
  ## Input:
  ## repli: simulation replication times
  ## n: sample size
  ## p: number of features
  ## verbose: whether to print the result table for the selected model
  ## default is True.
  ## save_result: whether to save result table to be a CSV file
  ## default is True.
  ##
  ## Output:
  ## result: result table for the selected model
  n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
  for (nums in n_p_list) {
    Source_Table6(repli_table6,n=nums[1],p=nums[2],verbose,save_result)
  }
  
  stopCluster(cl)
}


#################
#### Table 7 ####
#################
if(run_table7==TRUE){
  clnum<-24
  cl <- makeCluster(getOption("cl.cores", clnum))
  registerDoParallel(cl)
  
  source("Source_Table7.R")
  dir.create("./Table7/") # create a new folder to store results
  
  ## Function: Source_Table7(repli=5,n=200,p=160,verbose=T,save_result=T)
  ## Input:
  ## repli: simulation replication times
  ## n: sample size
  ## p: number of features
  ## verbose: whether to print the result table for the selected model
  ## default is True.
  ## save_result: whether to save result table to be a CSV file
  ## default is True.
  ##
  ## Output:
  ## result: result table for the selected model
  n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
  for (nums in n_p_list) {
    Source_Table7(repli_table7,n=nums[1],p=nums[2],verbose,save_result)
  }
  
  stopCluster(cl)
}


##################
#### Figure 1 ####
##################
if(run_figure1==TRUE){
  source("Source_Figure1.R")
  dir.create("./Figure1/") # create a new folder to store results
  
  ## Function: Source_Figure1(save_result=T)
  ## Input:
  ## save_result: whether to save result Figure to be a png file
  ## default is True.
  ##
  ## Output:
  ## Figure 1 in the paper
  Source_Figure1(save_result)
}


##################
#### Figure 2 ####
##################
if(run_figure2==TRUE){
  source("Source_Figure2.R")
  dir.create("./Figure2/") # create a new folder to store results
  
  ## Function: Source_Figure2(save_result=T)
  ## Input:
  ## save_result: whether to save result Figure to be a png file
  ## default is True.
  ##
  ## Output:
  ## Figure 2 in the paper
  Source_Figure2(save_result)
}


#################
#### Table 8 ####
#################
if(run_table8==TRUE){
  clnum<-24
  cl <- makeCluster(getOption("cl.cores", clnum))
  registerDoParallel(cl)
  
  source("Source_Table8.R")
  dir.create("./Table8/") # create a new folder to store results
  
  ## Function: Source_Table8(repli=5,verbose=T,save_result=T)
  ## Input:
  ## repli: simulation replication times
  ## verbose: whether to print the result table for the selected model
  ## default is True.
  ## save_result: whether to save result table to be a CSV file
  ## default is True.
  ##
  ## Output:
  ## result: result table
  Source_Table8(repli_table8,verbose,save_result)
  
  
  stopCluster(cl)
}
