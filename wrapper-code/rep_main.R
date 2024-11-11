#################
#### Table 1 ####
#################
source("Source_Table1.R")
dir.create("./Table1/") # create a new folder to store results
set.seed(1234)

repli=400 # simulation replication times

Method_list=c("Model1","Model2","Model3","Model4","Model5")
for(method in Method_list){
  Source_Table1(method,repli,verbose=F,save_result=T)
}

#################
#### Table 2 ####
#################
library(foreach)
library(doParallel)
clnum<-22
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

source("Source_Table2.R")
dir.create("./Table2/") # create a new folder to store results

repli=3 # simulation replication times

n_p_list=list(c(200,120),c(200,160),c(120,200),c(160,200))
for (nums in n_p_list) {
  Source_Table2(repli,n=nums[1],p=nums[2],verbose=T,save_result=T)
}

stopCluster(cl)

#################
#### Table 3 ####
#################
library(foreach)
library(doParallel)
clnum<-24
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

source("Source_Table3.R")
dir.create("./Table3/") # create a new folder to store results

repli=3 # simulation replication times

n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
for (nums in n_p_list) {
  Source_Table3(repli,n=nums[1],p=nums[2],verbose=T,save_result=T)
}

stopCluster(cl)

#################
#### Table 4 ####
#################
library(foreach)
library(doParallel)
clnum<-24
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

source("Source_Table4.R")
dir.create("./Table4/") # create a new folder to store results

repli=3 # simulation replication times

n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
for (nums in n_p_list) {
  Source_Table4(repli,n=nums[1],p=nums[2],verbose=T,save_result=T)
}

stopCluster(cl)

#################
#### Table 5 ####
#################
library(foreach)
library(doParallel)
clnum<-24
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

source("Source_Table5.R")
dir.create("./Table5/") # create a new folder to store results

repli=3 # simulation replication times

n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
for (nums in n_p_list) {
  Source_Table5(repli,n=nums[1],p=nums[2],verbose=T,save_result=T)
}

stopCluster(cl)

#################
#### Table 6 ####
#################
library(foreach)
library(doParallel)
clnum<-24
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

source("Source_Table6.R")
dir.create("./Table6/") # create a new folder to store results

repli=3 # simulation replication times

n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
for (nums in n_p_list) {
  Source_Table6(repli,n=nums[1],p=nums[2],verbose=T,save_result=T)
}

stopCluster(cl)

#################
#### Table 7 ####
#################
library(foreach)
library(doParallel)
clnum<-24
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

source("Source_Table7.R")
dir.create("./Table7/") # create a new folder to store results

repli=3 # simulation replication times

n_p_list=list(c(200,120),c(200,160),c(50,200),c(120,200),c(160,200))
for (nums in n_p_list) {
  Source_Table7(repli,n=nums[1],p=nums[2],verbose=T,save_result=T)
}

stopCluster(cl)

##################
#### Figure 1 ####
##################
set.seed(1234)
source("Source_Figure1.R")
dir.create("./Figure1/") # create a new folder to store results
Source_Figure1(save_result=T)

##################
#### Figure 2 ####
##################
set.seed(1234)
source("Source_Figure2.R")
dir.create("./Figure2/") # create a new folder to store results
Source_Figure2(save_result=T)

#################
#### Table 8 ####
#################
library(foreach)
library(doParallel)
clnum<-24
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

source("Source_Table8.R")
dir.create("./Table8/") # create a new folder to store results

repli=3 # simulation replication times

Source_Table8(repli,verbose=T,save_result=T)


stopCluster(cl)
