require(MASS)
require(glmnet)
require(glasso)
require(expm)
require(flare) # CLIME, TIGER
require(POET)
require(Matrix)
require(ggplot2)
require(patchwork)

Source_Figure1=function(save_result=T){
  ## Input:
  ## save_result: whether to save result Figure to be a png file
  ## default is True.
  ##
  ## Output:
  ## Figure 1 in the paper
  
  ############
  ### Data ###
  ############
  data_breast=read.csv("data_breast.csv",header = T)
  data_label=read.csv("label.csv",header = T)
  data_label_RD=which(data_label$pCR==0)
  data_label_pCR=which(data_label$pCR==1)
  data_breast_RD=data_breast[,data_label_RD]
  data_breast_pCR=data_breast[,data_label_pCR]
  
  
  n=ncol(data_breast)
  n_RD=length(data_label_RD)
  n_pCR=length(data_label_pCR)
  p=nrow(data_breast)
  
  RD_test_index=sample(1:n_RD,size = 25,replace = FALSE)
  pCR_test_index=sample(1:n_pCR,size = 8,replace = FALSE)
  RD_train_data=data_breast_RD[,-RD_test_index]
  pCR_train_data=data_breast_pCR[,-pCR_test_index]
  
  true_lable=c(rep(0,8),rep(1,25))
  
  ####################
  ### Select Genes ###
  ####################
  p_value_seq=numeric(p)
  for (i in 1:p) {
    p_value_seq[i]=t.test(RD_train_data[i,], pCR_train_data[i,], var.equal=TRUE)$p.value
  }
  significant_genes=order(p_value_seq,decreasing=FALSE)[1:101]
  
  RD_train_data=RD_train_data[significant_genes,]
  pCR_train_data=pCR_train_data[significant_genes,]
  train_data=cbind(pCR_train_data,RD_train_data)
  sd_train=apply(train_data,1,sd)
  for (i in 1:ncol(train_data)) {
    train_data[,i]=train_data[,i]/sd_train
  }
  
  X=t(train_data)
  
  #############
  ### Plots ###
  #############
  ## Screen Plot and Projection
  pca <- prcomp(t(X), scale = TRUE)
  variance = pca $sdev^2 / sum(pca $sdev^2)
  p1=ggplot() + 
    aes(c(1:10), variance[1:10])+
    geom_line() + 
    geom_point(size=4)+
    xlab("Principal Components") + 
    ylab("Eigenvalues") +
    ylim(0, 0.7)+ scale_x_continuous(breaks=1:10)
  
  pcaData <- as.data.frame(pca$x[, 1:2])
  colnames(pcaData) <- c("PC1", "PC2") 
  p2=ggplot(pcaData) +theme(aspect.ratio=1)+
    aes(PC1, PC2) + # define plot area
    geom_point(size = 2) + # adding data points
    coord_fixed()+ # fixing coordinates
    xlab("First Principle Component") + 
    ylab("Second Principle Component")
  
  p=p1+p2
  print(p)
  
  if(save_result==TRUE){
    ggsave("./Figure1/Figure1.png")
  }
}
