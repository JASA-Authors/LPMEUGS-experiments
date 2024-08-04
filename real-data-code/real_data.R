library (glmnet)
library(glasso)
library(expm)
library(flare) # CLIME, TIGER
library(POET)
library(foreach)
library(doParallel)
clnum<-24
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

set.seed(1234)

PCAkmeans=function(X, K, lambda1, lambda2, theta){
  # OUR method
  # This function is slight different with Model 1.R and add cluster information and delete the input of true covariance matrix
  # input:
  # X: data
  # K: number of groups (must be 4)
  # theta: judge when to stop
  # lambda1: parameter for adapted huber regression
  # lambda2: parameter for glasso regression
  
  # output:
  # Omega: estimated inversed covariance matrix
  # orderB: estimated cluster information
  
  n=nrow(X)
  p=ncol(X)
  
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf=kmeans(V, K)$cluster
  
  B_count=c()
  for (i in 1:K) {
    B_count=c(B_count,length(which(clusterinf==i)))
  }
  seq1=order(B_count,decreasing = FALSE)
  
  Bcluster=list()
  C=c()
  for (i in seq1) {
    Bcluster=c(Bcluster,list(which(clusterinf==i)))
    C=c(C,length(which(clusterinf==i)))
  }
  B=matrix(0,nrow = n,ncol = p)
  B[,1:C[1]]=X[,Bcluster[[1]]]
  for (i in 2:K) {
    B[, (sum(C[1:(i-1)])+1) : sum(C[1:i])] = X[,Bcluster[[i]]]
  }
  
  inverD=matrix(0,nrow=p,ncol=p)
  A=matrix(0,nrow=p,ncol=p)
  
  # for block 1
  X1=B[,1:C[1]]
  S1=t(X1)%*%(X1)
  S1=S1/n
  inverD1=glasso::glasso(S1,rho = lambda2)$wi
  inverD[1:C[1],1:C[1]]=inverD1
  
  # for other blocks
  for (i in 2:K) {
    numrepeat=1 # number of repeat times
    Zi=B[,1:sum(C[1:(i-1)])]
    Xi=B[,(sum(C[1:(i-1)])+1):(sum(C[1:i]))]
    
    # initial Ai0 (Di0^-1 is indentity matrix)
    wildX=c(Xi%*%diag(C[i]))
    wildX=as.matrix(wildX)
    
    wildZ=as.matrix(kronecker(diag(C[i]),Zi))
    
    if(ncol(wildZ)==1){
      lasso.mod=lm( wildX ~ 0+ wildZ)$coefficients
      newAi=t(matrix(lasso.mod,ncol=C[i]))
    }else{
      lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
      newAi=t(matrix(coef(lasso.mod)[-1],ncol=C[i]))
    }
    
    # initial Di1^-1 by Ai0
    Si=t(Xi-Zi%*%t(newAi))%*%(Xi-Zi%*%t(newAi))
    Si=Si/n
    newinverDi=glasso::glasso(Si,rho = lambda2)$wi
    
    # second time to calculate Ai and Di
    numrepeat=numrepeat+1
    oldAi=newAi
    oldinverDi=newinverDi
    wildX=c(Xi%*%expm::sqrtm(oldinverDi))
    wildX=as.matrix(wildX)
    
    wildZ=as.matrix(kronecker(expm::sqrtm(oldinverDi),Zi))
    
    if(ncol(wildZ)==1){
      lasso.mod=lm( wildX ~ 0+ wildZ)$coefficients
      newAi=t(matrix(lasso.mod,ncol=C[i]))
    }else{
      lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
      newAi=t(matrix(coef(lasso.mod)[-1],ncol=C[i]))
    }
    
    Si=t(Xi-Zi%*%t(newAi))%*%(Xi-Zi%*%t(newAi))
    Si=Si/n
    newinverDi=glasso::glasso(Si,rho = lambda2)$wi
    
    # judge whether converge
    conAi=sqrt(sum((oldAi-newAi)^2))
    coninverDi=sqrt(sum((oldinverDi-newinverDi)^2))
    
    while (numrepeat<100) {
      numrepeat=numrepeat+1
      oldAi=newAi
      oldinverDi=newinverDi
      wildX=c(Xi%*%expm::sqrtm(oldinverDi))
      wildX=as.matrix(wildX)
      
      wildZ=as.matrix(kronecker(expm::sqrtm(oldinverDi),Zi))
      
      if(ncol(wildZ)==1){
        lasso.mod=lm( wildX ~ 0+ wildZ)$coefficients
        newAi=t(matrix(lasso.mod,ncol=C[i]))
      }else{
        lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
        newAi=t(matrix(coef(lasso.mod)[-1],ncol=C[i]))
      }
      
      Si=t(Xi-Zi%*%t(newAi))%*%(Xi-Zi%*%t(newAi))
      Si=Si/n
      newinverDi=glasso::glasso(Si,rho = lambda2)$wi
      
      # judge whether converge
      conAi=sqrt(sum((oldAi-newAi)^2))
      coninverDi=sqrt(sum((oldinverDi-newinverDi)^2))
      
      if(conAi<theta && coninverDi<theta){
        break
      }
      
    }
    
    # insert Ai to A and inverDi to inverD
    A[(sum(C[1:(i-1)])+1):(sum(C[1:i])),1:(sum(C[1:(i-1)]))]=newAi
    inverD[(sum(C[1:(i-1)])+1):(sum(C[1:i])),(sum(C[1:(i-1)])+1):(sum(C[1:i]))]=newinverDi
  }
  
  hatT=diag(p)-A
  newOmega=t(hatT)%*%inverD%*%hatT
  
  return(list(Omega=newOmega,orderB=Bcluster))
}

BCDnormal=function(X,K,lambda1,lambda2,theta){
  # K=p is NO-GROUP method 
  n=nrow(X)
  p=ncol(X)
  numblock=K
  sizeblock=p/numblock
  
  inverD=matrix(0,nrow=p,ncol=p)
  A=matrix(0,nrow=p,ncol=p)
  # for block 1
  X1=X[,1:sizeblock]
  S1=t(X1)%*%(X1)
  S1=S1/n
  inverD1=glasso::glasso(S1,rho = lambda2)$wi
  inverD[1:sizeblock,1:sizeblock]=inverD1
  
  
  # for other blocks
  for (i in 2:numblock) {
    numrepeat=1 # number of repeat times
    Zi=X[,1:((i-1)*sizeblock)]
    Xi=X[,((i-1)*sizeblock+1):(i*sizeblock)]
    
    # initial Ai0 (Di0^-1 is indentity matrix)
    wildX=c(Xi%*%diag(sizeblock))
    wildX=as.vector(wildX)
    
    wildZ=kronecker(diag(sizeblock),Zi)
    
    if(ncol(wildZ)==1){
      lasso.mod=lm( wildX ~ 0+ wildZ)$coefficients
      newAi=t(matrix(lasso.mod,ncol=sizeblock))
    }else{
      lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
      newAi=t(matrix(coef(lasso.mod)[-1],ncol=sizeblock))
    }
    
    # initial Di1^-1 by Ai0
    Si=t(Xi-Zi%*%t(newAi))%*%(Xi-Zi%*%t(newAi))
    Si=Si/n
    newinverDi=glasso::glasso(Si,rho = lambda2)$wi
    
    # second time to calculate Ai and Di
    numrepeat=numrepeat+1
    oldAi=newAi
    oldinverDi=newinverDi
    wildX=c(Xi%*%expm::sqrtm(oldinverDi))
    wildX=as.vector(wildX)
    
    wildZ=kronecker(expm::sqrtm(oldinverDi),Zi)
    
    if(ncol(wildZ)==1){
      lasso.mod=lm( wildX ~ 0+ wildZ)$coefficients
      newAi=t(matrix(lasso.mod,ncol=sizeblock))
    }else{
      lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
      newAi=t(matrix(coef(lasso.mod)[-1],ncol=sizeblock))
    }
    
    Si=t(Xi-Zi%*%t(newAi))%*%(Xi-Zi%*%t(newAi))
    Si=Si/n
    newinverDi=glasso::glasso(Si,rho = lambda2)$wi
    
    # judge whether converge
    conAi=sqrt(sum((oldAi-newAi)^2))
    coninverDi=sqrt(sum((oldinverDi-newinverDi)^2))
    
    while (numrepeat<100) {
      if(conAi<theta && coninverDi<theta){
        break
      }
      
      numrepeat=numrepeat+1
      oldAi=newAi
      oldinverDi=newinverDi
      wildX=c(Xi%*%expm::sqrtm(oldinverDi))
      wildX=as.vector(wildX)
      
      wildZ=kronecker(expm::sqrtm(oldinverDi),Zi)
      
      if(ncol(wildZ)==1){
        lasso.mod=lm( wildX ~ 0+ wildZ)$coefficients
        newAi=t(matrix(lasso.mod,ncol=sizeblock))
      }else{
        lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
        newAi=t(matrix(coef(lasso.mod)[-1],ncol=sizeblock))
      }
      
      Si=t(Xi-Zi%*%t(newAi))%*%(Xi-Zi%*%t(newAi))
      Si=Si/n
      newinverDi=glasso::glasso(Si,rho = lambda2)$wi
      
      conAi=sqrt(sum((oldAi-newAi)^2))
      coninverDi=sqrt(sum((oldinverDi-newinverDi)^2))
      
    }
    
    # insert Ai to A and inverDi to inverD
    A[((i-1)*sizeblock+1):(i*sizeblock),1:((i-1)*sizeblock)]=newAi
    inverD[((i-1)*sizeblock+1):(i*sizeblock),((i-1)*sizeblock+1):(i*sizeblock)]=newinverDi
  }
  
  hatT=diag(p)-A
  newOmega=t(hatT)%*%inverD%*%hatT
  
  return(list(Omega=newOmega))
}

judge_label=function(mu0,mu1,pi0,pi1,test_data,Omega,true_lable){
  # this function is calculate Specificity, Sensitivity, and Matthews Correlation Coefficien (MCC)
  value0=t(test_data)%*%Omega%*%mu0
  value0=value0-matrix(rep(t(mu0)%*%Omega%*%mu0/2+log(pi0),33),ncol=1)
  value1=t(test_data)%*%Omega%*%mu1
  value1=value1-matrix(rep(t(mu1)%*%Omega%*%mu1/2+log(pi1),33),ncol=1)
  final=cbind(value0,value1)
  label_pred=apply(final,1,which.max)-1
  
  precision_table=table(label_pred,true_lable)
  Specificity=precision_table[1,1]/(precision_table[1,1]+precision_table[2,1])
  Sensitivity=precision_table[2,2]/(precision_table[2,2]+precision_table[1,2])
  MCC=(precision_table[1,1]*precision_table[2,2]-precision_table[1,2]*precision_table[2,1])/(
    sqrt(precision_table[1,1]+precision_table[1,2])*sqrt(precision_table[1,1]+precision_table[2,1])*
      sqrt(precision_table[2,2]+precision_table[2,1])*sqrt(precision_table[2,2]+precision_table[1,2])
  )
  
  return(c(Specificity,Sensitivity,MCC))
}

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

MDDA=function(x1){
  RD_test_index=sample(1:n_RD,size = 25,replace = FALSE)
  pCR_test_index=sample(1:n_pCR,size = 8,replace = FALSE)
  RD_train_data=data_breast_RD[,-RD_test_index]
  pCR_train_data=data_breast_pCR[,-pCR_test_index]
  
  true_lable=c(rep(0,25),rep(1,8))
  
  
  p_value_seq=numeric(p)
  for (i in 1:p) {
    p_value_seq[i]=t.test(RD_train_data[i,], pCR_train_data[i,], var.equal=TRUE)$p.value
  }
  significant_genes=order(p_value_seq,decreasing=FALSE)[1:101]
  
  RD_train_data=RD_train_data[significant_genes,]
  pCR_train_data=pCR_train_data[significant_genes,]
  train_data=cbind(RD_train_data,pCR_train_data)
  sd_train=apply(train_data,1,sd)
  for (i in 1:ncol(train_data)) {
    train_data[,i]=train_data[,i]/sd_train
  }
  
  X=t(train_data)
  test_data=cbind(data_breast_RD[,RD_test_index],data_breast_pCR[,pCR_test_index])
  test_data=test_data[significant_genes,]
  for (i in 1:ncol(test_data)) {
    test_data[,i]=test_data[,i]/sd_train
  }
  
  pi0=ncol(RD_train_data)/(ncol(RD_train_data)+ncol(pCR_train_data))
  pi1=ncol(pCR_train_data)/(ncol(RD_train_data)+ncol(pCR_train_data))
  mu0=matrix(rowMeans(train_data[,1:ncol(RD_train_data)]),ncol=1)
  mu1=matrix(rowMeans(train_data[,(ncol(RD_train_data)+1):ncol(train_data)]),ncol=1)
  
  est1=PCAkmeans(X,2,0.01,0.01,0.01)
  Omegahat1=est1$Omega
  ttt=unlist(est1$orderB)
  
  result_1=judge_label(as.matrix(mu0[ttt],ncol=1),as.matrix(mu1[ttt],ncol=1),pi0,pi1,test_data[ttt,],Omegahat1,true_lable)
  
  Omegahat2=BCDnormal(X,ncol(X),0.01,0.01,0.01)$Omega
  
  result_2=judge_label(mu0,mu1,pi0,pi1,test_data,Omegahat2,true_lable)
  
  Omegahat3=glasso::glasso(cov(X), rho=.1)$wi
  
  result_3=judge_label(mu0,mu1,pi0,pi1,test_data,Omegahat3,true_lable)
  
  Omegahat4=flare::sugm(X,nlambda=1,method = "tiger")$icov[[1]]
  
  result_4=judge_label(mu0,mu1,pi0,pi1,test_data,Omegahat4,true_lable)
  
  Omegahat5=flare::sugm(X,nlambda=1,method = "clime")$icov[[1]]
  
  result_5=judge_label(mu0,mu1,pi0,pi1,test_data,Omegahat5,true_lable)
  
  Omegahat6=solve(POET::POET(t(X),2)$SigmaY)
  
  result_6=judge_label(mu0,mu1,pi0,pi1,test_data,Omegahat6,true_lable)
  
  result=c(result_1,result_2,result_3,result_4,result_5,result_6)
  result
}



result1=foreach(i=1:100,.combine='rbind') %dopar% MDDA(i)
write.csv(result1,"real_data_result.csv")

stopCluster(cl)
