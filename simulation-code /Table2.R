library(MASS)
library (glmnet)
library(glasso)
library(expm)
library(flare) # CLIME, TIGER
library(POET)
library(foreach)
library(doParallel)

clnum<-detectCores()-10 # should be 22
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

set.seed(1234)

################
### Function ###
################
PCAkmeans=function(X, K, lambda1, lambda2, theta, Sigmatrue){
  # OUR method
  # input:
  # X: data
  # K: number of groups 
  # theta: judge when to stop
  # lambda1: parameter for adapted huber regression
  # lambda2: parameter for glasso regression
  
  # output:
  # Omega: estimated inversed covariance matrix
  # error: estimation error
  
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
  
  Bcluster=c()
  for (i in seq1) {
    Bcluster=c(Bcluster,which(clusterinf==i))
  }
  new_Sigmatrue=Sigmatrue[,Bcluster]
  new_Sigmatrue=new_Sigmatrue[Bcluster,]
  
  error=sqrt(sum((solve(new_Sigmatrue)-newOmega)^2))
  
  
  return(list(Omega=newOmega, error=error))
}

PCAkmeans_random=function(X, K, lambda1, lambda2, theta, Sigmatrue){
  # OUR with random method
  # input:
  # X: data
  # K: number of groups 
  # theta: judge when to stop
  # lambda1: parameter for adapted huber regression
  # lambda2: parameter for glasso regression
  
  # output:
  # Omega: estimated inversed covariance matrix
  # error: estimation error
  
  n=nrow(X)
  p=ncol(X)
  
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf=kmeans(V, K)$cluster
  
  Bcluster=list()
  C=c()
  for (i in 1:K) {
    Bcluster=c(Bcluster,list(which(clusterinf==i)))
    C=c(C,length(Bcluster[[i]]))
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
  
  Bcluster=c()
  for (i in 1:K) {
    Bcluster=c(Bcluster,which(clusterinf==i))
  }
  new_Sigmatrue=Sigmatrue[,Bcluster]
  new_Sigmatrue=new_Sigmatrue[Bcluster,]
  
  error=sqrt(sum((solve(new_Sigmatrue)-newOmega)^2))
  
  
  return(list(Omega=newOmega, error=error))
}

PCAkmeans_decrease=function(X, K, lambda1, lambda2, theta, Sigmatrue){
  # OUR with decrease method
  # input:
  # X: data
  # K: number of groups (must be 4)
  # theta: judge when to stop
  # lambda1: parameter for adapted huber regression
  # lambda2: parameter for glasso regression
  
  # output:
  # Omega: estimated inversed covariance matrix
  # error: estimation error
  
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
  seq1=order(B_count,decreasing = TRUE)
  
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
  
  Bcluster=c()
  for (i in seq1) {
    Bcluster=c(Bcluster,which(clusterinf==i))
  }
  new_Sigmatrue=Sigmatrue[,Bcluster]
  new_Sigmatrue=new_Sigmatrue[Bcluster,]
  
  error=sqrt(sum((solve(new_Sigmatrue)-newOmega)^2))
  
  
  return(list(Omega=newOmega, error=error))
}

BCDnormal2=function(X,K,lambda1,lambda2,theta, Sigmatrue){
  # ORACLE method
  n=nrow(X)
  p=ncol(X)
  numblock=K
  sizeblock=c(0.1*p,0.2*p,0.3*p,0.4*p)
  
  inverD=matrix(0,nrow=p,ncol=p)
  A=matrix(0,nrow=p,ncol=p)
  # for block 1
  X1=X[,1:sizeblock[1]]
  S1=t(X1)%*%(X1)
  S1=S1/n
  inverD1=glasso::glasso(S1,rho = lambda2)$wi
  inverD[1:sizeblock[1],1:sizeblock[1]]=inverD1
  
  
  # for other blocks
  for (i in 2:numblock) {
    numrepeat=1 # number of repeat times
    Zi=X[,1:sum(sizeblock[1:(i-1)])]
    Xi=X[,(sum(sizeblock[1:(i-1)])+1):(sum(sizeblock[1:i]))]
    
    # initial Ai0 (Di0^-1 is indentity matrix)
    wildX=c(Xi%*%diag(sizeblock[i]))
    wildX=as.vector(wildX)
    
    wildZ=kronecker(diag(sizeblock[i]),Zi)
    
    if(ncol(wildZ)==1){
      lasso.mod=lm( wildX ~ 0+ wildZ)$coefficients
      newAi=t(matrix(lasso.mod,ncol=sizeblock[i]))
    }else{
      lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
      newAi=t(matrix(coef(lasso.mod)[-1],ncol=sizeblock[i]))
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
      newAi=t(matrix(lasso.mod,ncol=sizeblock[i]))
    }else{
      lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
      newAi=t(matrix(coef(lasso.mod)[-1],ncol=sizeblock[i]))
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
        newAi=t(matrix(lasso.mod,ncol=sizeblock[i]))
      }else{
        lasso.mod=glmnet::glmnet (wildZ, wildX, alpha=1, lambda=lambda1,intercept=FALSE)
        newAi=t(matrix(coef(lasso.mod)[-1],ncol=sizeblock[i]))
      }
      
      Si=t(Xi-Zi%*%t(newAi))%*%(Xi-Zi%*%t(newAi))
      Si=Si/n
      newinverDi=glasso::glasso(Si,rho = lambda2)$wi
      
      conAi=sqrt(sum((oldAi-newAi)^2))
      coninverDi=sqrt(sum((oldinverDi-newinverDi)^2))
      
    }
    
    # insert Ai to A and inverDi to inverD
    A[(sum(sizeblock[1:(i-1)])+1):(sum(sizeblock[1:i])),1:(sum(sizeblock[1:(i-1)]))]=newAi
    inverD[(sum(sizeblock[1:(i-1)])+1):(sum(sizeblock[1:i])),(sum(sizeblock[1:(i-1)])+1):(sum(sizeblock[1:i]))]=newinverDi
  }
  
  hatT=diag(p)-A
  newOmega=t(hatT)%*%inverD%*%hatT
  serror=sqrt(sum((solve(Sigmatrue)-newOmega)^2))
  
  return(list(Omega=newOmega,error=serror))
}


##########################
### n=160 p=200        ###
##########################
n=160
p=200
K=4 # number of groups 

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
library(Matrix)
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

compare_par=function(i, Sigmatrue, n, p){
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=PCAkmeans(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error_random=PCAkmeans_random(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error_decrease=PCAkmeans_decrease(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=BCDnormal2(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  
  result_error=c(error1,error_random,error_decrease,error2)
  result_error
}

x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue4, n, p)

write.csv(x,file=paste0("table2_", n,"_", p, ".csv"),quote=F,row.names = F)

##########################
### n=120 p=200        ###
##########################
n=120
p=200
K=4 # number of groups 

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
library(Matrix)
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

compare_par=function(i, Sigmatrue, n, p){
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=PCAkmeans(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error_random=PCAkmeans_random(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error_decrease=PCAkmeans_decrease(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=BCDnormal2(X,K,lambda1,lambda2,theta, Sigmatrue)$error
   
  result_error=c(error1,error_random,error_decrease,error2)
  result_error
}

x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue4, n, p)

write.csv(x,file=paste0("table2_", n,"_", p, ".csv"),quote=F,row.names = F)

##########################
### n=200 p=120        ###
##########################
n=200
p=120
K=4 # number of groups 

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
library(Matrix)
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

compare_par=function(i, Sigmatrue, n, p){
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=PCAkmeans(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error_random=PCAkmeans_random(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error_decrease=PCAkmeans_decrease(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=BCDnormal2(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  
  result_error=c(error1,error_random,error_decrease,error2)
  result_error
}

x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue4, n, p)

write.csv(x,file=paste0("table2_", n,"_", p, ".csv"),quote=F,row.names = F)

##########################
### n=200 p=160        ###
##########################
n=200
p=160
K=4 # number of groups 

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
library(Matrix)
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

compare_par=function(i, Sigmatrue, n, p){
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=PCAkmeans(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error_random=PCAkmeans_random(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error_decrease=PCAkmeans_decrease(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=BCDnormal2(X,K,lambda1,lambda2,theta, Sigmatrue)$error
   
  result_error=c(error1,error_random,error_decrease,error2)
  result_error
}

x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue4, n, p)

write.csv(x,file=paste0("table2_", n,"_", p, ".csv"),quote=F,row.names = F)

