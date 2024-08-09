library(MASS)
library (glmnet)
library(glasso)
library(expm)
library(flare) # CLIME, TIGER
library(POET)
library(foreach)
library(doParallel)

# Specify the number of cores to be used for parallel computing
clnum<-24
cl <- makeCluster(getOption("cl.cores", clnum))
registerDoParallel(cl)

################
### Function ###
################
Group_Detect=function(X, K, lambda1, lambda2, theta, Sigmatrue){
  # OUR method
  # input:
  # X: data
  # K: number of groups
  # theta: threshold for stop
  # lambda1: parameter for glmnet
  # lambda2: parameter for glasso regression
  
  # output:
  # Omega: estimated precision matrix
  # error: estimation error
  
  n=nrow(X)
  p=ncol(X)
  
  # centralize each row
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf=kmeans(V, K)$cluster
  
  # detect orders based on cluster size
  B_count=c()
  for (i in 1:K) {
    B_count=c(B_count,length(which(clusterinf==i)))
  }
  seq1=order(B_count,decreasing = FALSE)
  
  # reorder sample X
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

Group_Precision=function(X,K,lambda1,lambda2,theta, Sigmatrue){
  # ORACLE method
  # when K=P, it becomes NO-GROUP method 
  # input:
  # X: data
  # K: number of groups
  # theta: judge when to stop
  # lambda1: parameter for glmnet
  # lambda2: parameter for glasso regression
  
  # output:
  # Omega: estimated precision matrix
  # error: estimation error
  
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
  serror=sqrt(sum((solve(Sigmatrue)-newOmega)^2))
  
  return(list(Omega=newOmega,error=serror))
}

##########################
### n=160 p=200        ###
##########################
n=160
p=200
K=4 # number of groups 
sizeblock=p/K # size of blocks

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

compare_par=function(i, Sigmatrue, n, p){
  set.seed(i)
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=Group_Detect(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=Group_Precision(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error3=Group_Precision(X,p,lambda1,lambda2,theta, Sigmatrue)$error
  error5=sqrt(sum((solve(Sigmatrue)-glasso::glasso(cov(X), rho=.1)$wi)^2))
  error6=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "tiger")$icov[[1]])^2))
  error7=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "clime")$icov[[1]])^2))
  error8=sqrt(sum((solve(Sigmatrue)-solve(POET::POET(t(X))$SigmaY))^2))
  
  result_error=c(error1,error2,error3,error5,error6,error7,error8)
  result_error
}

x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue1, n, p)
colnames(x)=c("OUR","ORACLE","NO-GROUP","G-LASSO","TIGER","CLIME","POET")
x=rbind(x,apply(x,2,mean))
x=rbind(x,apply(x,2,sd))

write.csv(x,file=paste0("model1_", n,"_", p, ".csv"),quote=F,row.names = F)

##########################
### n=120 p=200        ###
##########################
n=120
p=200
K=4 # number of groups 
sizeblock=p/K # size of blocks

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

compare_par=function(i, Sigmatrue, n, p){
  set.seed(i)
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=Group_Detect(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=Group_Precision(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error3=Group_Precision(X,p,lambda1,lambda2,theta, Sigmatrue)$error
  error5=sqrt(sum((solve(Sigmatrue)-glasso::glasso(cov(X), rho=.1)$wi)^2))
  error6=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "tiger")$icov[[1]])^2))
  error7=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "clime")$icov[[1]])^2))
  error8=sqrt(sum((solve(Sigmatrue)-solve(POET::POET(t(X))$SigmaY))^2))
  
  result_error=c(error1,error2,error3,error5,error6,error7,error8)
  result_error
}

x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue1, n, p)
colnames(x)=c("OUR","ORACLE","NO-GROUP","G-LASSO","TIGER","CLIME","POET")
x=rbind(x,apply(x,2,mean))
x=rbind(x,apply(x,2,sd))

write.csv(x,file=paste0("model1_", n,"_", p, ".csv"),quote=F,row.names = F)

##########################
### n=200 p=120        ###
##########################
n=200
p=120
K=4 # number of groups 
sizeblock=p/K # size of blocks

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

compare_par=function(i, Sigmatrue, n, p){
  set.seed(i)
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=Group_Detect(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=Group_Precision(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error3=Group_Precision(X,p,lambda1,lambda2,theta, Sigmatrue)$error
  error4=sqrt(sum((solve(Sigmatrue)-solve(cov(X)))^2))
  error5=sqrt(sum((solve(Sigmatrue)-glasso::glasso(cov(X), rho=.1)$wi)^2))
  error6=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "tiger")$icov[[1]])^2))
  error7=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "clime")$icov[[1]])^2))
  error8=sqrt(sum((solve(Sigmatrue)-solve(POET::POET(t(X))$SigmaY))^2))
  
  result_error=c(error1,error2,error3,error4,error5,error6,error7,error8)
  result_error
}

x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue1, n, p)
colnames(x)=c("OUR","ORACLE","NO-GROUP","SAMPLE","G-LASSO","TIGER","CLIME","POET")
x=rbind(x,apply(x,2,mean))
x=rbind(x,apply(x,2,sd))

write.csv(x,file=paste0("model1_", n,"_", p, ".csv"),quote=F,row.names = F)

##########################
### n=200 p=160        ###
##########################
n=200
p=160
K=4 # number of groups 
sizeblock=p/K # size of blocks

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

compare_par=function(i, Sigmatrue, n, p){
  set.seed(i)
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=Group_Detect(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=Group_Precision(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error3=Group_Precision(X,p,lambda1,lambda2,theta, Sigmatrue)$error
  error4=sqrt(sum((solve(Sigmatrue)-solve(cov(X)))^2))
  error5=sqrt(sum((solve(Sigmatrue)-glasso::glasso(cov(X), rho=.1)$wi)^2))
  error6=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "tiger")$icov[[1]])^2))
  error7=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "clime")$icov[[1]])^2))
  error8=sqrt(sum((solve(Sigmatrue)-solve(POET::POET(t(X))$SigmaY))^2))
  
  result_error=c(error1,error2,error3,error4,error5,error6,error7,error8)
  result_error
}

x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue1, n, p)
colnames(x)=c("OUR","ORACLE","NO-GROUP","SAMPLE","G-LASSO","TIGER","CLIME","POET")
x=rbind(x,apply(x,2,mean))
x=rbind(x,apply(x,2,sd))

write.csv(x,file=paste0("model1_", n,"_", p, ".csv"),quote=F,row.names = F)


##########################
### n=50 p=200        ###
##########################
n=50
p=200
K=4 # number of groups 
sizeblock=p/K # size of blocks

lambda1=0.01 
lambda2=0.01 
theta=0.01 # judge when to stop the loop

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

compare_par=function(i, Sigmatrue, n, p){
  set.seed(i)
  X=as.matrix(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue))
  
  error1=Group_Detect(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error2=Group_Precision(X,K,lambda1,lambda2,theta, Sigmatrue)$error
  error3=Group_Precision(X,p,lambda1,lambda2,theta, Sigmatrue)$error
  error5=sqrt(sum((solve(Sigmatrue)-glasso::glasso(cov(X), rho=.1)$wi)^2))
  error6=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "tiger")$icov[[1]])^2))
  error7=sqrt(sum((solve(Sigmatrue)-flare::sugm(X,nlambda=1,method = "clime")$icov[[1]])^2))
  error8=sqrt(sum((solve(Sigmatrue)-solve(POET::POET(t(X))$SigmaY))^2))
  
  result_error=c(error1,error2,error3,error5,error6,error7,error8)
  result_error
}

# The environment and setup of the computer cluster may have an impact on the random number generating process of R.
x <- foreach(i=1:200,.combine='rbind') %dopar% compare_par(i, Sigmatrue1, n, p)
colnames(x)=c("OUR","ORACLE","NO-GROUP","G-LASSO","TIGER","CLIME","POET")
x=rbind(x,apply(x,2,mean))
x=rbind(x,apply(x,2,sd))

write.csv(x,file=paste0("model1_", n,"_", p, ".csv"),quote=F,row.names = F)

stopCluster(cl)
