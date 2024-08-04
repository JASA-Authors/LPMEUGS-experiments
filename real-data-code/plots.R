set.seed(1234)

################
### Function ###
################
Group_Detect=function(X, K, lambda1, lambda2, theta){
  # OUR method
  # input:
  # X: data
  # K: number of groups
  # theta: threshold for stop
  # lambda1: parameter for glmnet
  # lambda2: parameter for glasso regression
  
  # output:
  # Omega: estimated precision matrix
  # order1: estimated group information
  
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
  
  
  return(list(Omega=newOmega,order1=Bcluster))
}

Group_Precision=function(X,K,lambda1,lambda2,theta){
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

############
### Data ###
############
data_breast=read.csv("D:/data_breast.csv",header = T)
data_label=read.csv("D:/label.csv",header = T)
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
library(ggplot2)
library(patchwork)
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

p1+p2

## Heat Map
library(reshape2)
a1=Group_Detect(X,2,0.01,0.01,0.01)
Omegahat=a1$Omega 
order1=unlist(a1$order1)
label=colnames(X)
label=label[order1]
data=expand.grid(X=label, Y=label)
data$Z=melt(Omegahat)[,3]
data$Z[which(data$Z>0.2)]=0.2
data$Z[which(data$Z< -0.2)]=-0.2
g1=ggplot(data)+ aes(X, Y, fill= Z) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",#colors in the scale
                       midpoint=0,breaks=c(-0.2,-0.1,0,0.1,0.2),
                       labels=c("< -0.2","-0.1","0","0.1",
                                ">0.2")) +
  ggtitle("OUR")+geom_tile()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
    axis.text.x=element_blank(),axis.ticks.x=element_blank(),  
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(legend.title=element_blank())

Omegahat2=flare::sugm(X,nlambda=1,method = "clime")$icov[[1]]
Omegahat2=Omegahat2[order1,]
Omegahat2=Omegahat2[,order1]
data2=expand.grid(X=label, Y=label)
data2$Z=melt(Omegahat2)[,3]
data2$Z[which(data2$Z>0.2)]=0.2
data2$Z[which(data2$Z< -0.2)]=-0.2
g2=ggplot(data2)+ aes(X, Y, fill= Z) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",#colors in the scale
                       midpoint=0,breaks=c(-0.2,-0.1,0,0.1,0.2),
                       labels=c("< -0.2","-0.1","0","0.1",
                                ">0.2")) +
  ggtitle("CLIME")+geom_tile()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),  
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(legend.title=element_blank(),legend.position = "none")

Omegahat3=flare::sugm(X,nlambda=1,method = "tiger")$icov[[1]]
Omegahat3=Omegahat3[order1,]
Omegahat3=Omegahat3[,order1]
data3=expand.grid(X=label, Y=label)
data3$Z=melt(Omegahat3)[,3]
data3$Z[which(data3$Z>0.2)]=0.2
data3$Z[which(data3$Z< -0.2)]=-0.2
g3=ggplot(data3)+ aes(X, Y, fill= Z) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",#colors in the scale
                       midpoint=0,breaks=c(-0.2,-0.1,0,0.1,0.2),
                       labels=c("< -0.2","-0.1","0","0.1",
                                ">0.2")) +
  ggtitle("TIGER")+geom_tile()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),  
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(legend.title=element_blank())

Omegahat4=glasso::glasso(cov(X), rho=.1)$wi
Omegahat4=Omegahat4[order1,]
Omegahat4=Omegahat4[,order1]
data4=expand.grid(X=label, Y=label)
data4$Z=melt(Omegahat4)[,3]
data4$Z[which(data4$Z>0.2)]=0.2
data4$Z[which(data4$Z< -0.2)]=-0.2
g4=ggplot(data4)+ aes(X, Y, fill= Z) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",#colors in the scale
                       midpoint=0,breaks=c(-0.2,-0.1,0,0.1,0.2),
                       labels=c("< -0.2","-0.1","0","0.1",
                                ">0.2")) +
  ggtitle("G-LASSO")+geom_tile()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),  
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(legend.title=element_blank())

Omegahat5=Group_Precision(X,ncol(X),0.01,0.01,0.1)$Omega
Omegahat5=Omegahat5[order1,]
Omegahat5=Omegahat5[,order1]
data5=expand.grid(X=label, Y=label)
data5$Z=melt(Omegahat5)[,3]
data5$Z[which(data5$Z>0.2)]=0.2
data5$Z[which(data5$Z< -0.2)]=-0.2
g5=ggplot(data5)+ aes(X, Y, fill= Z) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",#colors in the scale
                       midpoint=0,breaks=c(-0.2,-0.1,0,0.1,0.2),
                       labels=c("< -0.2","-0.1","0","0.1",
                                ">0.2")) +
  ggtitle("NO-GROUP")+geom_tile()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),  
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(legend.title=element_blank())

Omegahat6=solve(POET::POET(t(X),2)$SigmaY)
Omegahat6=Omegahat6[order1,]
Omegahat6=Omegahat6[,order1]
data6=expand.grid(X=label, Y=label)
data6$Z=melt(Omegahat6)[,3]
data6$Z[which(data6$Z>0.2)]=0.2
data6$Z[which(data6$Z< -0.2)]=-0.2
g6=ggplot(data6)+ aes(X, Y, fill= Z) + 
  scale_fill_gradient2(low="blue",mid="white", high="red",#colors in the scale
                       midpoint=0,breaks=c(-0.2,-0.1,0,0.1,0.2),
                       labels=c("< -0.2","-0.1","0","0.1",
                                ">0.2")) +
  ggtitle("POET")+geom_tile()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),  
        axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(legend.title=element_blank())

g1+g5+g4+g2+g3+g6+
  plot_layout(ncol=3, guides = "collect")
