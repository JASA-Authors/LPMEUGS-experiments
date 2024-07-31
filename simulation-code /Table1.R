library(MASS)
library(pdfCluster) # adj.rand.index
library(Matrix)

set.seed(1234)

###############
### Model 1 ###
###############
n=100
p=120
K=4 # number of groups 
sizeblock=p/K # size of blocks

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue1))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=160
K=4 # number of groups 
sizeblock=p/K # size of blocks

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue1))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=200
K=4 # number of groups 
sizeblock=p/K # size of blocks

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue1))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=120
p=100
K=4 # number of groups 
sizeblock=p/K # size of blocks

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue1))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=160
p=100
K=4 # number of groups 
sizeblock=p/K # size of blocks

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue1))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=200
p=100
K=4 # number of groups 
sizeblock=p/K # size of blocks

# Model 1
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue1))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)


###############
### Model 2 ###
###############
n=100
p=120
K=4 # number of groups 
sizeblock=p/K # size of blocks
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue2=kronecker(diag(K),ar1_cor(p/K,0.9)) 

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue2))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=160
K=4 # number of groups 
sizeblock=p/K # size of blocks
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue2=kronecker(diag(K),ar1_cor(p/K,0.9)) 

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue2))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=200
K=4 # number of groups 
sizeblock=p/K # size of blocks
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue2=kronecker(diag(K),ar1_cor(p/K,0.9)) 

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue2))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=200
p=120
K=4 # number of groups 
sizeblock=p/K # size of blocks
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue2=kronecker(diag(K),ar1_cor(p/K,0.9)) 

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue2))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=200
p=160
K=4 # number of groups 
sizeblock=p/K # size of blocks
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue2=kronecker(diag(K),ar1_cor(p/K,0.9)) 

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue2))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=200
p=160
K=4 # number of groups 
sizeblock=p/K # size of blocks
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue2=kronecker(diag(K),ar1_cor(p/K,0.9)) 

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue2))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

###############
### Model 3 ###
###############
n=100
p=120
K=4 # number of groups 
sizeblock=p/K # size of blocks
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  a_new=sample(1:p,replace = FALSE)
  Sigmatrue3=Sigmatrue1[a_new,]
  Sigmatrue3=Sigmatrue3[,a_new]
  
  truelable=c()
  for (i in 1:K) {
    truelable=c(truelable,replicate(sizeblock,i))
  }
  
  truelable=truelable[a_new]
  
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue3))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=160
K=4 # number of groups 
sizeblock=p/K # size of blocks
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  a_new=sample(1:p,replace = FALSE)
  Sigmatrue3=Sigmatrue1[a_new,]
  Sigmatrue3=Sigmatrue3[,a_new]
  
  truelable=c()
  for (i in 1:K) {
    truelable=c(truelable,replicate(sizeblock,i))
  }
  
  truelable=truelable[a_new]
  
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue3))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=200
K=4 # number of groups 
sizeblock=p/K # size of blocks
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  a_new=sample(1:p,replace = FALSE)
  Sigmatrue3=Sigmatrue1[a_new,]
  Sigmatrue3=Sigmatrue3[,a_new]
  
  truelable=c()
  for (i in 1:K) {
    truelable=c(truelable,replicate(sizeblock,i))
  }
  
  truelable=truelable[a_new]
  
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue3))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=200
p=120
K=4 # number of groups 
sizeblock=p/K # size of blocks
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  a_new=sample(1:p,replace = FALSE)
  Sigmatrue3=Sigmatrue1[a_new,]
  Sigmatrue3=Sigmatrue3[,a_new]
  
  truelable=c()
  for (i in 1:K) {
    truelable=c(truelable,replicate(sizeblock,i))
  }
  
  truelable=truelable[a_new]
  
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue3))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=200
p=160
K=4 # number of groups 
sizeblock=p/K # size of blocks
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  a_new=sample(1:p,replace = FALSE)
  Sigmatrue3=Sigmatrue1[a_new,]
  Sigmatrue3=Sigmatrue3[,a_new]
  
  truelable=c()
  for (i in 1:K) {
    truelable=c(truelable,replicate(sizeblock,i))
  }
  
  truelable=truelable[a_new]
  
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue3))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=200
K=4 # number of groups 
sizeblock=p/K # size of blocks
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock) # true covariance matrix

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  a_new=sample(1:p,replace = FALSE)
  Sigmatrue3=Sigmatrue1[a_new,]
  Sigmatrue3=Sigmatrue3[,a_new]
  
  truelable=c()
  for (i in 1:K) {
    truelable=c(truelable,replicate(sizeblock,i))
  }
  
  truelable=truelable[a_new]
  
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue3))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

###############
### Model 4 ###
###############
n=100
p=120
K=4 # number of groups 
CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

truelable=c(replicate(0.1*p,1),replicate(0.2*p,2),replicate(0.3*p,3),replicate(0.4*p,4))

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue4))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=160
K=4 # number of groups 
CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

truelable=c(replicate(0.1*p,1),replicate(0.2*p,2),replicate(0.3*p,3),replicate(0.4*p,4))

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue4))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=200
K=4 # number of groups 
CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

truelable=c(replicate(0.1*p,1),replicate(0.2*p,2),replicate(0.3*p,3),replicate(0.4*p,4))

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue4))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=120
p=100
K=4 # number of groups 
CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

truelable=c(replicate(0.1*p,1),replicate(0.2*p,2),replicate(0.3*p,3),replicate(0.4*p,4))

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue4))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)


n=160
p=100
K=4 # number of groups 
CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

truelable=c(replicate(0.1*p,1),replicate(0.2*p,2),replicate(0.3*p,3),replicate(0.4*p,4))

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue4))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=200
p=100
K=4 # number of groups 
CSblock1=matrix(0.9,nrow = 0.2*p,ncol = 0.2*p)+diag(0.2*p)*0.1
CSblock2=matrix(0.9,nrow = 0.4*p,ncol = 0.4*p)+diag(0.4*p)*0.1
ar1_cor <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}
Sigmatrue4=bdiag(ar1_cor(0.1*p,0.9),CSblock1,ar1_cor(0.3*p,0.9),CSblock2)

truelable=c(replicate(0.1*p,1),replicate(0.2*p,2),replicate(0.3*p,3),replicate(0.4*p,4))

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue4))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

###############
### Model 5 ###
###############
n=100
p=120
K=4 # number of groups 
sizeblock=p/K
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock)

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  Sigmatrue5=Sigmatrue1
  for (t in 1:(p-1)) {
    for (j in (t+1):p) {
      t=rbinom(1,1,0.02)*0.01
      Sigmatrue5[t,j]=Sigmatrue5[t,j]+r
      Sigmatrue5[j,t]=Sigmatrue5[j,t]+r
    }
  }
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue5))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=160
K=4 # number of groups 
sizeblock=p/K
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock)

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  Sigmatrue5=Sigmatrue1
  for (t in 1:(p-1)) {
    for (j in (t+1):p) {
      t=rbinom(1,1,0.02)*0.01
      Sigmatrue5[t,j]=Sigmatrue5[t,j]+r
      Sigmatrue5[j,t]=Sigmatrue5[j,t]+r
    }
  }
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue5))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=100
p=200
K=4 # number of groups 
sizeblock=p/K
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock)

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  Sigmatrue5=Sigmatrue1
  for (t in 1:(p-1)) {
    for (j in (t+1):p) {
      t=rbinom(1,1,0.02)*0.01
      Sigmatrue5[t,j]=Sigmatrue5[t,j]+r
      Sigmatrue5[j,t]=Sigmatrue5[j,t]+r
    }
  }
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue5))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=120
p=100
K=4 # number of groups 
sizeblock=p/K
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock)

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  Sigmatrue5=Sigmatrue1
  for (t in 1:(p-1)) {
    for (j in (t+1):p) {
      t=rbinom(1,1,0.02)*0.01
      Sigmatrue5[t,j]=Sigmatrue5[t,j]+r
      Sigmatrue5[j,t]=Sigmatrue5[j,t]+r
    }
  }
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue5))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=160
p=100
K=4 # number of groups 
sizeblock=p/K
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock)

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  Sigmatrue5=Sigmatrue1
  for (t in 1:(p-1)) {
    for (j in (t+1):p) {
      t=rbinom(1,1,0.02)*0.01
      Sigmatrue5[t,j]=Sigmatrue5[t,j]+r
      Sigmatrue5[j,t]=Sigmatrue5[j,t]+r
    }
  }
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue5))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)

n=200
p=100
K=4 # number of groups 
sizeblock=p/K
CSblock=matrix(0.9,nrow = sizeblock,ncol = sizeblock)+diag(sizeblock)*0.1
Sigmatrue1=kronecker(diag(K),CSblock)

truelable=c()
for (i in 1:K) {
  truelable=c(truelable,replicate(sizeblock,i))
}

repli=200
error1=c()
error2=c()

for (i in 1:repli) {
  Sigmatrue5=Sigmatrue1
  for (t in 1:(p-1)) {
    for (j in (t+1):p) {
      t=rbinom(1,1,0.02)*0.01
      Sigmatrue5[t,j]=Sigmatrue5[t,j]+r
      Sigmatrue5[j,t]=Sigmatrue5[j,t]+r
    }
  }
  X=as.matrix(mvrnorm(n=n,mu=rep(0,p),Sigma = Sigmatrue5))
  # each row of X minuses its mean
  Y=sweep(X,1,rowMeans(X))
  
  # PCA
  eiv=eigen(t(Y)%*%Y/p)$vectors
  V=as.matrix(eiv[,1:(K-1)],ncol=K-1)
  
  # Kmeans
  clusterinf1=kmeans(V, K)$cluster
  clusterinf2=kmeans(t(X), K)$cluster
  
  error1=c(error1,adj.rand.index(truelable, clusterinf1))
  error2=c(error2,adj.rand.index(truelable, clusterinf2))
}

mean(error1)
mean(error2)
