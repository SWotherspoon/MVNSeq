## -----------------------------------------------------------------------------
set.seed(31)
library(MVNSeq)

## -----------------------------------------------------------------------------
y <- as.matrix(log(iris[,1:4]))

## -----------------------------------------------------------------------------
library(MVNSeq)
km <- kmeans(y,3)
fit <- mvnMix(y,km$cluster)

## -----------------------------------------------------------------------------
library(mclust)
fit.mclust <- Mclust(y,G=3,modelNames=c("VVV"))
summary(fit.mclust,parameters=T)

## -----------------------------------------------------------------------------
table(max.col(fit$pars$P),fit.mclust$classification)

## -----------------------------------------------------------------------------
j <- order(max.col(table(max.col(fit$pars$P),fit.mclust$classification)))
sapply(fit$pars$mvn[j],function(p) p$mu)
lapply(fit$pars$mvn[j],function(p) p$Sigma)

## -----------------------------------------------------------------------------
## Number of components K, responses q and groups m
library(mvtnorm)
sim <- function(K,q,m,uscale=0,vscale=1,minl=100,meanl=500,hmm=T) {

  ## Vector of sequence lengths
  ns <- minl+rpois(m,max(meanl-minl,0))
  n <- sum(ns)

  ## Group structure
  gr <- rep(1:m,ns)
  J <- c(1,cumsum(ns))

  ## Responses, group means, population means, and covariances
  Y <- array(0,c(n,q,K))
  A <- array(0,c(m,q,K))
  Mu <- array(0,c(K,q))
  U <- vector(mode="list",K)
  V <- vector(mode="list",K)

  ## Simulate a response for each state and group
  for(k in 1:K) {
    ## Covariance of random effects
    U[[k]] <- crossprod(matrix(rnorm(q*q,0,uscale),q,q))
    ## Covariance of response about group means
    V[[k]] <- crossprod(matrix(rnorm(q*q,0,vscale),q,q))
    ## Population mean
    Mu[k,] <- runif(q,0,10)
    ## Group means
    A[,,k] <- rmvnorm(m,Mu[k,],U[[k]])
    ## Responses
    Y[,,k] <- A[gr,,k] + rmvnorm(n,rep(0,q),V[[k]])
  }

  ## state prior
  p <- runif(K)
  p <- p/sum(p)


  if(hmm) {
    ## Random transition matrix
    Q <- matrix(runif(K*K),K,K)
    Q <- Q/rowSums(Q)
  } else {
    Q <- matrix(p,K,K,byrow=T)
  }

  ## Simulate (Markovian) sequence of class membership and observations
  cl <- integer(n)
  y <- array(0,c(n,q))
  for(i in 1:m) {
    j <- J[i]
    cl[j] <- sample(K,1,prob=p)
    y[j,] <- Y[j,,cl[j]]
    for(j in (J[i]+1):J[i+1]) {
      cl[j] <- sample(K,1,prob=Q[cl[j-1],])
      y[j,] <- Y[j,,cl[j]]
    }
  }
  list(y=y,cl=cl,gr=gr,A=A,Mu=Mu,U=U,V=V,p=p,Q=Q)
}

## -----------------------------------------------------------------------------
s <- sim(K=3,q=3,m=1,hmm=F)

## -----------------------------------------------------------------------------
km <- kmeans(s$y,3)
fit <- mvnMix(s$y,km$cluster)

## -----------------------------------------------------------------------------
library(mclust)
fit.mclust <- Mclust(s$y,G=3,modelNames=c("VVV"))
summary(fit.mclust,parameters=T)

## -----------------------------------------------------------------------------
table(max.col(fit$pars$P),fit.mclust$classification)

## -----------------------------------------------------------------------------
j <- order(max.col(table(max.col(fit$pars$P),fit.mclust$classification)))
sapply(fit$pars$mvn[j],function(p) p$mu)
lapply(fit$pars$mvn[j],function(p) p$Sigma)

## -----------------------------------------------------------------------------
table(max.col(fit$pars$P),s$cl)

## -----------------------------------------------------------------------------
j <- order(max.col(table(max.col(fit$pars$P),s$cl)))
t(s$Mu)
sapply(fit$pars$mvn[j],function(p) p$mu)

## -----------------------------------------------------------------------------
s$V
lapply(fit$pars$mvn[j],function(p) p$Sigma)

## -----------------------------------------------------------------------------
km <- kmeans(s$y,s$Mu)
fit <- gmvnHMM(s$y,km$cluster,s$gr,common.transition = TRUE)

## -----------------------------------------------------------------------------
table(unlist(lapply(fit$pars,function(p) max.col(p$P))),s$cl)

## -----------------------------------------------------------------------------
## Check parameter estimates match simulation
opar <- par(mfrow=c(2,2),mar=c(1,1,1,1)+0.1)
matplot(s$A[,,1],t(sapply(fit$pars,function(p) p$mvn[[1]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,2],t(sapply(fit$pars,function(p) p$mvn[[2]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,3],t(sapply(fit$pars,function(p) p$mvn[[3]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
par(opar)

## -----------------------------------------------------------------------------
fit$pars[[1]]$mvn[[1]]$Sigma
s$V[[1]]
fit$pars[[1]]$mvn[[2]]$Sigma
s$V[[2]]
fit$pars[[1]]$mvn[[3]]$Sigma
s$V[[3]]

## -----------------------------------------------------------------------------
fit$pars[[1]]$Q
s$Q

## -----------------------------------------------------------------------------
km <- kmeans(s$y,3)
fit <- mvnMix(s$y,km$cluster)

## -----------------------------------------------------------------------------
table(max.col(fit$pars$P),s$cl)

## -----------------------------------------------------------------------------
j <- order(max.col(table(max.col(fit$pars$P),s$cl)))
t(s$Mu)
sapply(fit$pars$mvn[j],function(p) p$mu)

## -----------------------------------------------------------------------------
s$V
lapply(fit$pars$mvn[j],function(p) p$Sigma)

## -----------------------------------------------------------------------------
km <- kmeans(s$y,3)
fit <- mvnHMM(s$y,km$cluster)

## -----------------------------------------------------------------------------
table(max.col(fit$pars$P),s$cl)

## -----------------------------------------------------------------------------
j <- order(max.col(table(max.col(fit$pars$P),s$cl)))
t(s$Mu)
sapply(fit$pars$mvn[j],function(p) p$mu)

## -----------------------------------------------------------------------------
s$V
lapply(fit$pars$mvn[j],function(p) p$Sigma)

## -----------------------------------------------------------------------------
fit$Q[j,j]
s$Q

## -----------------------------------------------------------------------------
km <- kmeans(s$y,s$Mu)
fit <- grmvnMix(s$y,km$cluster,s$gr,common.fractions = TRUE)

## -----------------------------------------------------------------------------
table(unlist(lapply(fit$pars,function(p) max.col(p$P))),s$cl)

## -----------------------------------------------------------------------------
## Check parameter estimates match simulation
opar <- par(mfrow=c(2,2),mar=c(1,1,1,1)+0.1)
matplot(s$A[,,1],t(sapply(fit$pars,function(p) p$mvn[[1]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,2],t(sapply(fit$pars,function(p) p$mvn[[2]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,3],t(sapply(fit$pars,function(p) p$mvn[[3]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
par(opar)

## -----------------------------------------------------------------------------
fit$pars[[1]]$mvn[[1]]$Sigma
s$V[[1]]
fit$pars[[1]]$mvn[[2]]$Sigma
s$V[[2]]
fit$pars[[1]]$mvn[[3]]$Sigma
s$V[[3]]

## -----------------------------------------------------------------------------
fit$muv[[1]]$U
s$U[[1]]
fit$muv[[2]]$U
s$U[[2]]
fit$muv[[3]]$U
s$U[[3]]

## -----------------------------------------------------------------------------
s <- sim(K=3,q=3,m=100,minl=100,meanl=500,hmm=T,uscale = 0.8)
pairs(s$y,pch=".",col=s$cl)

## -----------------------------------------------------------------------------
km <- kmeans(s$y,s$Mu)
fit <- gmvnMix(s$y,km$cluster,s$gr,common.fractions = TRUE)

## -----------------------------------------------------------------------------
table(unlist(lapply(fit$pars,function(p) max.col(p$P))),s$cl)

## -----------------------------------------------------------------------------
## Check parameter estimates match simulation
opar <- par(mfrow=c(2,2),mar=c(1,1,1,1)+0.1)
matplot(s$A[,,1],t(sapply(fit$pars,function(p) p$mvn[[1]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,2],t(sapply(fit$pars,function(p) p$mvn[[2]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,3],t(sapply(fit$pars,function(p) p$mvn[[3]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
par(opar)

## -----------------------------------------------------------------------------
fit$pars[[1]]$mvn[[1]]$Sigma
s$V[[1]]
fit$pars[[1]]$mvn[[2]]$Sigma
s$V[[2]]
fit$pars[[1]]$mvn[[3]]$Sigma
s$V[[3]]

## -----------------------------------------------------------------------------
km <- kmeans(s$y,s$Mu)
fit <- gmvnHMM(s$y,km$cluster,s$gr,common.transition = TRUE)

## -----------------------------------------------------------------------------
table(unlist(lapply(fit$pars,function(p) max.col(p$P))),s$cl)

## -----------------------------------------------------------------------------
## Check parameter estimates match simulation
opar <- par(mfrow=c(2,2),mar=c(1,1,1,1)+0.1)
matplot(s$A[,,1],t(sapply(fit$pars,function(p) p$mvn[[1]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,2],t(sapply(fit$pars,function(p) p$mvn[[2]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,3],t(sapply(fit$pars,function(p) p$mvn[[3]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
par(opar)

## -----------------------------------------------------------------------------
fit$pars[[1]]$mvn[[1]]$Sigma
s$V[[1]]
fit$pars[[1]]$mvn[[2]]$Sigma
s$V[[2]]
fit$pars[[1]]$mvn[[3]]$Sigma
s$V[[3]]

## -----------------------------------------------------------------------------
fit$pars[[1]]$Q
s$Q

## -----------------------------------------------------------------------------
km <- kmeans(s$y,s$Mu)
fit <- grmvnHMM(s$y,km$cluster,s$gr,common.transition = TRUE)

## -----------------------------------------------------------------------------
table(unlist(lapply(fit$pars,function(p) max.col(p$P))),s$cl)

## -----------------------------------------------------------------------------
## Check parameter estimates match simulation
opar <- par(mfrow=c(2,2),mar=c(1,1,1,1)+0.1)
matplot(s$A[,,1],t(sapply(fit$pars,function(p) p$mvn[[1]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,2],t(sapply(fit$pars,function(p) p$mvn[[2]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,3],t(sapply(fit$pars,function(p) p$mvn[[3]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
par(opar)

## -----------------------------------------------------------------------------
fit$pars[[1]]$mvn[[1]]$Sigma
s$V[[1]]
fit$pars[[1]]$mvn[[2]]$Sigma
s$V[[2]]
fit$pars[[1]]$mvn[[3]]$Sigma
s$V[[3]]

## -----------------------------------------------------------------------------
fit$muv[[1]]$U
s$U[[1]]
fit$muv[[2]]$U
s$U[[2]]
fit$muv[[3]]$U
s$U[[3]]

## -----------------------------------------------------------------------------
fit$pars[[1]]$Q
s$Q

## -----------------------------------------------------------------------------
km <- kmeans(s$y,s$Mu)
fit <- grmvnMix(s$y,km$cluster,s$gr,common.fractions = TRUE)

## -----------------------------------------------------------------------------
table(unlist(lapply(fit$pars,function(p) max.col(p$P))),s$cl)

## -----------------------------------------------------------------------------
## Check parameter estimates match simulation
opar <- par(mfrow=c(2,2),mar=c(1,1,1,1)+0.1)
matplot(s$A[,,1],t(sapply(fit$pars,function(p) p$mvn[[1]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,2],t(sapply(fit$pars,function(p) p$mvn[[2]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,3],t(sapply(fit$pars,function(p) p$mvn[[3]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
par(opar)

## -----------------------------------------------------------------------------
fit$pars[[1]]$mvn[[1]]$Sigma
s$V[[1]]
fit$pars[[1]]$mvn[[2]]$Sigma
s$V[[2]]
fit$pars[[1]]$mvn[[3]]$Sigma
s$V[[3]]

## -----------------------------------------------------------------------------
fit$muv[[1]]$U
s$U[[1]]
fit$muv[[2]]$U
s$U[[2]]
fit$muv[[3]]$U
s$U[[3]]

## -----------------------------------------------------------------------------
km <- kmeans(s$y,s$Mu)
fit <- grmvnHMM(s$y,km$cluster,s$gr,common.transition = TRUE)

## -----------------------------------------------------------------------------
table(unlist(lapply(fit$pars,function(p) max.col(p$P))),s$cl)

## -----------------------------------------------------------------------------
## Check parameter estimates match simulation
opar <- par(mfrow=c(2,2),mar=c(1,1,1,1)+0.1)
matplot(s$A[,,1],t(sapply(fit$pars,function(p) p$mvn[[1]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,2],t(sapply(fit$pars,function(p) p$mvn[[2]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
matplot(s$A[,,3],t(sapply(fit$pars,function(p) p$mvn[[3]]$mu)),pch=16,cex=0.5,xlab="",ylab="")
par(opar)

## -----------------------------------------------------------------------------
fit$pars[[1]]$mvn[[1]]$Sigma
s$V[[1]]
fit$pars[[1]]$mvn[[2]]$Sigma
s$V[[2]]
fit$pars[[1]]$mvn[[3]]$Sigma
s$V[[3]]

## -----------------------------------------------------------------------------
fit$muv[[1]]$U
s$U[[1]]
fit$muv[[2]]$U
s$U[[2]]
fit$muv[[3]]$U
s$U[[3]]

## -----------------------------------------------------------------------------
fit$pars[[1]]$Q
s$Q

