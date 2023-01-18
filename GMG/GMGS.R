#-*- coding:utf-8 -*-
# given values
n_dim = 20
n_particles = 500
m = 480
threshold = seq(0.1,1,by=0.1)

library(mvtnorm)
library(parallel)
library(foreach)
library(doSNOW)

numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

## density
ldens = function(t,x){
  dens = 0.5*dmvnorm(x[1:t],3*rep(1,t),diag(rep(1,t)))+0.5*dmvnorm(x[1:t],-3*rep(1,t),diag(rep(1,t)))
  return(log(dens))
}

# resampling scheme using systematic resampling
inv_cdf = function(su,w){
  j = 1
  s = w[1]
  m = length(su)
  a = rep(0,m)
  for (i in 1:m){
    while (su[i]>s) {
      j = j+1
      s = s+w[j]
    }
    a[i] = j
  }
  return(a)
}

systematic = function(w){
  m = length(w)
  su = (rep(runif(1),m)+0:(m-1))/m
  return(inv_cdf(su,w))
}

## SMC
doit <- function(l){
  set.seed(l)
  X = array(rep(0,n_particles*n_dim*length(threshold)),c(n_particles,n_dim,length(threshold)))
  rejuvs = rep(0,length(threshold))
  x.estimate = matrix(rep(0,n_dim*length(threshold)),nrow=length(threshold))
  for (j in 1:length(threshold)) {
    X[,1,j] = rnorm(n_particles,0,3)
    W = exp(apply(X[,,j], 1, ldens,t=1)-log(dnorm(X[,1,j],0,3)))
    w = W/sum(W)
    if (1/sum(w^2) < threshold[j]*n_particles){
      idx = systematic(w)
      X = X[idx,,]
      W = rep(1/n_particles,n_particles)
      rejuvs[j] = rejuvs[j]+1
    }
    for (i in 2:n_dim) {
      X[,i,j] = rnorm(n_particles,0,3)
      W = exp(log(W)+apply(X[,,j], 1, ldens,t=i)-log(dnorm(X[,i,j],0,3)))
      w = W/sum(W)
      if (1/sum(w^2) < threshold[j]*n_particles){
        idx = systematic(w)
        X = X[idx,,]
        W = rep(1/n_particles,n_particles)
        rejuvs[j] = rejuvs[j] + 1
      }
    }
    w = W/sum(W)
    x.estimate[j,] = as.vector(w%*%X[,,j])
  }
  return(list(estimate = x.estimate,rejuvs=rejuvs))
}

res = foreach (l = 1:m,.combine = rbind,.packages = "mvtnorm") %dopar% {
  return(doit(l))
}

load("GMGS.RData")
X = array(rep(0,m*n_dim*length(threshold)),c(m,n_dim,length(threshold)))
rejuvs = matrix(rep(0,m*length(threshold)),nrow=m)
for (i in 1:m){
  X[i,,] = t(res[i,]$estimate)
  rejuvs[i,] = res[i,]$rejuvs
}
colMeans(rejuvs)
plot(X[,19,1],X[,20,1])
mse = rep(0,length(threshold))
for (j in 1:length(threshold)) {
  mse[j] = sum(X[,,j]^2)
}
mse = mse/m
plot(ess,mse,xlab="ESS Threshold",ylab="MSE")

save.image("/public1/home/scf0347/ResampFreq/GMG/GMGS.RData")
