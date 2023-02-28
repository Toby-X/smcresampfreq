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

uniform_spacing = function(N){
  a = -log(runif(N+1))
  b = cumsum(a)
  c = sort(b,decreasing = F)
  return(c/b[length(b)])
}

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

multinomial = function(w){
  return(inv_cdf(uniform_spacing(length(w)),w))
}

# stratified resampling
stratified = function(w){
  m = length(w)
  su = (runif(m)+0:(m-1))/m
  return(inv_cdf(su,w))
}

numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

## density
ldens = function(t,x){
  dens = 0.5*dmvnorm(x[1:t],3*rep(1,t),diag(rep(1,t)))+0.5*dmvnorm(x[1:t],-3*rep(1,t),diag(rep(1,t)))
  return(log(dens))
}

## SMC
ess <- function(l){
  set.seed(l)
  X = array(rep(0,n_particles*n_dim),c(n_particles,n_dim))
  rejuvs = 0
  x.estimate = rep(0,n_dim)
  X[,1] = rnorm(n_particles,0,3)
  W = exp(apply(X, 1, ldens,t=1)-log(dnorm(X[,1],0,3)))
  w = W/sum(W)
  idx = stratified(w)
  X = X[idx,]
  W = rep(1/n_particles,n_particles)
  rejuvs = rejuvs+1
  for (i in 2:(n_dim-1)) {
    X[,i] = rnorm(n_particles,0,3)
    W = exp(log(W)+apply(X, 1, ldens,t=i)-log(dnorm(X[,i],0,3)))
    idx = stratified(w)
    X = X[idx,]
    W = rep(1/n_particles,n_particles)
    rejuvs = rejuvs + 1
  }
  i = n_dim
  X[,i] = rnorm(n_particles,0,3)
  W = exp(log(W)+apply(X, 1, ldens,t=i)-log(dnorm(X[,i],0,3)))
  w = W/sum(W)
  x.estimate = as.vector(w%*%X)
  return(list(estimate = x.estimate,rejuvs=rejuvs))
}

seq_dens <- function(X,current){
  if (current == 1){
    res = 0.5*dnorm(X[,current+1],3,1)+0.5*dnorm(X[,current+1],-3,1)
  } else {
    res = 0.5*dnorm(X[,current+1],3,1)+0.5*dnorm(X[,current+1],-3,1)+
      0.5*dnorm(X[,current-1],3,1)+0.5*dnorm(X[,current-1],-3,1)
  }
}

new <- function(l){
  set.seed(l)
  X = array(rep(0,n_particles*n_dim),c(n_particles,n_dim))
  rejuvs = 0
  x.estimate = rep(0,n_dim)
  cri = rep(0,n_dim)
  
  X[,1] = rnorm(n_particles,0,3)
  W = exp(apply(X, 1, ldens,t=1)-log(dnorm(X[,1],0,3)))
  w = W/sum(W)
  # Resampling
  X[,2] = rnorm(n_particles,0,3)
  W.tmp = exp(apply(X,1,ldens,t=2))-log(dnorm(X[,2],0,3))
  cri[1] = sum(w*seq_dens(X,1)*(W/mean(W)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)))
  if (cri[1]>0){
    idx = stratified(w)
    X = X[idx,]
    W = rep(1/n_particles,n_particles)
    rejuvs = rejuvs+1
  }
  
  for (i in 2:(n_dim-1)) {
    X[,i] = rnorm(n_particles,0,3)
    W1 = exp(apply(X, 1, ldens,t=i)-log(dnorm(X[,i],0,3)))
    W = W*W1
    w = W/sum(W)
    #resampling
    X[,i+1] = rnorm(n_particles,0,3)
    W.tmp = exp(apply(X,1,ldens,t=i+1))-log(dnorm(X[,i+1],0,3))
    cri[i] = sum(w*seq_dens(X,1)*(W1/mean(W1)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)-1))# the same
    if (cri[i]>0){
      idx = stratified(w)
      X = X[idx,]
      W = rep(1/n_particles,n_particles)
      rejuvs = rejuvs+1
    }
  }
  i = n_dim
  X[,i] = rnorm(n_particles,0,3)
  W = exp(log(W)+apply(X, 1, ldens,t=i)-log(dnorm(X[,i],0,3)))
  w = W/sum(W)
  x.estimate = as.vector(w%*%X)
  return(list(estimate = x.estimate,rejuvs=rejuvs))
}

res1 = foreach (l = 1:m,.combine = rbind,.packages = "mvtnorm") %dopar% {
  return(ess(l))
}

res1 = foreach (l = 1:m,.combine = rbind,.packages = "mvtnorm") %dopar% {
  return(new(l))
}

# load("GMGM.RData")
# X = array(rep(0,m*n_dim*length(threshold)),c(m,n_dim,length(threshold)))
# rejuvs = matrix(rep(0,m*length(threshold)),nrow=m)
# for (i in 1:m){
#   X[i,,] = t(res[i,]$estimate)
#   rejuvs[i,] = res[i,]$rejuvs
# }
# colMeans(rejuvs)
# plot(X[,19,1],X[,20,1])
# mse = rep(0,length(threshold))
# for (j in 1:length(threshold)) {
#   mse[j] = sum(X[,,j]^2)
# }
# mse = mse/m
# plot(ess,mse,xlab="ESS Threshold",ylab="MSE")

# mse = matrix(rep(0,m*length(threshold)),nrow = m)
# llikelihood = matrix(rep(0,m*length(threshold)),nrow = m)
# rejuvs = matrix(rep(0,m*length(threshold)),nrow=m)
# for (i in 1:m){
#   mse[i,] = res[i,]$mse
#   llikelihood[i,] = res[i,]$llikelihood
#   rejuvs[i,] = res[i,]$rejuvs
# }
# mse.mean = colMeans(mse)
# llikelihood.mean = colMeans(llikelihood)
# rejuvs.mean = colMeans(rejuvs)
# plot(mse.mean)
# plot(llikelihood.mean)

save.image("/public1/home/scf0347/ResampFreq/GMG/GMGNew.RData")
