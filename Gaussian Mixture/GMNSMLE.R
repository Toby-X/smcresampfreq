#-*- coding:utf-8 -*-
# change into .1 variance proposal for every random walk
# path1 = file.choose()
path1 = "/public1/home/scf0347/ResampFreq/GaussianMixture/mixture.dat"
mixture.dat = read.table(path1,header=TRUE)
y = mixture.dat$y
# set.seed(109)
# y = rep(0,400)
# for (i in 1:500) {
#   if(runif(1)<0.7){
#     y[i] = rnorm(1,7,0.5)
#   }else{
#     y[i] = rnorm(1,10,0.5)
#   }
# }

## libraries
library(Boom)
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

## Given Value
n = 200#number of particles
p = 50#number of distributions to sample from using smc
m = 96#number of repetitive experiments

## doSNOW progress bar
pb = txtProgressBar(max = m,style=3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress=progress)

# prior choices from Richardson and Green P735
set.seed(5001)
kexi = mean(y)
R = max(y)-min(y)
K = 1/R^2
set.seed(101)
beta = rgamma(1,0.2,10/R^2)
alpha = 2
delta = 1

uniform_spacing = function(N){
  a = -log(runif(N+1))
  b = cumsum(a)
  return(b[-length(b)]/b[length(b)])
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

stratified = function(w){
  m = length(w)
  su = (runif(m)+0:(m-1))/m
  return(inv_cdf(su,w))
}

##the data is generated from .7*dnorm(x,7,.5) + .3*dnorm(x,10,.5)
## pin is proportional to likelihood(y;theta_r)^phi_n * f(theta_r)
## update mu via additive normal random-walk proposal
## lambda via multiplicative log-normal random-walk
## omega via additive normal random-walk
log.likelihood <- function(mu,lambda,omega){
  sum(log(omega[1]*dnorm(y,mu[1],lambda^(-1/2))+omega[2]*dnorm(y,mu[2],lambda^(-1/2))))
}

## using one iteration of MH kernel
## using first 50 steps of annealing
## using dnorm(0,1) as proposal in initial state and in MCMC
phi <- function(n){
  if (n<=0.2*p){
    return(n*0.15/(0.2*p))
  }else if(n<=0.6*p){
    return(n*(0.4-0.15)/(0.4*p)+0.15)
  }else{
    return(n*(1-0.4)/(0.4*p)+0.4)
  }
}

log.fn <- function(n,mu,lambda,omega){
  res = log.likelihood(mu,lambda,omega)*phi(n)+sum(log(dnorm(mu,kexi,K^(-1/2))))+
    sum(log(dgamma(lambda,alpha,beta)))+sum(log(ddirichlet(omega,c(delta,delta))))
  return(res)
}

## 这里lognormal multiplicative还要再推导一下，现在是推导出来的结果
# acceptance <- function(n,mu,lambda,omega,mut,lambdat,omegat){
#   u = exp(omega[1])/(1+exp(omega[1]))
#   ut = exp(omegat[1])/(1+exp(omegat[1]))
#   log.res = log.fn(n,mu,lambda,omega)-log.fn(n,mut,lambdat,omegat)+
#     log(omegat[1]*(1-omegat[1])/(omega[1]*(1-omega[1])))+
#     sum(log(lambda)-log(lambdat))
#   res = exp(log.res)
#   return(min(res,1))
# }

acceptance.mu <- function(n,mu,lambda,omega,mut,lambdat,omegat){
  log.res = log.fn(n,mu,lambda,omega)-log.fn(n,mut,lambdat,omegat)
  res = exp(log.res)
  return(min(res,1))
}

acceptance.lambda1 <- function(n,mu,lambda,omega,mut,lambdat,omegat){
  log.res = log.fn(n,mu,lambda,omega)-log.fn(n,mut,lambdat,omegat)+log(lambda[1])-log(lambdat[1])
  res = exp(log.res)
  return(min(res,1))
}

acceptance.lambda2 <- function(n,mu,lambda,omega,mut,lambdat,omegat){
  log.res = log.fn(n,mu,lambda,omega)-log.fn(n,mut,lambdat,omegat)+log(lambda[2])-log(lambdat[2])
  res = exp(log.res)
  return(min(res,1))
}

acceptance.omega <- function(n,mu,lambda,omega,mut,lambdat,omegat){
  u = exp(omega[1])/(1+exp(omega[1]))
  ut = exp(omegat[1])/(1+exp(omegat[1]))
  log.res = log.fn(n,mu,lambda,omega)-log.fn(n,mut,lambdat,omegat)+
    log(omegat[1]*(1-omegat[1])/(omega[1]*(1-omega[1])))
  res = exp(log.res)
  return(min(res,1))
}

## Next step is to write the whole process including weight
Kn <- function(mu,lambda,omega,mut,lambdat,omegat){
  u = exp(omega[1])/(1+exp(omega[1]))
  ut = exp(omegat[1])/(1+exp(omegat[1]))
  res = sum(log(dnorm(mu-mut,0,.1)))+sum(log(dlnorm(lambda/lambdat,0,.1)))+
    log(dnorm(u-ut,0,.1))-log(omega[1])-log(1-omega[1])
  res = exp(res)
  return(res)
}

gmm <- function(){
  mu7.estimate = 0
  mu10.estimate = 0
  omega7.estimate = 0
  lambda.estimate = 0
  resamp = rep(0,p)
  mu = array(rep(0,n*p*2),c(n,p,2))# the first layer is mu1, the second layer is mu2, the same is as follows
  lambda = array(rep(0,n*p*2),c(n,p,2))
  omega = array(rep(0,n*p*2),c(n,p,2))
  rejuvs = 0
  w_u = rep(1,n)
  w_u.tmp = rep(1,n)
  
  for (i in 1:n){
    mu[i,1,] = rnorm(2,kexi,K^(-1/2))
    lambda[i,1,] = rgamma(2,alpha,beta)
    omega[i,1,] = rdirichlet(1,c(delta,delta))
  }
  
  ## estimate weights for step 1
  for (i in 1:n) {
    w_u[i] = exp(log.fn(1,mu[i,1,],lambda[i,1,],omega[i,1,]))/prod(dnorm(mu[i,1,],kexi,K^(-1/2)))/prod(dgamma(lambda[i,1,],alpha,beta))/ddirichlet(omega[i,1,],c(delta,delta))
  }
  if (any(is.na(w_u))){
    idx.na = is.na(w_u)
    w_u[idx.na]=0
  }
  w_n = w_u/sum(w_u)
  if (1/sum(w_n^2)<n){
    idx = stratified(w_n)
    lambda = lambda[idx,,]
    mu = mu[idx,,]
    omega = omega[idx,,]
    w_n = rep(1/n,n)
    rejuvs = rejuvs + 1
    resamp[1] = 1
  }
  
  ## MAIN LOOP 
  for (i in 2:(p-1)){
    ## MCMC Move
    mu.update = array(rnorm(n*2,0,.1),c(n,2))
    mu.tmp = mu[,i-1,]+mu.update
    lambda.update = array(rlnorm(n*2,0,.1),c(n,2))
    lambda.tmp = lambda[,i-1,]*lambda.update
    u.update = rnorm(n,0,.1)
    u.tmp = log(omega[,i-1,1]/(1-omega[,i-1,1]))+u.update
    omega.tmp = omega[,i-1,]
    omega.tmp[,1] = exp(u.tmp)/(1+exp(u.tmp))
    omega.tmp[,2] = 1-omega.tmp[,1]
    # 加个循环把每个样本都MCMC update了
    ## 这里由于三个一起update所以效益可能会低，要不要查查它update一共多少次
    for (k in 1:n) {
      acrate.mu = acceptance.mu(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
      acrate.lambda1 = acceptance.lambda1(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
      acrate.lambda2 = acceptance.lambda2(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
      acrate.omega = acceptance.omega(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
      rand = runif(1)
      if (is.na(acrate.mu)){
        mu[k,i,] = mu.tmp[k,]
      }else if (rand<=acrate.mu){
        mu[k,i,] = mu.tmp[k,]
      }else{
        mu[k,i,] = mu[k,i-1,]
      }
      if (is.na(acrate.lambda1)){
        lambda[k,i,1] = lambda.tmp[k,1]
      }else if (rand<=acrate.lambda1){
        lambda[k,i,1] = lambda.tmp[k,1]
      }else{
        lambda[k,i,1] = lambda[k,i-1,1]
      }
      if (is.na(acrate.lambda2)){
        lambda[k,i,2] = lambda.tmp[k,2]
      }else if (rand<=acrate.lambda2){
        lambda[k,i,2] = lambda.tmp[k,2]
      }else{
        lambda[k,i,2] = lambda[k,i-1,2]
      }
      if (is.na(acrate.omega)){
        omega[k,i,] = omega.tmp[k,]
      }else if (rand<=acrate.omega){
        omega[k,i,] = omega.tmp[k,]
      }else{
        omega[k,i,] = omega[k,i-1,]
      }
      divident=0
      for (o in 1:n) {
        divident = divident + exp(log.fn(i-1,mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,]))*Kn(mu[k,i,],lambda[k,i,],omega[k,i,],mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,])
      }
      w_u.tmp[k] = exp(log.fn(i,mu[k,i,],lambda[k,i,],omega[k,i,]))/divident
    }
    w_u = w_n*w_u.tmp
    if (any(is.na(w_u))){
      idx.na = is.na(w_u)
      w_u[idx.na]=0
    }
    w_n = w_u/sum(w_u)
    
    if (1/sum(w_n^2)<n){
      idx = stratified(w_n)
      lambda = lambda[idx,,]
      mu = mu[idx,,]
      omega = omega[idx,,]
      w_n = rep(1/n,n)
      rejuvs = rejuvs + 1
      resamp[i] = 1
    }
  }
  i = p
  ## MCMC Move
  mu.update = array(rnorm(n*2,0,.1),c(n,2))
  mu.tmp = mu[,i-1,]+mu.update
  lambda.update = array(rlnorm(n*2,0,.1),c(n,2))
  lambda.tmp = lambda[,i-1,]*lambda.update
  u.update = rnorm(n,0,.1)
  u.tmp = log(omega[,i-1,1]/(1-omega[,i-1,1]))+u.update
  omega.tmp = omega[,i-1,]
  omega.tmp[,1] = exp(u.tmp)/(1+exp(u.tmp))
  omega.tmp[,2] = 1-omega.tmp[,1]
  # 加个循环把每个样本都MCMC update了
  ## 这里由于三个一起update所以效益可能会低，要不要查查它update一共多少次
  for (k in 1:n) {
    acrate.mu = acceptance.mu(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
    acrate.lambda1 = acceptance.lambda1(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
    acrate.lambda2 = acceptance.lambda2(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
    acrate.omega = acceptance.omega(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
    rand = runif(1)
    if (is.na(acrate.mu)){
      mu[k,i,] = mu.tmp[k,]
    }else if (rand<=acrate.mu){
      mu[k,i,] = mu.tmp[k,]
    }else{
      mu[k,i,] = mu[k,i-1,]
    }
    if (is.na(acrate.lambda1)){
      lambda[k,i,1] = lambda.tmp[k,1]
    }else if (rand<=acrate.lambda1){
      lambda[k,i,1] = lambda.tmp[k,1]
    }else{
      lambda[k,i,1] = lambda[k,i-1,1]
    }
    if (is.na(acrate.lambda2)){
      lambda[k,i,2] = lambda.tmp[k,2]
    }else if (rand<=acrate.lambda2){
      lambda[k,i,2] = lambda.tmp[k,2]
    }else{
      lambda[k,i,2] = lambda[k,i-1,2]
    }
    if (is.na(acrate.omega)){
      omega[k,i,] = omega.tmp[k,]
    }else if (rand<=acrate.omega){
      omega[k,i,] = omega.tmp[k,]
    }else{
      omega[k,i,] = omega[k,i-1,]
    }
    divident=0
    for (o in 1:n) {
      divident = divident + exp(log.fn(i-1,mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,]))*Kn(mu[k,i,],lambda[k,i,],omega[k,i,],mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,])
    }
    w_u.tmp[k] = exp(log.fn(i,mu[k,i,],lambda[k,i,],omega[k,i,]))/divident
  }
  w_u = w_n*w_u.tmp
  w_n = w_u/sum(w_u)
  
  idx1 = mu[,p,1]>8
  idx2 = mu[,p,2]>8
  # !idx2 == idx1 always holds in this situation
  weight = c(w_n[idx1],w_n[idx2])
  mu10 = c(mu[idx1,p,1],mu[idx2,p,2])
  mu7 = c(mu[!idx2,p,2],mu[!idx1,p,1])
  
  omega10 = c(omega[idx1,p,1],omega[idx2,p,2])
  omega7 = c(omega[!idx2,p,2],omega[!idx1,p,1])
  
  lambda10 = c(lambda[idx1,p,1],lambda[idx2,p,2])
  lambda7 = c(lambda[!idx2,p,2],lambda[!idx1,p,1])
  
  mu7.estimate = sum(weight*mu7)
  mu10.estimate = sum(weight*mu10)
  omega7.estimate = sum(weight*omega7)
  lambda.estimate = (sum(weight*lambda7)+sum(weight*lambda10))/2
  
  return(list(mu7=mu7.estimate,mu10=mu10.estimate,lambda=lambda.estimate,omega7=omega7.estimate, rejuvs = rejuvs, time=resamp))
}

seq_dens <- function(mu,omega,lambda,current){
  res = rep(0,n)
  if (current == 1){
    for (i in 1:n) {
      res[i] = log.fn(current+1,mu[i,current+1,],lambda[i,current+1,],omega[i,current+1,])
    }
  } else{
    for (i in 1:n) {
      res[i] = log.fn(current-1,mu[i,current-1,],lambda[i,current-1,],omega[i,current-1,])+log.fn(current+1,mu[i,current+1,],lambda[i,current+1,],omega[i,current+1,])
    }
  }
  return(exp(res))
}

newm <- function(){
  mu7.estimate = 0
  mu10.estimate = 0
  omega7.estimate = 0
  lambda.estimate = 0
  mu = array(rep(0,n*p*2),c(n,p,2))# the first layer is mu1, the second layer is mu2, the same is as follows
  lambda = array(rep(0,n*p*2),c(n,p,2))
  omega = array(rep(0,n*p*2),c(n,p,2))
  rejuvs = 0
  w_u = rep(1,n)
  w_u.tmp = rep(1,n)
  w_u.tmp1 = rep(1,n)
  cri = rep(0,p)
  resamp = rep(0,p)
  for (i in 1:n){
    mu[i,1,] = rnorm(2,kexi,K^(-1/2))
    lambda[i,1,] = rgamma(2,alpha,beta)
    omega[i,1,] = rdirichlet(1,c(delta,delta))
  }
  
  ## estimate weights for step 1
  for (i in 1:n) {
    w_u[i] = exp(log.fn(1,mu[i,1,],lambda[i,1,],omega[i,1,]))/prod(dnorm(mu[i,1,],kexi,K^(-1/2)))/prod(dgamma(lambda[i,1,],alpha,beta))/ddirichlet(omega[i,1,],c(delta,delta))
  }
  if (any(is.na(w_u))){
    idx.na = is.na(w_u)
    w_u[idx.na]=0
  }
  w_n = w_u/sum(w_u)
  
  ## RESAMPLING
  i = 2
  mu[,i,] = mu[,i-1,]
  lambda[,i,] = lambda[,i-1,]
  omega[,i,] = omega[,i-1,]
  # 加个循环把每个样本都MCMC update了
  ## 这里由于三个一起update所以效益可能会低，要不要查查它update一共多少次
  for (k in 1:n) {
    divident=0
    for (o in 1:n) {
      divident = divident + exp(log.fn(i-1,mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,]))*Kn(mu[k,i,],lambda[k,i,],omega[k,i,],mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,])
    }
    w_u.tmp[k] = exp(log.fn(i,mu[k,i,],lambda[k,i,],omega[k,i,]))/divident
  }
  if (any(is.na(w_u.tmp))){
    idx.na = is.na(w_u.tmp)
    w_u.tmp[idx.na] = 0
  }
  cri[1] = sum(w_n*seq_dens(mu,omega,lambda,1)*(w_u/mean(w_u)*w_u.tmp/mean(w_u.tmp)-w_u.tmp/mean(w_u.tmp)))
  if (!is.na(cri[1]) & cri[1]>0){
    idx = stratified(w_n)
    lambda = lambda[idx,,]
    mu = mu[idx,,]
    omega = omega[idx,,]
    w_n = rep(1/n,n)
    rejuvs = rejuvs + 1
    resamp[1] = 1
  }
  
  ## MAIN LOOP
  for (i in 2:(p-1)){
    ## MCMC Move
    mu.update = array(rnorm(n*2,0,.1),c(n,2))
    mu.tmp = mu[,i-1,]+mu.update
    lambda.update = array(rlnorm(n*2,0,.1),c(n,2))
    lambda.tmp = lambda[,i-1,]*lambda.update
    u.update = rnorm(n,0,.1)
    u.tmp = log(omega[,i-1,1]/(1-omega[,i-1,1]))+u.update
    omega.tmp = omega[,i-1,]
    omega.tmp[,1] = exp(u.tmp)/(1+exp(u.tmp))
    omega.tmp[,2] = 1-omega.tmp[,1]
    # 加个循环把每个样本都MCMC update了
    ## 这里由于三个一起update所以效益可能会低，要不要查查它update一共多少次
    for (k in 1:n) {
      acrate.mu = acceptance.mu(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
      acrate.lambda1 = acceptance.lambda1(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
      acrate.lambda2 = acceptance.lambda2(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
      acrate.omega = acceptance.omega(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
      rand = runif(1)
      if (is.na(acrate.mu)){
        mu[k,i,] = mu.tmp[k,]
      }else if (rand<=acrate.mu){
        mu[k,i,] = mu.tmp[k,]
      }else{
        mu[k,i,] = mu[k,i-1,]
      }
      if (is.na(acrate.lambda1)){
        lambda[k,i,1] = lambda.tmp[k,1]
      }else if (rand<=acrate.lambda1){
        lambda[k,i,1] = lambda.tmp[k,1]
      }else{
        lambda[k,i,1] = lambda[k,i-1,1]
      }
      if (is.na(acrate.lambda2)){
        lambda[k,i,2] = lambda.tmp[k,2]
      }else if (rand<=acrate.lambda2){
        lambda[k,i,2] = lambda.tmp[k,2]
      }else{
        lambda[k,i,2] = lambda[k,i-1,2]
      }
      if (is.na(acrate.omega)){
        omega[k,i,] = omega.tmp[k,]
      }else if (rand<=acrate.omega){
        omega[k,i,] = omega.tmp[k,]
      }else{
        omega[k,i,] = omega[k,i-1,]
      }
      divident=0
      for (o in 1:n) {
        divident = divident + exp(log.fn(i-1,mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,]))*Kn(mu[k,i,],lambda[k,i,],omega[k,i,],mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,])
      }
      w_u.tmp[k] = exp(log.fn(i,mu[k,i,],lambda[k,i,],omega[k,i,]))/divident
    }
    w_u = w_n*w_u.tmp
    if (any(is.na(w_u))){
      idx.na = is.na(w_u)
      w_u[idx.na]=0
    }
    w_n = w_u/sum(w_u)
    
    ## RESAMPLING
    i = i+1
    mu[,i,] = mu[,i-1,]
    lambda[,i,] = lambda[,i-1,]
    omega[,i,] = omega[,i-1,]
    # 加个循环把每个样本都MCMC update了
    ## 这里由于三个一起update所以效益可能会低，要不要查查它update一共多少次
    for (k in 1:n) {
      divident=0
      for (o in 1:n) {
        divident = divident + exp(log.fn(i-1,mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,]))*Kn(mu[k,i,],lambda[k,i,],omega[k,i,],mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,])
      }
      w_u.tmp1[k] = exp(log.fn(i,mu[k,i,],lambda[k,i,],omega[k,i,]))/divident
    }
    if (any(is.na(w_u.tmp1))){
      idx.na = is.na(w_u.tmp1)
      w_u.tmp1[idx.na] = 0
    }
    if (any(is.na(w_u.tmp))){
      idx.na = is.na(w_u.tmp)
      w_u.tmp[idx.na] = 0
    }
    cri[i-1] = sum(w_n*seq_dens(mu,omega,lambda,1)*(w_u.tmp1/mean(w_u.tmp1)*w_u.tmp/mean(w_u.tmp)-w_u.tmp1/mean(w_u.tmp1)))
    if (!is.na(cri[i-1]) & cri[i-1]>0){
      idx = stratified(w_n)
      lambda = lambda[idx,,]
      mu = mu[idx,,]
      omega = omega[idx,,]
      w_n = rep(1/n,n)
      rejuvs = rejuvs + 1
      resamp[i-1]=1
    }
  }
  
  i = p
  mu.update = array(rnorm(n*2,0,.1),c(n,2))
  mu.tmp = mu[,i-1,]+mu.update
  lambda.update = array(rlnorm(n*2,0,.1),c(n,2))
  lambda.tmp = lambda[,i-1,]*lambda.update
  u.update = rnorm(n,0,.1)
  u.tmp = log(omega[,i-1,1]/(1-omega[,i-1,1]))+u.update
  omega.tmp = omega[,i-1,]
  omega.tmp[,1] = exp(u.tmp)/(1+exp(u.tmp))
  omega.tmp[,2] = 1-omega.tmp[,1]
  # 加个循环把每个样本都MCMC update了
  ## 这里由于三个一起update所以效益可能会低，要不要查查它update一共多少次
  for (k in 1:n) {
    acrate.mu = acceptance.mu(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
    acrate.lambda1 = acceptance.lambda1(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
    acrate.lambda2 = acceptance.lambda2(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
    acrate.omega = acceptance.omega(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,],lambda[k,i-1,],omega[k,i-1,])
    rand = runif(1)
    if (is.na(acrate.mu)){
      mu[k,i,] = mu.tmp[k,]
    }else if (rand<=acrate.mu){
      mu[k,i,] = mu.tmp[k,]
    }else{
      mu[k,i,] = mu[k,i-1,]
    }
    if (is.na(acrate.lambda1)){
      lambda[k,i,1] = lambda.tmp[k,1]
    }else if (rand<=acrate.lambda1){
      lambda[k,i,1] = lambda.tmp[k,1]
    }else{
      lambda[k,i,1] = lambda[k,i-1,1]
    }
    if (is.na(acrate.lambda2)){
      lambda[k,i,2] = lambda.tmp[k,2]
    }else if (rand<=acrate.lambda2){
      lambda[k,i,2] = lambda.tmp[k,2]
    }else{
      lambda[k,i,2] = lambda[k,i-1,2]
    }
    if (is.na(acrate.omega)){
      omega[k,i,] = omega.tmp[k,]
    }else if (rand<=acrate.omega){
      omega[k,i,] = omega.tmp[k,]
    }else{
      omega[k,i,] = omega[k,i-1,]
    }
    divident=0
    for (o in 1:n) {
      divident = divident + exp(log.fn(i-1,mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,]))*Kn(mu[k,i,],lambda[k,i,],omega[k,i,],mu[o,i-1,],lambda[o,i-1,],omega[o,i-1,])
    }
    w_u.tmp[k] = exp(log.fn(i,mu[k,i,],lambda[k,i,],omega[k,i,]))/divident
  }
  w_u = w_n*w_u.tmp
  if (any(is.na(w_u))){
    idx.na = is.na(w_u)
    w_u[idx.na]=0
  }
  w_n = w_u/sum(w_u)
  
  idx1 = mu[,p,1]>8
  idx2 = mu[,p,2]>8
  # !idx2 == idx1 always holds in this situation
  weight = c(w_n[idx1],w_n[idx2])
  mu10 = c(mu[idx1,p,1],mu[idx2,p,2])
  mu7 = c(mu[!idx2,p,2],mu[!idx1,p,1])
  
  omega10 = c(omega[idx1,p,1],omega[idx2,p,2])
  omega7 = c(omega[!idx2,p,2],omega[!idx1,p,1])
  
  lambda10 = c(lambda[idx1,p,1],lambda[idx2,p,2])
  lambda7 = c(lambda[!idx2,p,2],lambda[!idx1,p,1])
  
  mu7.estimate = sum(weight*mu7)
  mu10.estimate = sum(weight*mu10)
  omega7.estimate = sum(weight*omega7)
  lambda.estimate = (sum(weight*lambda7)+sum(weight*lambda10))/2
  
  return(list(mu7=mu7.estimate,mu10=mu10.estimate,lambda=lambda.estimate,omega7=omega7.estimate, rejuvs = rejuvs, criterion = cri, time = resamp))
}

ess = foreach (l = 1:m,.combine = rbind,.packages = "Boom",
               .options.snow=opts) %dopar% {
                 # if (l == 38+48 || l == 48+36){
                 #   set.seed(l^2)
                 # } else {
                 #   set.seed(l)
                 # }
                 set.seed(l)
                 return(gmm())
               }

new = foreach (l = 1:m,.combine = rbind,.packages = "Boom",
               .options.snow=opts) %dopar% {
                 set.seed(l)
                 return(newm())
               }
close(pb)

# load("p50MCMC1StrMLE.RData")
# head(ess)
# head(new)
# 
# ess.resamp = matrix(rep(0,p*m),nrow=m)
# new.resamp = matrix(rep(0,p*m),nrow=m)
# for (i in 1:m) {
#   ess.resamp[i,] = ess[i,]$time
#   new.resamp[i,] = new[i,]$time
# }
# 
# plot(1:50,colMeans(ess.resamp[nzero,]),"l",ylim=c(0,1))
# par(new = T)
# plot(1:50,colMeans(new.resamp[nzero,]),"l",col="red",ylim = c(0,1))
# 
# plot(colMeans(ess.resamp[nzero,])-colMeans(new.resamp[nzero,]))
# 
# diff.resamp = abs(ess.resamp-new.resamp)
# plot(diff.resamp)
# 
# new.rejuvs = rep(0,m)
# for (i in 1:m){
#   new.rejuvs[i] = new[i,]$rejuvs
# }
# 
# nzero = new.rejuvs>=30
# 
# ess.mu7 = rep(0,m)
# ess.mu10 = rep(0,m)
# ess.omega7 = rep(0,m)
# ess.lambda = rep(0,m)
# 
# for (i in 1:m){
#   ess.mu7[i] = ess[i,]$mu7
#   ess.mu10[i] = ess[i,]$mu10
#   ess.lambda[i] = ess[i,]$lambda
#   ess.omega7[i] = ess[i,]$omega7
# }
# 
# mse.mu = (ess.mu7-7)^2+(ess.mu10-10)^2
# mse.lambda = (ess.lambda-4)^2
# mse.omega = (ess.omega7-.7)^2
# 
# mean(mse.mu)
# mean(mse.lambda)
# mean(mse.omega)
# 
# new.mu7 = rep(0,m)
# new.mu10 = rep(0,m)
# new.omega7 = rep(0,m)
# new.lambda = rep(0,m)
# new.rejuvs = rep(0,m)
# for (i in 1:m){
#   new.mu7[i] = new[i,]$mu7
#   new.mu10[i] = new[i,]$mu10
#   new.lambda[i] = new[i,]$lambda
#   new.omega7[i] = new[i,]$omega7
#   new.rejuvs[i] = new[i,]$rejuvs
# }
# rejvs.mean = mean(new.rejuvs)
# mse.mu.n = (new.mu7[nzero]-7)^2+(new.mu10[nzero]-10)^2
# mse.lambda.n = (new.lambda[nzero]-4)^2
# mse.omega.n = (new.omega7[nzero]-.7)^2
# 
# mean(mse.lambda)
# mean(mse.lambda.n)
# mean(mse.mu)
# mean(mse.mu.n)
# mean(mse.omega)
# mean(mse.omega.n)
# 
# boxplot(mse.lambda.n)
# boxplot(mse.lambda)

save.image("/public1/home/scf0347/ResampFreq/GaussianMixture/p50MCMC1StrMLE.RData")