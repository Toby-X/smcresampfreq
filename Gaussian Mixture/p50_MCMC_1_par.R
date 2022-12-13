#-*- coding:utf-8 -*-
# change into .1 variance proposal for every random walk
path1 = "/public1/home/scf0347/ResampFreq/GaussianMixture/mixture.dat"
mixture.dat = read.table(path1,header=TRUE)
y = mixture.dat$y

## libraries
library(Boom)
library(parallel)
library(foreach)
library(doParallel)
numCores = detectCores()
registerDoParallel(numCores)

## Given Value
n = 500#number of particles
p = 50#number of distributions to sample from using smc
m = 200#number of repetitive experiments
threshold = seq(0.4,1,length=7)
mu = array(rep(0,n*p*2*length(threshold)),c(n,p,2,length(threshold),m))# the first layer is mu1, the second layer is mu2, the same is as follows
lambda = array(rep(0,n*p*2*length(threshold)),c(n,p,2,length(threshold),m))
omega = array(rep(0,n*p*2*length(threshold)*m),c(n,p,2,length(threshold),m))
wu = matrix(rep(1,n*m),c(n,m))
mse.mu = array(rep(0,m*length(threshold)),c(m,length(threshold)))
mse.lambda = array(rep(0,m*length(threshold)),c(m,length(threshold)))
mse.omega = array(rep(0,m*length(threshold)),c(m,length(threshold)))


# prior choices from Richardson and Green P735
kexi = mean(y)
R = max(y)-min(y)
K = 1/R^2
set.seed(101)
beta = rgamma(1,0.2,10/R^2)
alpha = 2
delta = 1

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
  if (n<=10){
    return(n*0.01/10)
  }else{
    r = (1/0.01)^(1/40)
    return(0.01*r^(n-10))
  }
}

log.fn <- function(n,mu,lambda,omega){
  res = log.likelihood(mu,lambda,omega)*phi(n)+sum(log(dnorm(mu,kexi,K^(-1/2))))+
    sum(log(dgamma(lambda,alpha,beta)))+sum(log(ddirichlet(omega,c(delta,delta))))
  return(res)
}

## 这里lognormal multiplicative还要再推导一下，现在是推导出来的结果
acceptance <- function(n,mu,lambda,omega,mut,lambdat,omegat){
  u = exp(omega[1])/(1+exp(omega[1]))
  ut = exp(omegat[1])/(1+exp(omegat[1]))
  log.res = log.fn(n,mu,lambda,omega)-log.fn(n,mut,lambdat,omegat)+
    log(omegat[1]*(1-omegat[1])/(omega[1]*(1-omega[1])))+
    sum(log(lambda)-log(lambdat))
  res = exp(log.res)
  return(min(res,1))
}

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

## Using suboptimal backward kernels
## extract the first sample

# change for parallel settings
w_u = matrix(rep(1,n*m),nrow=n)
w_n = matrix(rep(1,n*m),nrow=n)
idx = matrix(rep(0,n*m),nrow=n)
w_u.tmp = matrix(rep(1,n*m),nrow=n)

mu.update = array(rep(0,n*2*m),c(n,2,m))
mu.tmp = mu.update
lambda.update = array(rep(0,n*2*m),c(n,2,m))
lambda.tmp = lambda.update
u.update = array(rep(0,n*m),c(n,m))
u.tmp = u.update
omega.tmp = array(rep(0,n*m*2),c(n,2,m))

acrate.mu = rep(0,m)
acrate.lambda1 = rep(0,m)
acrate.lambda2 = rep(0,m)
acrate.omega = rep(0,m)
rand = rep(0,m)
divident = rep(0,m)

idx1 = array(rep(0,n*m),c(n,m))
idx2 = array(rep(0,n*m),c(n,m))
mu10 = array(rep(0,n*m),c(n,m))
mu7 = array(rep(0,n*m),c(n,m))
omega10 = array(rep(0,n*m),c(n,m))
omega7 = array(rep(0,n*m),c(n,m))
lambda_all = array(rep(0,2*n*m),c(2*n,m))

foreach (l = 1:m) %do% {
for (j in 1:length(threshold)) {
set.seed(l)

for (i in 1:n){
  mu[i,1,,j,l] = rnorm(2,kexi,K^(-1/2))
  lambda[i,1,,j,l] = rgamma(2,alpha,beta)
  omega[i,1,,j,l] = rdirichlet(1,c(delta,delta))
}

## estimate weights for step 1
for (i in 1:n) {
  w_u[i,l] = exp(log.fn(1,mu[i,1,,j,l],lambda[i,1,,j,l],omega[i,1,,j,l]))/prod(dnorm(mu[i,1,,j,l],kexi,K^(-1/2)))/prod(dgamma(lambda[i,1,,j,l],alpha,beta))/ddirichlet(omega[i,1,,j,l],c(delta,delta))
}
w_n[,l] = w_u[,l]/sum(w_u[,l])
if (1/sum(w_n^2) < n){
  idx[,l] = sample(1:n,n,replace=T,prob=w_n)
  lambda[,,,,l] = lambda[idx,,,,l]
  mu[,,,,l] = mu[idx,,,,l]
  omega[,,,,l] = omega[idx,,,,l]
  w_n[,l] = rep(1/n,n)
}

## MAIN LOOP
for (i in 2:p){
  if (any(is.na(w_n[,l]))){
    lambda[,,1,j,l] = 9999
    lambda[,,2,j,l] = 0
    omega[,,1,j,l] = 1
    omega[,,2,j,l] = 0
    mu[,,1,j,l] = 9999
    mu[,,2,j,l] = 0
    break
  }
  ## MCMC Move
  mu.update[,,l] = array(rnorm(n*2,0,.1),c(n,2))
  mu.tmp[,,l] = mu[,i-1,,j,l]+mu.update[,,l]
  lambda.update[,,l] = array(rlnorm(n*2,0,.1),c(n,2))
  lambda.tmp[,,l] = lambda[,i-1,,j,l]*lambda.update[,,l]
  u.update[,l] = rnorm(n,0,.1)
  u.tmp[,l] = log(omega[,i-1,1,j,l]/(1-omega[,i-1,1,j,l]))+u.update[,l]
  omega.tmp[,,l] = omega[,i-1,,j,l]
  omega.tmp[,1,l] = exp(u.tmp[,l])/(1+exp(u.tmp[,l]))
  omega.tmp[,2,l] = 1-omega.tmp[,1,l]
  # 加个循环把每个样本都MCMC update了
  ## 这里由于三个一起update所以效益可能会低，要不要查查它update一共多少次
  for (k in 1:n) {
    acrate.mu[l] = acceptance.mu(i,mu.tmp[k,,l],lambda.tmp[k,,l],omega.tmp[k,,l],mu[k,i-1,,j,l],lambda[k,i-1,,j,l],omega[k,i-1,,j,l])
    acrate.lambda1[l] = acceptance.lambda1(i,mu.tmp[k,,l],lambda.tmp[k,,l],omega.tmp[k,,l],mu[k,i-1,,j,l],lambda[k,i-1,,j,l],omega[k,i-1,,j,l])
    acrate.lambda2[l] = acceptance.lambda2(i,mu.tmp[k,,l],lambda.tmp[k,,l],omega.tmp[k,,l],mu[k,i-1,,j,l],lambda[k,i-1,,j,l],omega[k,i-1,,j,l])
    acrate.omega[l] = acceptance.omega(i,mu.tmp[k,,l],lambda.tmp[k,,l],omega.tmp[k,,l],mu[k,i-1,,j,l],lambda[k,i-1,,j,l],omega[k,i-1,,j,l])
    rand[l] = runif(1)
    if (rand[l]<=acrate.mu[l]){
      mu[k,i,,j,l] = mu.tmp[k,,l]
    }else{
      mu[k,i,,j,l] = mu[k,i-1,,j,l]
    }
    if (rand[l]<=acrate.lambda1[l]){
      lambda[k,i,1,j,l] = lambda.tmp[k,1,l]
    }else{
      lambda[k,i,1,j,l] = lambda[k,i-1,1,j,l]
    }
    if (rand[l]<=acrate.lambda2[l]){
      lambda[k,i,2,j,l] = lambda.tmp[k,2,l]
    }else{
      lambda[k,i,2,j,l] = lambda[k,i-1,2,j,l]
    }
    if (rand[l]<=acrate.omega[l]){
      omega[k,i,,j,l] = omega.tmp[k,,l]
    }else{
      omega[k,i,,j,l] = omega[k,i-1,,j,l]
    }
    for (o in 1:n) {
      divident[l] = divident[l] + exp(log.fn(i-1,mu[o,i-1,,j,l],lambda[o,i-1,,j,l],omega[o,i-1,,j,l]))*Kn(mu[k,i,,j,l],lambda[k,i,,j,l],omega[k,i,,j,l],mu[o,i-1,,j,l],lambda[o,i-1,,j,l],omega[o,i-1,,j,l])
    }
    w_u.tmp[k,l] = exp(log.fn(i,mu[k,i,,j,l],lambda[k,i,,j,l],omega[k,i,,j,l]))/divident[l]
  }
  w_u[,l] = w_n[,l]*w_u.tmp[,l]
  w_n[,l] = w_u[,l]/sum(w_u[,l])
  ## Resampling
  if (any(is.na(w_n[,l]))){
    lambda[,,1,j,l] = 9999
    lambda[,,2,j,l] = 0
    omega[,,1,j,l] = 1
    omega[,,2,j,l] = 0
    mu[,,1,j,l] = 9999
    mu[,,2,j,l] = 0
    break
  }
  if (1/sum(w_n[,l]^2)<threshold[j]){
    idx[,l] = sample(1:n,n,replace=T,prob=w_n[,l])
    lambda[,,,,l] = lambda[idx,,,,l]
    mu[,,,,l] = mu[idx,,,,l]
    omega[,,,,l] = omega[idx,,,,l]
    w_n[,l] = rep(1/n,n)
  }
}


idx1[,l] = mu[,50,1,j,l]>8
idx2[,l] = mu[,50,2,j,l]>8
mu10[,l] = c(mu[idx1,50,1,j,l],mu[idx2,50,2,j,l])
mu7[,l] = c(mu[!idx1,50,1,j,l],mu[!idx2,50,2,j,l])

omega10[,l] = c(omega[idx1,50,1,j,l],omega[idx2,50,2,j,l])
omega7[,l] = c(omega[!idx1,50,1,j,l],omega[!idx2,50,2,j,l])

lambda_all[,l] = c(lambda[,50,1,j,l],lambda[,50,2,j,l])

mse.mu[l,j] = mse.mu[l,j] + (mean(mu7[,l])-7)^2+var(mu7[,l])+(mean(mu10[,l])-10)^2+var(mu10[,l])
mse.omega[l,j] = mse.omega[l,j] +  (mean(omega10[,l])-0.3)^2+var(omega10[,l])
mse.lambda[l,j] = mse.lambda[l,j] + (mean(lambda_all[,l])-4)^2+var(lambda_all[,l])
}
}

mse.mu = mse.mu/n
mse.lambda = mse.lambda/n
mse.omega = mse.omega/n

#mse.mu = mse.mu/m
#mse.lambda = mse.lambda/m
#mse.omega = mse.omega/m

save.image("/public1/home/scf0347/ResampFreq/GaussianMixture/p50MCMC1par.RData")