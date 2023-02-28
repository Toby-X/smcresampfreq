#-*- coding:utf-8 -*-
# change into .25 variance proposal for every random walk
path1 = "/public1/home/scf0347/ResampFreq/GaussianMixture/mixture.dat"
mixture.dat = read.table(path1,header=TRUE)
y = mixture.dat$y

## libraries
library(Boom)

## Given Value
n = 500
p = 50
mu = array(rep(0,n*p*2),c(n,p,2))# the first layer is mu1, the second layer is mu2, the same is as follows
lambda = array(rep(0,n*p*2),c(n,p,2))
omega = array(rep(0,n*p*2),c(n,p,2))
wu = rep(0,n)
threshold = n

# prior choices from Richardson and Green P735
kexi = mean(y)
R = max(y)-min(y)
K = 1/R^2
set.seed(101)
beta = rgamma(1,0.2,10/R^2)
alpha = 2
delta = 1

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



for (i in 1:n){
  mu[i,1,] = rnorm(2,kexi,K^(-1/2))
  lambda[i,1,] = rgamma(2,alpha,beta)
  omega[i,1,] = rdirichlet(1,c(delta,delta))
}

## estimate weights for step 1
w_u = rep(1,n)
for (i in 1:n) {
  w_u[i] = exp(log.fn(1,mu[i,1,],lambda[i,1,],omega[i,1,]))/prod(dnorm(mu[i,1,],kexi,K^(-1/2)))/prod(dgamma(lambda[i,1,],alpha,beta))/ddirichlet(omega[i,1,],c(delta,delta))
}
if (any(is.na(w_u))){
  idx.na = is.na(w_u)
  w_u[idx.na]=0
}
w_n = w_u/sum(w_u)
if (1/sum(w_n^2) < threshold){
  idx = systematic(w_n)
  lambda = lambda[idx,,]
  mu = mu[idx,,]
  omega = omega[idx,,]
  w_n = rep(1/n,n)
}

##monitor accepted times
ac.mu = array(rep(0,n*p),c(n,p))
ac.lambda = array(rep(0,n*p*2),c(n,p,2))
ac.omega = array(rep(0,n*p),c(n,p))
##这里是不是实际上应该lambda1和lambda2分别进行Update呢(mu也是类似的)
##mu的acceptance rate对于1和2是一样的
##但是lambda的是不一样的

## MAIN LOOP
star_time <- Sys.time()
w_u.tmp = rep(1,n)
for (i in 2:p){
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
    divident = 0
    for (j in 1:n) {
      divident = divident + exp(log.fn(i-1,mu[j,i-1,],lambda[j,i-1,],omega[j,i-1,]))*Kn(mu[k,i,],lambda[k,i,],omega[k,i,],mu[j,i-1,],lambda[j,i-1,],omega[j,i-1,])
    }
    w_u.tmp[k] = exp(log.fn(i,mu[k,i,],lambda[k,i,],omega[k,i,]))/divident
  }
  w_u = w_n*w_u.tmp
  if (any(is.na(w_u))){
    idx.na = is.na(w_u)
    w_u[idx.na]=0
  }
  w_n = w_u/sum(w_u)
  ## Resampling
  if (1/sum(w_n^2)<threshold){
    idx = systematic(w_n)
    lambda = lambda[idx,,]
    mu = mu[idx,,]
    omega = omega[idx,,]
    w_n = rep(1/n,n)
  }
  sprintf("已完成 %d%%", round(i*100/p))
}

end_time <- Sys.time()
total_time = end_time-star_time

idx1 = mu[,50,1]>8
idx2 = mu[,50,2]>8
mu10 = c(mu[idx1,50,1],mu[idx2,50,2])
mu7 = c(mu[!idx1,50,1],mu[!idx2,50,2])
# (mean(mu7)-7)^2+var(mu7)
#str(mu7)
#str(mu10)
#hist(mu7)
#mean(mu7)
# hist(mu10)
#mean(mu7)
#mean(mu10)
##稍微好一点点，但是还是应该看看acceptance rate再做结论
#mean(rowSums(ac.mu))

omega10 = c(omega[idx1,50,1],omega[idx2,50,2])
omega7 = c(omega[!idx1,50,1],omega[!idx2,50,2])
# hist(omega10)
# hist(omega7)
#omega_all = c(omega10,omega7)
#hist(omega_all)
#mean(omega10)
#mean(omega7)
#mean(rowSums(ac.omega))

lambda_all = c(lambda[,50,1],lambda[,50,2])
# hist(lambda_all)
# mean(lambda_all)

#mean(lambda_all)
#mean(rowSums(ac.lambda))/2

##但是其实还有一种可能就是Resampling做太多了，导致了最后一直有degeneracy
save.image("/public1/home/scf0347/ResampFreq/GaussianMixture/Base/Base_sys.RData")