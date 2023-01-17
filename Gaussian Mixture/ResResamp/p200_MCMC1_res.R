#-*- coding:utf-8 -*-
# change into .1 variance proposal for every random walk
path1 = "/public1/home/scf0347/ResampFreq/GaussianMixture/mixture.dat"
mixture.dat = read.table(path1,header=TRUE)
y = mixture.dat$y

## libraries
library(Boom)
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

## Given Value
n = 500#number of particles
p = 200#number of distributions to sample from using smc
m = 96#number of repetitive experiments
threshold = seq(0.1,1,length=10)

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

# resampling scheme using residual resampling
resresample <- function(w){
  n = length(w)
  a = rep(0,n)
  nw = n*w
  intpart = floor(nw)
  sip = sum(intpart)
  res = nw - intpart
  sres = n - sip
  a[1:sip] = rep(1:n,intpart)
  if (sres > 0){
    a[(sip+1):n] = sample(1:n,sres,prob=res/sres)
  }
  return(a)
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

doit <- function(l){
  set.seed(l)
  mu = array(rep(0,n*p*2*length(threshold)),c(n,p,2,length(threshold)))# the first layer is mu1, the second layer is mu2, the same is as follows
  lambda = array(rep(0,n*p*2*length(threshold)),c(n,p,2,length(threshold)))
  omega = array(rep(0,n*p*2*length(threshold)),c(n,p,2,length(threshold)))
  mse.mu = rep(0,length(threshold))
  mse.omega = rep(0,length(threshold))
  mse.omega = rep(0,length(threshold))
  mse.lambda = rep(0,length(threshold))
  for (j in 1:length(threshold)) {
    w_u = rep(1,n)
    w_u.tmp = rep(1,n)
    for (i in 1:n){
      mu[i,1,,j] = rnorm(2,kexi,K^(-1/2))
      lambda[i,1,,j] = rgamma(2,alpha,beta)
      omega[i,1,,j] = rdirichlet(1,c(delta,delta))
    }
    
    ## estimate weights for step 1
    for (i in 1:n) {
      w_u[i] = exp(log.fn(1,mu[i,1,,j],lambda[i,1,,j],omega[i,1,,j]))/prod(dnorm(mu[i,1,,j],kexi,K^(-1/2)))/prod(dgamma(lambda[i,1,,j],alpha,beta))/ddirichlet(omega[i,1,,j],c(delta,delta))
    }
    if (any(is.na(w_u))){
      idx.na = is.na(w_u)
      w_u[idx.na]=0
    }
    w_n = w_u/sum(w_u)
    if (1/sum(w_n^2) < threshold[j]*n){
      idx = resresample(w_n)
      lambda = lambda[idx,,,]
      mu = mu[idx,,,]
      omega = omega[idx,,,]
      w_n = rep(1/n,n)
    }
    
    ## MAIN LOOP
    for (i in 2:p){
      ## MCMC Move
      mu.update = array(rnorm(n*2,0,.1),c(n,2))
      mu.tmp = mu[,i-1,,j]+mu.update
      lambda.update = array(rlnorm(n*2,0,.1),c(n,2))
      lambda.tmp = lambda[,i-1,,j]*lambda.update
      u.update = rnorm(n,0,.1)
      u.tmp = log(omega[,i-1,1,j]/(1-omega[,i-1,1,j]))+u.update
      omega.tmp = omega[,i-1,,j]
      omega.tmp[,1] = exp(u.tmp)/(1+exp(u.tmp))
      omega.tmp[,2] = 1-omega.tmp[,1]
      # 加个循环把每个样本都MCMC update了
      ## 这里由于三个一起update所以效益可能会低，要不要查查它update一共多少次
      for (k in 1:n) {
        acrate.mu = acceptance.mu(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,,j],lambda[k,i-1,,j],omega[k,i-1,,j])
        acrate.lambda1 = acceptance.lambda1(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,,j],lambda[k,i-1,,j],omega[k,i-1,,j])
        acrate.lambda2 = acceptance.lambda2(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,,j],lambda[k,i-1,,j],omega[k,i-1,,j])
        acrate.omega = acceptance.omega(i,mu.tmp[k,],lambda.tmp[k,],omega.tmp[k,],mu[k,i-1,,j],lambda[k,i-1,,j],omega[k,i-1,,j])
        rand = runif(1)
        if (is.na(acrate.mu)){
          mu[k,i,,j] = mu.tmp[k,]
        }else if (rand<=acrate.mu){
          mu[k,i,,j] = mu.tmp[k,]
        }else{
          mu[k,i,,j] = mu[k,i-1,,j]
        }
        if (is.na(acrate.lambda1)){
          lambda[k,i,1,j] = lambda.tmp[k,1]
        }else if (rand<=acrate.lambda1){
          lambda[k,i,1,j] = lambda.tmp[k,1]
        }else{
          lambda[k,i,1,j] = lambda[k,i-1,1,j]
        }
        if (is.na(acrate.lambda2)){
          lambda[k,i,2,j] = lambda.tmp[k,2]
        }else if (rand<=acrate.lambda2){
          lambda[k,i,2,j] = lambda.tmp[k,2]
        }else{
          lambda[k,i,2,j] = lambda[k,i-1,2,j]
        }
        if (is.na(acrate.omega)){
          omega[k,i,,j] = omega.tmp[k,]
        }else if (rand<=acrate.omega){
          omega[k,i,,j] = omega.tmp[k,]
        }else{
          omega[k,i,,j] = omega[k,i-1,,j]
        }
        divident=0
        for (o in 1:n) {
          divident = divident + exp(log.fn(i-1,mu[o,i-1,,j],lambda[o,i-1,,j],omega[o,i-1,,j]))*Kn(mu[k,i,,j],lambda[k,i,,j],omega[k,i,,j],mu[o,i-1,,j],lambda[o,i-1,,j],omega[o,i-1,,j])
        }
        w_u.tmp[k] = exp(log.fn(i,mu[k,i,,j],lambda[k,i,,j],omega[k,i,,j]))/divident
      }
      w_u = w_n*w_u.tmp
      if (any(is.na(w_u))){
        idx.na = is.na(w_u)
        w_u[idx.na]=0
      }
      w_n = w_u/sum(w_u)
      if (1/sum(w_n^2)<threshold[j]*n){
        idx = resresample(w_n)
        lambda = lambda[idx,,,]
        mu = mu[idx,,,]
        omega = omega[idx,,,]
        w_n = rep(1/n,n)
      }
    }
    
    idx1 = mu[,50,1,j]>8
    idx2 = mu[,50,2,j]>8
    mu10 = c(mu[idx1,50,1,j],mu[idx2,50,2,j])
    mu7 = c(mu[!idx1,50,1,j],mu[!idx2,50,2,j])
    
    omega10 = c(omega[idx1,50,1,j],omega[idx2,50,2,j])
    omega7 = c(omega[!idx1,50,1,j],omega[!idx2,50,2,j])
    
    lambda_all = c(lambda[,50,1,j],lambda[,50,2,j])
    
    mse.mu[j] = (mean(mu7)-7)^2+var(mu7)+(mean(mu10)-10)^2+var(mu10)
    mse.omega[j] = (mean(omega10)-0.3)^2+var(omega10)
    mse.lambda[j] = (mean(lambda_all)-4)^2+var(lambda_all)
  }
  return(list(mse.mu=mse.mu,mse.lambda=mse.lambda,mse.omega=mse.omega,mean.mu7=mean(mu7),mean.omega10=mean(omega10),mean.lambda=mean(lambda_all)))
}

res = foreach (l = 1:m,.combine = rbind,.packages = "Boom",
               .options.snow=opts) %dopar% {
                 return(doit(l))
               }
close(pb)

# load("p50MCMC1par.RData")
mse.mu = rep(0,length(threshold))
mse.lambda = rep(0,length(threshold))
mse.omega = rep(0,length(threshold))
for (i in 1:m){
  mse.mu = mse.mu + res[i,]$mse.mu
  mse.lambda = mse.lambda + res[i,]$mse.lambda
  mse.omega = mse.omega + res[i,]$mse.omega
}
mse.mu = mse.mu/m
mse.omega = mse.omega/m
mse.lambda = mse.lambda/m
# plot(mse.mu)
# plot(mse.lambda)
# plot(mse.omega)
mse.tot = mse.mu+mse.lambda+mse.omega
# plot(mse.tot)

save.image("/public1/home/scf0347/ResampFreq/GaussianMixture/SysResample/p2000MCMC1res.RData")
