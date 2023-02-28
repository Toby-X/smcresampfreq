#-*- coding:utf-8 -*-
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

## doSNOW progress bar
m = 48
pb = txtProgressBar(max = m,style=3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress=progress)

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

##INITIAL VALUES
t = 50##dimension
n = 500##number of particles
N = 1e6
m = 200##first parallel to yield mse
p = 48## parallel experiments
x = rep(0,t)
y = rep(0,t)
alpha=0.91
sigma=1.0
beta=0.5
W=rep(1,n)/n
threshold = seq(.1,1,by=0.1)
X = matrix(rep(0,n*t),nrow=n)
X0 = matrix(rep(0,N*t),nrow=N)
x.estimate.ss = array(rep(0,t*m*length(threshold)),c(m,t,length(threshold)))
x.estimate.ori = rep(0,t)

doit <- function(l){
  set.seed(l)
  ##Generate Samples
  x = rep(0,t)
  y = rep(0,t)
  X = matrix(rep(0,n*t),nrow=n)
  X0 = matrix(rep(0,N*t),nrow=N)
  x.estimate.ss = array(rep(0,t*m*length(threshold)),c(m,t,length(threshold)))
  x.estimate.ori = rep(0,t)
  rejuvs = rep(0,length(threshold))
  
  v = rnorm(t)
  u = rnorm(t)
  x[1] = rnorm(1,0,sigma^2/(1-alpha^2))
  y[1] = beta*exp(x[1]/2)*u[1]
  for (i in 2:t){
    x[i] = alpha*x[i-1]+sigma*v[i]
    y[i] = beta*exp(x[i]/2)*u[i]
  }
  
  ## create estimate mean
  X0[,1] = rnorm(N,0,sqrt(sigma^2/(1-alpha^2)))
  W0 = dnorm(y[1],0,beta*exp(X0[,1]/2))
  w0 = W0/sum(W0)
  x.estimate.ori[1] = sum(w0*X0[,1])
  if (1/sum(w0^2) < N){
    idx = sample(1:N,N,prob = w0, replace = T)
    X0 = X0[idx,]
    W0 = rep(1/N,N)
  }
  for (i in 2:t) {
    X0[,i]=rnorm(N,alpha*X0[,i-1],sigma)
    W0 = W0*dnorm(y[i],0,beta*exp(X0[,i]/2))
    w0 = W0/sum(W0)
    x.estimate.ori[i] = sum(w0*X0[,i])
    if (1/sum(w0^2)<N){
      idx = sample(1:N,N,prob = w0, replace = T)
      X0 = X0[idx,]
      W0 = rep(1/n,n)
    }
  }
  
  
  for (l in 1:m) {
    for (j in 1:length(threshold)){#change of series
      X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
      W = dnorm(y[1],0,beta*exp(X[,1]/2))
      if (any(is.na(W))){
        idx = is.na(W)
        W[idx] = 0
      }
      w = W/sum(W)
      x.estimate.ss[l,1,j] = sum(w*X[,1])
      if (1/sum(w^2) < threshold[j]*n){
        idx = resresample(w)
        X = X[idx,]
        W = rep(1/n,n)
        rejuvs[j] = rejuvs[j] + 1
      }
      for (i in 2:t) {
        X[,i]=rnorm(n,alpha*X[,i-1],sigma)
        W = W*dnorm(y[i],0,beta*exp(X[,i]/2))
        w = W/sum(W)
        x.estimate.ss[l,i,j] = sum(w*X[,i])
        if (1/sum(w^2)<threshold[j]*n){
          idx = resresample(w)
          X = X[idx,]
          W = rep(1/n,n)
          rejuvs[j] = rejuvs+1
        }
      }
    }
  }
  mse.ss <- matrix(rep(0,length(threshold)*t),nrow=length(threshold))
  for (k in 1:length(threshold)){
    for (i in 1:m)
    {
      mse.ss[k,] = mse.ss[k,]+(x.estimate.ori-x.estimate.ss[i,,k])^2
    }
  }
  mse.ss = mse.ss/m
  mse.sum.ss <- rowSums(mse.ss)
  rejuvs = rejuvs/m
  
  return(list(mse = mse.sum.ss, rejuvs = rejuvs))
}

res = foreach (l = 1:p,.combine = rbind,.packages = "Boom",
               .options.snow=opts) %dopar% {
                 return(doit(l))
               }
close(pb)

save.image("/public1/home/scf0347/ResampFreq/SV/SVRP.RData")