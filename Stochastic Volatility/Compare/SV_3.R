#-*- coding:utf-8 -*-
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

##INITIAL VALUES
t = 50##dimension
n = 500##number of particles
m = 200##first parallel to yield mse
p = 48## parallel experiments
A = 50# number of lookaheads
alpha=0.91
sigma=1.0
beta=0.5
W=rep(1,n)/n

## doSNOW progress bar
pb = txtProgressBar(max = p,style=3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress=progress)

## FUNCTIONS
sv_gen <- function(t,alpha,beta,sigma){
  x = rep(0,t)
  y = rep(0,t)
  ##Generate Samples
  v = rnorm(t)
  u = rnorm(t)
  x[1] = rnorm(1,0,sigma^2/(1-alpha^2))
  y[1] = beta*exp(x[1]/2)*u[1]
  for (i in 2:t){
    x[i] = alpha*x[i-1]+sigma*v[i]
    y[i] = beta*exp(x[i]/2)*u[i]
  }
  dat = data.frame(x,y)
  return(dat)
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

systematic = function(w){# actually stratified
  m = length(w)
  su = (runif(m)+0:(m-1))/m
  return(inv_cdf(su,w))
}

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

ess_smc_s <- function(n,t,alpha,beta,sigma,y){
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = rep(0,t)
  rejuvs = 0
  
  X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
  W = dnorm(y[1],0,beta*exp(X[,1]/2))
  w = W/sum(W)
  x.estimate.ss[1] = sum(w*X[,1])
  # Resampling
  idx = systematic(w)
  X = X[idx,]
  W = rep(1/n,n)
  rejuvs = rejuvs+1
  
  for (i in 2:(t-1)) {
    X[,i]=rnorm(n,alpha*X[,i-1],sigma)
    W = W * dnorm(y[i],0,beta*exp(X[,i]/2))
    w = W/sum(W)
    x.estimate.ss[i] = sum(w*X[,i])
    idx = systematic(w)
    X = X[idx,]
    W = rep(1/n,n)
    rejuvs = rejuvs + 1
  }
  i=t
  X[,i]=rnorm(n,alpha*X[,i-1],sigma)
  W = W * dnorm(y[i],0,beta*exp(X[,i]/2))
  w = W/sum(W)
  x.estimate.ss[i] = sum(w*X[,i])
  return(list(estimate=x.estimate.ss,rejuvs=rejuvs))
}

var_smc_s <- function(n,t,alpha,beta,sigma,y,A){
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = rep(0,t)
  rejuvs = 0
  cri = rep(0,t)
  
  X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
  W = dnorm(y[1],0,beta*exp(X[,1]/2)) ## 这个是G而不是Gbar
  W.u = W
  w = W.u/sum(W.u)
  x.estimate.ss[1] = sum(w*X[,1])
  # Resampling
  x.tmp = rnorm(n,alpha*X[,1],sigma)
  W.tmp = dnorm(y[2],0,beta*exp(x.tmp/2))
  cri[1] = sum(w*W.tmp*(W/mean(W)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)))
  if (cri[1] > 0){
    idx = systematic(w)
    X = X[idx,]
    W.u = rep(1/n,n)
    rejuvs = rejuvs+1
  }
  
  for (i in 2:(t-1)) {
    X[,i] = rnorm(n,alpha*X[,i-1],sigma)
    W = dnorm(y[i],0,beta*exp(X[,i]/2))
    W.u = W.u*W
    w = W.u/sum(W.u)
    x.estimate.ss[i] = sum(w*X[,i])
    # Resampling
    x.tmp = rnorm(n,alpha*X[,i],sigma)
    W.tmp = dnorm(y[i+1],0,beta*exp(x.tmp/2))
    cri[i] = sum(w*W.tmp*(W/mean(W)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)))
    if (cri[i] > 0){
      idx = systematic(w)
      X = X[idx,]
      W.u = rep(1/n,n)
      rejuvs = rejuvs+1
    }
  }
  
  X[,t] = rnorm(n,alpha*X[,t-1],sigma)
  W = dnorm(y[t],0,beta*exp(X[,t]/2))
  W.u = W.u*W
  w = W.u/sum(W.u)
  x.estimate.ss[t] = sum(w*X[,t])
  
  return(list(estimate=x.estimate.ss,rejuvs=rejuvs))
}

multi_s <- function(n,t,alpha,beta,sigma,y,A){
  
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = rep(0,t)
  rejuvs = 0
  cri = rep(0,t)
  
  X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
  W = dnorm(y[1],0,beta*exp(X[,1]/2)) ## 这个是G而不是Gbar
  W.u = W
  w = W.u/sum(W.u)
  x.estimate.ss[1] = sum(w*X[,1])
  # Resampling
  W.tmp = rep(0, n)
  for (k in 1:n) {
    x.tmp = rnorm(A,alpha*X[k,1],sigma)
    W.tmp[k] = mean(dnorm(y[2],0,beta*exp(x.tmp/2)))
  }
  cri[1] = sum(w*W.tmp*(W/mean(W)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)))
  if (cri[1] > 0){
    idx = systematic(w)
    X = X[idx,]
    W.u = rep(1/n,n)
    rejuvs = rejuvs+1
  }
  
  for (i in 2:(t-1)) {
    X[,i] = rnorm(n,alpha*X[,i-1],sigma)
    W = dnorm(y[i],0,beta*exp(X[,i]/2))
    W.u = W.u*W
    w = W.u/sum(W.u)
    x.estimate.ss[i] = sum(w*X[,i])
    # Resampling
    for (k in 1:n) {
      x.tmp = rnorm(A,alpha*X[,i],sigma)
      W.tmp[k] = mean(dnorm(y[i+1],0,beta*exp(x.tmp/2)))
    }
    cri[i] = sum(w*W.tmp*(W/mean(W)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)))
    if (cri[i] > 0){
      idx = systematic(w)
      X = X[idx,]
      W.u = rep(1/n,n)
      rejuvs = rejuvs+1
    }
  }
  
  X[,t] = rnorm(n,alpha*X[,t-1],sigma)
  W = dnorm(y[t],0,beta*exp(X[,t]/2))
  W.u = W.u*W
  w = W.u/sum(W.u)
  x.estimate.ss[t] = sum(w*X[,t])
  
  return(list(estimate=x.estimate.ss,rejuvs=rejuvs))
}

pri_multi_s <- function(n,t,alpha,beta,sigma,y,A){
  
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = rep(0,t)
  rejuvs = 0
  cri = rep(0,t)
  
  X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
  W = dnorm(y[1],0,beta*exp(X[,1]/2)) ## 这个是G而不是Gbar
  W.u = W
  w = W.u/sum(W.u)
  x.estimate.ss[1] = sum(w*X[,1])
  # Resampling
  W.tmp = rep(0, n)
  eta = rep(0,n)
  for (k in 1:n) {
    x.tmp = rnorm(A,alpha*X[k,1],sigma)
    W.tmp[k] = mean(dnorm(y[2],0,beta*exp(x.tmp/2)))
    eta[k] = mean(dnorm(y[2],0,beta*exp(x.tmp/2))^2)
  }
  cri[1] = sum(w*W.tmp*(W/mean(W)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)))
  if (cri[1] > 0){
    b = W*eta^(1/2)
    b.n = b/sum(b)
    idx = systematic(b.n)
    X = X[idx,]
    W.u = W[idx]/b[idx]
    rejuvs = rejuvs+1
  }
  
  for (i in 2:(t-1)) {
    X[,i] = rnorm(n,alpha*X[,i-1],sigma)
    W = dnorm(y[i],0,beta*exp(X[,i]/2))
    W.u = W.u*W
    w = W.u/sum(W.u)
    x.estimate.ss[i] = sum(w*X[,i])
    # Resampling
    for (k in 1:n) {
      x.tmp = rnorm(A,alpha*X[,i],sigma)
      W.tmp[k] = mean(dnorm(y[i+1],0,beta*exp(x.tmp/2)))
      eta[k] = mean(dnorm(y[i+1],0,beta*exp(x.tmp/2))^2)
    }
    cri[i] = sum(w*W.tmp*(W/mean(W)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)))
    if (cri[i] > 0){
      b = W.u*eta^(1/2)
      b.n = b/sum(b)
      idx = systematic(b.n)
      X = X[idx,]
      W.u = W.u[idx]/b[idx]
      rejuvs = rejuvs+1
    }
  }
  
  X[,t] = rnorm(n,alpha*X[,t-1],sigma)
  W = dnorm(y[t],0,beta*exp(X[,t]/2))
  W.u = W.u*W
  w = W.u/sum(W.u)
  x.estimate.ss[t] = sum(w*X[,t])
  
  return(list(estimate=x.estimate.ss,rejuvs=rejuvs))
}

pri_s <- function(n,t,alpha,beta,sigma,y,A){
  
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = rep(0,t)
  rejuvs = 0
  cri = rep(0,t)
  eta = rep(0,n)
  
  X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
  W = dnorm(y[1],0,beta*exp(X[,1]/2)) ## 这个是G而不是Gbar
  W.u = W
  w = W.u/sum(W.u)
  x.estimate.ss[1] = sum(w*X[,1])
  # Resampling
  W.tmp = rep(0, n)
  for (k in 1:n) {
    x.tmp = rnorm(A,alpha*X[k,1],sigma)
    W.tmp[k] = mean(dnorm(y[2],0,beta*exp(x.tmp/2)))
    eta[k] = mean(dnorm(y[2],0,beta*exp(x.tmp/2))^2)
  }
  b = W*eta^(1/2)
  b.n = b/sum(b)
  idx = systematic(b.n)
  X = X[idx,]
  W.u = W[idx]/b[idx]
  rejuvs = rejuvs+1
  
  for (i in 2:(t-1)) {
    X[,i] = rnorm(n,alpha*X[,i-1],sigma)
    W = dnorm(y[i],0,beta*exp(X[,i]/2))
    W.u = W.u*W
    w = W.u/sum(W.u)
    x.estimate.ss[i] = sum(w*X[,i])
    # Resampling
    for (k in 1:n) {
      x.tmp = rnorm(A,alpha*X[,i],sigma)
      W.tmp[k] = mean(dnorm(y[i+1],0,beta*exp(x.tmp/2)))
      eta[k] = mean(dnorm(y[i+1],0,beta*exp(x.tmp/2))^2)
    }
    b = W.u*eta^(1/2)
    b.n = b/sum(b)
    idx = systematic(b.n)
    X = X[idx,]
    W.u = W.u[idx]/b[idx]
    rejuvs = rejuvs+1
  }
  
  X[,t] = rnorm(n,alpha*X[,t-1],sigma)
  W = dnorm(y[t],0,beta*exp(X[,t]/2))
  W.u = W.u*W
  w = W.u/sum(W.u)
  x.estimate.ss[t] = sum(w*X[,t])
  
  return(list(estimate=x.estimate.ss,rejuvs=rejuvs))
}

doit <- function(l){
  if (l == 14){
    set.seed(l^2)
  }else{
    set.seed(l)
  }
  data = sv_gen(t,alpha,beta,sigma)
  x = data$x
  y = data$y
  mse.s = 0
  mse.new.s = 0
  mse.multi.s = 0
  mse.pri.multi.s = 0
  mse.pri.s = 0
  rejuvs = rep(0,5)
  
  for (i in 1:m) {
    ess.s = ess_smc_s(n,t,alpha,beta,sigma,y)
    new.s = var_smc_s(n,t,alpha,beta,sigma,y)
    multi.s = multi_s(n,t,alpha,beta,sigma,y,A)
    pri.multi.s = pri_multi_s(n,t,alpha,beta,sigma,y,A)
    pri.s = pri_s(n,t,alpha,beta,sigma,y,A)

    mse.s = mse.s + (ess.s$estimate-x)^2
    mse.new.s = mse.new.s + (new.s$estimate-x)^2
    mse.multi.s = mse.multi.s + (multi.s$estimate-x)^2
    mse.pri.multi.s = mse.pri.multi.s + (pri.multi.s$estimate-x)^2
    mse.pri.s = mse.pri.s + (pri.s$estimate-x)^2

    rejuvs[1] = ess.s$rejuvs + rejuvs[1]
    rejuvs[2] = new.s$rejuvs + rejuvs[2]
    rejuvs[3] = multi.s$rejuvs + rejuvs[3]
    rejuvs[4] = pri.s$rejuvs + rejuvs[4]
    rejuvs[5] = pri.multi.s$rejuvs + rejuvs[5]
  }
  mse.s = sum(mse.s/m)
  mse.new.s = sum(mse.new.s/m)
  mse.multi.s = sum(mse.multi.s/m)
  mse.pri.multi.s = sum(mse.pri.multi.s/m)
  mse.pri.s = sum(mse.pri.s/m)
  rejuvs = rejuvs/m
  
  return(list(mse.s = mse.s, new.s = mse.new.s, multi.s = mse.multi.s, pri.s = mse.pri.s, multi.pri.s = mse.pri.multi.s, rejuvs = rejuvs))
}

res = foreach (l = 1:p,.combine = rbind,
               .options.snow=opts) %dopar% {
                 return(doit(l))
               }
close(pb)

# load("newmethod.RData")
# head(res)
# mse.m = rep(0,p)
# mse.r = rep(0,p)
# mse.s = rep(0,p)
# mse.new.m = rep(0,p)
# mse.new.r = rep(0,p)
# mse.new.s = rep(0,p)
# rejuvs = rep(0,6)
# for (i in 1:p) {
#   mse.m[i] = res[i,]$mse.m
#   mse.r[i] = res[i,]$mse.r
#   mse.s[i] = res[i,]$mse.s
#   mse.new.m[i] = res[i,]$new.m
#   mse.new.r[i] = res[i,]$new.r
#   mse.new.s[i] = res[i,]$new.s
#   rejuvs = rejuvs + res[i,]$rejuvs
# }
# rejuvs = rejuvs/p
# 
# sum(mse.new.m<mse.m)
# sum(mse.new.r<mse.r)
# sum(mse.new.s<mse.s)
# 
# mean(mse.m)
# mean(mse.new.m)
# mean(mse.r)
# mean(mse.new.r)
# mean(mse.s)
# mean(mse.new.s)

save.image("/public1/home/scf0347/ResampFreq/SV/Compare/newmethod.RData")
