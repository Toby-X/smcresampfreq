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
N = 2e5#control group
m = 200##first parallel to yield mse
p = 240## parallel experiments
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

ess_smc_m <- function(n,t,alpha,beta,sigma,y){
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = rep(0,t)
  rejuvs = 0
  
  X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
  W = dnorm(y[1],0,beta*exp(X[,1]/2))
  w = W/sum(W)
  x.estimate.ss[1] = sum(w*X[,1])
  # Resampling
  idx = sample(1:n,n,prob = w, replace = T)
  X = X[idx,]
  W = rep(1/n,n)
  rejuvs = rejuvs+1

  for (i in 2:t) {
    X[,i]=rnorm(n,alpha*X[,i-1],sigma)
    W = W * dnorm(y[i],0,beta*exp(X[,i]/2))
    w = W/sum(W)
    x.estimate.ss[i] = sum(w*X[,i])
    idx = sample(1:n,n,prob = w, replace = T)
    X = X[idx,]
    W = rep(1/n,n)
    rejuvs = rejuvs + 1
  }
  return(list(estimate=x.estimate.ss,rejuvs=rejuvs))
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
  
  for (i in 2:t) {
    X[,i]=rnorm(n,alpha*X[,i-1],sigma)
    W = W * dnorm(y[i],0,beta*exp(X[,i]/2))
    w = W/sum(W)
    x.estimate.ss[i] = sum(w*X[,i])
    idx = systematic(w)
    X = X[idx,]
    W = rep(1/n,n)
    rejuvs = rejuvs + 1
  }
  return(list(estimate=x.estimate.ss,rejuvs=rejuvs))
}

ess_smc_r <- function(n,t,alpha,beta,sigma,y){
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = rep(0,t)
  rejuvs = 0
  
  X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
  W = dnorm(y[1],0,beta*exp(X[,1]/2))
  w = W/sum(W)
  x.estimate.ss[1] = sum(w*X[,1])
  # Resampling
  idx = resresample(w)
  X = X[idx,]
  W = rep(1/n,n)
  rejuvs = rejuvs+1
  
  for (i in 2:t) {
    X[,i]=rnorm(n,alpha*X[,i-1],sigma)
    W = W * dnorm(y[i],0,beta*exp(X[,i]/2))
    w = W/sum(W)
    x.estimate.ss[i] = sum(w*X[,i])
    idx = resresample(w)
    X = X[idx,]
    W = rep(1/n,n)
    rejuvs = rejuvs + 1
  }
  return(list(estimate=x.estimate.ss,rejuvs=rejuvs))
}

# seq_dens <- function(X,alpha,beta,sigma,current,y){
#   ## t is the current step of the sampling procedure
#   ldens = rep(0,nrow(X))
#   if (current==1){
#       ldens = dnorm(X[,current+1],alpha*X[,current],sigma,log = T)+dnorm(y[current+1],0,beta*exp(X[,current+1]/2),log = T)
#   } else if(current == 2){
#     ldens = dnorm(X[,current-1],0,sigma/sqrt(1-alpha^2),log = T)+dnorm(X[,current+1],alpha*X[,current],sigma,log = T)+
#       dnorm(y[current-1],0,beta*exp(X[,current-1]/2),log = T)+dnorm(y[current+1],0,beta*exp(X[,current+1]/2),log = T)
#   } else{
#     ldens = dnorm(X[,current-1],alpha*X[,current-2],sigma,log = T)+dnorm(X[,current+1],alpha*X[,current],sigma,log = T)+
#       dnorm(y[current-1],0,beta*exp(X[,current-1]/2),log = T)+dnorm(y[current+1],0,beta*exp(X[,current+1]/2),log = T)
#   }
#   
#   return(exp(ldens))
# }

var_smc_m <- function(n,t,alpha,beta,sigma,y){
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
    idx = sample(1:n,n,prob = w, replace = T)
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
      idx = sample(1:n,n,prob = w, replace = T)
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

var_smc_s <- function(n,t,alpha,beta,sigma,y){
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

var_smc_r <- function(n,t,alpha,beta,sigma,y){
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
    idx = resresample(w)
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
    # 这里是要加上weight的
    cri[i] = sum(w*W.tmp*(W/mean(W)*W.tmp/mean(W.tmp)-W.tmp/mean(W.tmp)))
    if (cri[i] > 0){
      idx = resresample(w)
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

doit <- function(l){
  if (l == 14|| l == 63 || l == 80 || l==69|| l == 131|| l==224){
    set.seed(l^3)
  }else{
    set.seed(l)
  }
  data = sv_gen(t,alpha,beta,sigma)
  x = data$x
  y = data$y
  mse.m = 0
  mse.r = 0
  mse.s = 0
  mse.new.m = 0
  mse.new.r = 0
  mse.new.s = 0
  rejuvs = rep(0,6)
  
  con.group = ess_smc_s(N,t,alpha,beta,sigma,y)
  for (i in 1:m) {
    ess.m = ess_smc_m(n,t,alpha,beta,sigma,y)
    ess.r = ess_smc_r(n,t,alpha,beta,sigma,y)
    ess.s = ess_smc_s(n,t,alpha,beta,sigma,y)
    new.m = var_smc_m(n,t,alpha,beta,sigma,y)
    new.r = var_smc_r(n,t,alpha,beta,sigma,y)
    new.s = var_smc_s(n,t,alpha,beta,sigma,y)

    mse.m = mse.m + (ess.m$estimate-con.group$estimate)^2
    mse.r = mse.r + (ess.r$estimate-con.group$estimate)^2
    mse.s = mse.s + (ess.s$estimate-con.group$estimate)^2
    mse.new.m = mse.new.m + (new.m$estimate-con.group$estimate)^2
    mse.new.r = mse.new.r + (new.r$estimate-con.group$estimate)^2
    mse.new.s = mse.new.s + (new.s$estimate-con.group$estimate)^2

    rejuvs[1] = ess.m$rejuvs + rejuvs[1]
    rejuvs[2] = ess.r$rejuvs + rejuvs[2]
    rejuvs[3] = ess.s$rejuvs + rejuvs[3]
    rejuvs[4] = new.m$rejuvs + rejuvs[4]
    rejuvs[5] = new.r$rejuvs + rejuvs[5]
    rejuvs[6] = new.s$rejuvs + rejuvs[6]
  }
  mse.m = sum(mse.m/m)
  mse.r = sum(mse.r/m)
  mse.s = sum(mse.s/m)
  mse.new.m = sum(mse.new.m/m)
  mse.new.r = sum(mse.new.r/m)
  mse.new.s = sum(mse.new.s/m)
  rejuvs = rejuvs/m
  
  return(list(mse.m = mse.m,mse.r = mse.r, mse.s = mse.s, new.m = mse.new.m, new.r = mse.new.r, new.s = mse.new.s, rejuvs = rejuvs))
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

save.image("/public1/home/scf0347/ResampFreq/SV/con_group.RData")
