#-*- coding:utf-8 -*-
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

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

## doSNOW progress bar
m = 48
pb = txtProgressBar(max = m,style=3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress=progress)

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
x.estimate.ori = rep(0,t)

doit <- function(l){
  set.seed(l)
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
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = array(rep(0,t*m*length(threshold)),c(m,t,length(threshold)))
  rejuvs = rep(0,length(threshold))
  
  for (l in 1:m) {
    for (j in 1:length(threshold)){
      X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
      W = dnorm(y[1],0,beta*exp(X[,1]/2))
      w = W/sum(W)
      x.estimate.ss[l,1,j] = sum(w*X[,1]^3)
      if (1/sum(w^2) < threshold[j]*n){
        idx = systematic(w)
        X = X[idx,]
        W = rep(1/n,n)
        rejuvs[j] = rejuvs[j]+1
      }
      for (i in 2:t) {
        X[,i]=rnorm(n,alpha*X[,i-1],sigma)
        W = W * dnorm(y[i],0,beta*exp(X[,i]/2))
        w = W/sum(W)
        x.estimate.ss[l,i,j] = sum(w*X[,i]^3)
        if (1/sum(w^2)<threshold[j]*n){
          idx = systematic(w)
          X = X[idx,]
          W = rep(1/n,n)
          rejuvs[j] = rejuvs[j] + 1
        }
      }
    }
  }
  mse.ss <- matrix(rep(0,length(threshold)*t),nrow=length(threshold))
  for (k in 1:length(threshold)){
    for (i in 1:m)
    {
      mse.ss[k,] = mse.ss[k,]+(x^3-x.estimate.ss[i,,k])^2
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

load("SVSO3.RData")
mse = matrix(rep(0,p*length(threshold)),nrow=p)
rejuvs = matrix(rep(0,p*length(threshold)),nrow=p)
for (i in 1:p) {
  mse[i,] = res[i,]$mse
  rejuvs[i,] = res[i,]$rejuvs
}
mse.mean = colMeans(mse)
rejuvs.mean = colMeans(rejuvs)
ess = seq(0.1,1,by=.1)
plot(ess, mse.mean,ylab="MSE",xlab="ESS Threshold")

save.image("/public1/home/scf0347/ResampFreq/SV/SVSO3.RData")