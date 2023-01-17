#-*- coding:utf-8 -*-
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

m = 192
## doSNOW progress bar
pb = txtProgressBar(max = m,style=3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress=progress)

##INITIAL VALUES
t = 50##dimension
n = 500##number of particles
m = 200##first parallel to yield mse
p = 48## parallel experiments
alpha=0.91
sigma=1.0
beta=0.5
W=rep(1,n)/n
threshold = seq(.1,1,by=0.1)

doit <- function(l,x,y){
  X = matrix(rep(0,n*t),nrow=n)
  x.estimate.ss = array(rep(0,t*m*length(threshold)),c(m,t,length(threshold)))
  x.estimate.ori = rep(0,t)
  
  for (l in 1:m) {
    for (j in 1:length(threshold)){#change of series
      X[,1] = rnorm(n)
      W = dnorm(X[,1],0,sigma^2/(1-alpha^2))*dnorm(y[i],0,beta*exp(X[,i]/2))/dnorm(X[,1])
      if (any(is.na(W))){
        idx = is.na(W)
        W[idx] = 0
      }
      w = W/sum(W)
      x.estimate.ss[l,1,j] = sum(w*X[,1])
      if (1/sum(w^2) < threshold[j]*n){
        idx = sample(1:n,n,prob = w, replace = T)
        X = X[idx,]
        W = rep(1/n,n)
      }
      for (i in 2:t) {
        X[,i]=rnorm(n,X[,i-1]+y[i])
        lW.tmp = log(dnorm(X[,i],alpha*X[,i-1],sigma^2))+log(dnorm(y[i],0,beta*exp(X[,i]/2)))-log(dnorm(X[,i],X[,i-1]+y[i]))
        if (all(lW.tmp<log(1e-100))){
          lW.tmp = lW.tmp + log(1e200)
        }
        # if (all(exp(lW.tmp)==0)){
        #   cat("i=",i)
        # }
        W = exp(log(W)+lW.tmp)
        if (any(is.infinite(lW.tmp))){
          idx.inf = is.infinite(lW.tmp)
          W[idx.inf] = 0
        }
        # if (all(W<1e-50)){
        #   W = W*1e100
        # }
        # if (any(is.na(W))){
        #   idx = is.na(W)
        #   W[idx] = 0
        # }
        # if (all(W==0)){
        #   x.estimate.ss[l,i,j] = 9999
        #   next
        # }
        w = W/sum(W)
        x.estimate.ss[l,i,j] = sum(w*X[,i])
        if (1/sum(w^2)<threshold[j]*n){
          idx = sample(1:n,n,prob = w, replace = T)
          X = X[idx,]
          W = rep(1/n,n)
        }
      }
    }
  }
  mse.ss <- matrix(rep(0,length(threshold)*t),nrow=length(threshold))
  for (k in 1:length(threshold)){
    for (i in 1:m)
    {
      mse.ss[k,] = mse.ss[k,]+(x-x.estimate.ss[i,,k])^2
    }
  }
  mse.ss = mse.ss/m
  mse.sum.ss <- rowSums(mse.ss)
  
  return(list(mse = mse.sum.ss, x.estimate=x.estimate.ss, x=x))
}

doit1 <- function(x){
  
}

res = foreach (l = 1:p,.combine = rbind,.packages = "Boom",
               .options.snow=opts) %dopar% {
                 x = rep(0,t)
                 y = rep(0,t)
                 set.seed(l)
                 ##Generate Samples
                 v = rnorm(t)
                 u = rnorm(t)
                 x[1] = rnorm(1,0,sigma^2/(1-alpha^2))
                 y[1] = beta*exp(x[1]/2)*u[1]
                 for (i in 2:t){
                   x[i] = alpha*x[i-1]+sigma*v[i]
                   y[i] = beta*exp(x[i]/2)*u[i]
                 }
                 return(doit(l,x,y))
               }
close(pb)

save.image("/public1/home/scf0347/ResampFreq/SV/SVMO.RData")