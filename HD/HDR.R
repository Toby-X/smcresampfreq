#-*- coding:utf-8 -*-
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

##INITIAL VALUES
p=20##dimension
n=500##samples
m=192##parallel experiments
threshold = seq(.1,1,by=.1)

pb = txtProgressBar(max = m,style=3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress=progress)

##FUNCTIONS
dens.fun <- function(x){
  temp = function(x) {   exp(-(abs(sqrt(sum(x^2)))^3)/3) }
  c(apply(x,1,temp))
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

doit <- function(l){
  set.seed(l)
  rejuvs=rep(0,length(threshold))
  x.estimate = array(rep(0,p*length(threshold)),c(p,length(threshold)))
  X=array(rep(0,n*p*length(threshold)),c(n,p,length(threshold)))
  for (j in 1:length(threshold)){
    ##MAIN LOOP
    X[,1,j] = rnorm(n)
    W = dens.fun(cbind(X[,1,j]))/dnorm(X[,1,j])
    w = W/sum(W)
    if (1/sum(w^2) < threshold[j]*n){
      idx = resresample(w)
      X = X[idx,,]
      rejuvs[j] = rejuvs[j] + 1
      W = rep(1,n)
    }
    for (t in 2:p) {
      X[,t,j] = rnorm(n)
      u = dens.fun(cbind(X[,1:t,j]))/dens.fun(cbind(X[,1:(t-1),j]))/dnorm(X[,t,j])
      W = W*u
      w = W/sum(W)
      if (1/sum(w^2) < threshold[j]*n){
        idx = resresample(w)
        X = X[idx,,]
        rejuvs[j] = rejuvs[j] + 1
        W = rep(1,n)
      }
    }
    ##Result
    w = W/sum(W)
    x.estimate[,j] = as.vector(w%*%as.matrix(X[,,j]))
  }
  return(list(estimate=x.estimate,rejuvs=rejuvs))
}

res = foreach(l=1:m,.combine = rbind,
              .options.snow=opts) %dopar% {
                return(doit(l))
              }

load("HDR.RData")
X = array(rep(0,m*p*length(threshold)),c(m,p,length(threshold)))
rejuvs = matrix(rep(0,m*length(threshold)),nrow = m)
for (i in 1:m) {
  X[i,,] = res[i,]$estimate
  rejuvs[i,] = res[i,]$rejuvs
}
colMeans(rejuvs)
mse = rep(0,length(threshold))
for (j in 1:length(threshold)) {
  mse[j] = sum(X[,,j]^2)
}
mse = mse/m
ess = seq(.1,1,by=.1)
plot(ess,mse,ylab = "MSE",xlab = "ESS Threshold")

save.image("/public1/home/scf0347/ResampFreq/HD/HDR.RData")