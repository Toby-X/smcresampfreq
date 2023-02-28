#-*- coding:utf-8 -*-
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

##INITIAL VALUES
p=50##dimension
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
    x.estimate[1,j] = sum(w*X[,1,j])
    if (1/sum(w^2) < threshold[j]*n){
      idx = sample(1:n,n,replace = T,prob = w)
      X = X[idx,,]
      rejuvs[j] = rejuvs[j] + 1
      W = rep(1,n)
    }
    for (t in 2:p) {
      X[,t,j] = rnorm(n)
      u = dens.fun(cbind(X[,1:t,j]))/dens.fun(cbind(X[,1:(t-1),j]))/dnorm(X[,t,j])
      W = W*u
      w = W/sum(W)
      x.estimate[t,j] = sum(w*X[,t,j])
      if (1/sum(w^2) < threshold[j]*n){
        idx = sample(1:n,n,replace = T,prob = w)
        X = X[idx,,]
        rejuvs[j] = rejuvs[j] + 1
        W = rep(1,n)
      }
    }
  }
  return(list(estimate=x.estimate,rejuvs=rejuvs))
}

res = foreach(l=1:m,.combine = rbind,
              .options.snow=opts) %dopar% {
  return(doit(l))
}

# load("HDM1.RData")
# X = array(rep(0,m*p*length(threshold)),c(m,p,length(threshold)))
# rejuvs = matrix(rep(0,m*length(threshold)),nrow = m)
# for (i in 1:m) {
#   X[i,,] = res[i,]$estimate
#   rejuvs[i,] = res[i,]$rejuvs
# }
# colMeans(rejuvs)
# mse = rep(0,length(threshold))
# for (j in 1:length(threshold)) {
#   mse[j] = sum(X[,,j]^2)
# }
# mse = mse/m
# plot(ess,mse,xlab="ESS Threshold",ylab="MSE")
# boxplot(X[,,5])

save.image("/public1/home/scf0347/ResampFreq/HD/HDM.RData")