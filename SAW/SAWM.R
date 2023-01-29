#-*- coding: utf-8 -*-
# simulate SAW using growth method for a 2D template
t = 10 # length of SAW
n_particles = 200 # number of particles
N = 1e6
d = 2
m = 192 # parallel experiments
threshold = seq(0,1,by=.1)

# libraries
library(parallel)
library(foreach)
library(doSNOW)
numCores = detectCores()
cl = makeCluster(numCores)
registerDoSNOW(cl)

pb = txtProgressBar(max = m,style=3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress=progress)



find_next <- function(p, x){
  # x is a sequence of 2D, t is the current time
  nxt = rbind(c(1,0),c(-1,0),c(0,1),c(0,-1))
  nxt.r = matrix(c(0,0),nrow=1)
  for (i in 1:nrow(nxt)) {
    x.tmp = x[p-1,] + nxt[i,]
    if (!any(apply(x, 1, setequal, y = x.tmp))){
      nxt.r = rbind(nxt.r,nxt[i,])
    }
  }
  nxt.r = nxt.r[-1,]
  return(nxt.r)
}

doit <- function(l){
  set.seed(l)
  x.estimate = array(rep(0,t*length(threshold)*d),c(t,d,length(threshold)))
  rejuvs = rep(0,length(threshold))
  cn = rep(1,length(threshold))
  for (j in 1:length(threshold)) {
    X = array(rep(0,n_particles*t*d),c(n_particles,t,d))
    W = rep(1/n_particles,n_particles)
    W.count = rep(1,t)
    W.update = rep(1,n_particles)
    for (i in 2:t) {
      for (k in 1:n_particles) {
        avail = find_next(i,as.matrix(X[k,,]))
        nt = nrow(avail)
        if (is.null(nt)){
          W[k] = W[k]
          W.update[k] = 1
          X[k,i,] = X[k,i-1,]+avail
        }else if (nt != 0){
          W[k] = W[k] * nt
          W.update[k] = nt
          nid = sample(1:nt,1)
          X[k,i,] = X[k,i-1,]+avail[nid,]
        }else{
          W[k] = 0
          W.update[k] = 0
          X[k,i,] = X[k,i-1,]
        }
      }
      w = W/sum(W)
      W.count[i] = mean(W.update)
      if (1/sum(w^2)<threshold[j]*n_particles){
        idx = sample(1:n_particles,n_particles,replace = T,prob = w)
        X = X[idx,,]
        W = rep(1/n_particles,n_particles)
        rejuvs[j] = rejuvs[j]+1
      }
    }
    w = W/sum(W)
    x.estimate[,1,j] = w%*%as.matrix(X[,,1])
    x.estimate[,2,j] = w%*%as.matrix(X[,,2])
    cn[j] = prod(W.count)
  }
  return(list(estimate = x.estimate, rejuvs = rejuvs,cn=cn))
}

res = foreach(l=1:m,.combine = rbind,
              .options.snow=opts) %dopar% {
                return(doit(l))
              }

load("SAWM.RData")#cn0=11786.05
# cn0 = 1.17704242*2.638^t*t^(43/32-1)
cn0 = 44100
cn = matrix(rep(0,m*length(threshold)),nrow=m)
for (i in 1:m) {
  cn[i,] = res[i,]$cn
}
mse = rep(0,length(threshold))
for (j in 1:length(threshold)) {
  mse[j] = sum((cn[,j]-cn0)^2)
}
plot(mse)

X = array(rep(0,m*t*d*length(threshold)),c(m,t,d,length(threshold)))
rejuvs = matrix(rep(0,m*length(threshold)),nrow = m)
# # head(res)
# # str(res[1,]$estimate)
for (i in 1:m) {
  X[i,,,] = res[i,]$estimate
  rejuvs[i,] = res[i,]$rejuvs
}
colMeans(rejuvs)
# mse = rep(0,length(threshold))
# for (j in 1:length(threshold)) {
#   mse[j] = sum(X[,,,j]^2)/m
# }
# ess = seq(0,1,by=.1)
# plot(ess,mse,xlab = "ESS Threshold", ylab = "MSE")

save.image("/public1/home/scf0347/ResampFreq/SAW/SAWM.RData")