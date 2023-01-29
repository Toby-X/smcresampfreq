#-*- coding:utf-8 -*-
t = 20 # length of SAW
n_particles = 50 # number of particles
d = 2

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

X = array(rep(0,n_particles*t*d),c(n_particles,t,d))
W = rep(1,n_particles)
W.count = rep(1,t)
W.update = rep(1,t)
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
  ind = which(W.update==0)
  W.count[i] = mean(W.update)
  if (1/sum(w^2)<n_particles){
    idx = systematic(w)
    X = X[idx,,]
    W = rep(1,n_particles)
    # rejuvs[j] = rejuvs[j]+1
  }
}
# W.count
# prod(W.count)
W.count
cumprod(W.count)

load("SAW.RData")
cn0 = prod(W.count)

save.image("/public1/home/scf0347/ResampFreq/SAW/SAW.RData")