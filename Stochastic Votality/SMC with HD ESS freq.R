############################################################################
### EXAMPLE 1 SMC with rejuvenate for high-dimensional model
############################################################################
##INITIAL VALUES
p=20
n=1000
m=2000
X=matrix(NA,nrow=n,ncol=p)
W=rep(1,n)/n
N.eff=rep(NA,p)
N.eff.mean = NULL
ft1xt1=W
rejuvs=0
Rej=matrix(rep(0,p*m),nrow=m)
alpha=seq(from=0,to=1,by=0.1)
var.x=matrix(rep(0,length(alpha)*p),nrow=length(alpha))
x.estimate = array(rep(0,p*m*length(alpha)),c(m,p,length(alpha)))
var.tmp = 0
var.sum = rep(0,length(alpha))

##FUNCTIONS
dens.fun <- function(x){
  temp = function(x) {   exp(-(abs(sqrt(sum(x^2)))^3)/3) }
  c(apply(x,1,temp))
}

#w is the normalized weight
resampling <- function(w){
  u=rep(0,n)
  u[1] = runif(1,0,1/n)
  for (i in 2:n) {
    u[i] = u[1] + (i-1)/n
  }
  s = w[1]
  m = 1
  a = rep(0,n)
  for (i in 1:n) {
    while(s<u[i]){
      m = m+1
      s = s + w[m]
    }
    a[i] = m
  }
  return(a)
}

for (l in 1:m){
  for (k in 1:length(alpha)){
    ##MAIN LOOP
    X[,1] = rnorm(n)
    W = dens.fun(cbind(X[,1]))/dnorm(X[,1])
    w = W/sum(W)
    x.estimate[l,1,k] = sum(w*X[,1])
    N.eff[1] = 1/sum(w^2)
    if (N.eff[1] < alpha[k]*n){
      idx = sample(1:n,n,replace=T,prob=w)
      X = X[idx,]
      rejuvs = rejuvs + 1
      W = rep(1/n,n)
    }
    for (t in 2:p) {
      X[,t] = rnorm(n)
      u = dens.fun(cbind(X[,1:t]))/dens.fun(cbind(X[,1:(t-1)]))/dnorm(X[,t])
      W = W*u
      w = W/sum(W)
      x.estimate[l,t,k] = sum(w*X[,t])
      N.eff[t] = 1/sum(w^2)
      if (N.eff[t] < alpha[k]*n){
        idx = sample(1:n,n,replace=T,prob=w)
        X = X[idx,]
        rejuvs = rejuvs + 1
        W = rep(1/n,n)
      }
    }
    
    ##Result
    w = W/sum(W)
    Rej[l,k] = rejuvs
    rejuvs=0
  }
}

mse <- matrix(rep(0,length(alpha)*p),nrow=length(alpha))
for (k in 1:length(alpha)){
  for (i in 1:m)
  {
    mse[k,] = mse[k,]+x.estimate[i,,k]^2
  }
}
mse = mse/m
mse.sum = rowSums(mse)

plot(mse.sum)
for(i in 1:length(alpha)){
  boxplot(x.estimate[,,i])
}
