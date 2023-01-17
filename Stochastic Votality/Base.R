t = 50##dimension
n = 500##number of particles
N = 1e6
m = 200##first parallel to yield mse
p = 192## parallel experiments
alpha=0.91
sigma=1.0
beta=0.5
W=rep(1,n)/n
threshold = c(1)

x = rep(0,t)
y = rep(0,t)
X = matrix(rep(0,n*t),nrow=n)
X0 = matrix(rep(0,N*t),nrow=N)
x.estimate.ss = array(rep(0,t*m*length(threshold)),c(m,t,length(threshold)))
x.estimate.ori = rep(0,t)

l = 14
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
      idx = sample(1:n,n,prob = w, replace = T)
      X = X[idx,]
      W = rep(1/n,n)
    }
    for (i in 2:t) {
      X[,i]=rnorm(n,alpha*X[,i-1],sigma)
      lW.tmp = log(dnorm(y[i],0,beta*exp(X[,i]/2)))
      # if (all(lW.tmp<log(1e-200))){
      #   lW.tmp = lW.tmp + log(1e300)
      # }
      # if (all(exp(lW.tmp)==0)){
      #   cat("i=",i)
      # }
      W = W*exp(lW.tmp)
      # if (all(W<1e-50)){
      #   W = W*1e100
      # }
      # if (any(is.na(W))){
      #   idx = is.na(W)
      #   W[idx] = 0
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

save.image("Base.RData")

load("Base.RData")
mse.ss
mse.sum.ss
boxplot(mse.ss)