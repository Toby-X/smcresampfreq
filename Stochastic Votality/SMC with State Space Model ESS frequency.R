############################################################################
### EXAMPLE 2 SMC with rejuvenate for state space model
############################################################################
##INITIAL VALUES
t = 50##dimension
n = 500##number of particles
N = 1e6
m = 1e3##parallel experiments
x = rep(0,t)
y = rep(0,t)
alpha=0.91
sigma=1.0
beta=0.5
W=rep(1,n)/n
threshold = seq(.3,1,by=0.1)
X = matrix(rep(0,n*t),nrow=n)
X0 = matrix(rep(0,N*t),nrow=N)
rejuvs.ss = 0
x.estimate.ss = array(rep(0,t*m*length(threshold)),c(m,t,length(threshold)))
x.estimate.ori = rep(0,t)
Rej.ss=matrix(rep(0,length(threshold)*m),nrow=m)


set.seed(4)
##Generate Samples
v = rnorm(t)
u = rnorm(t)
x[1] = rnorm(1,0,sigma^2/(1-alpha^2))
y[1] = beta*exp(x[1]/2)*u[1]
for (i in 2:t){
  x[i] = alpha*x[i-1]+sigma*v[i]
  y[i] = beta*exp(x[i]/2)*u[i]
}

## create estimate mean
X0[,1] = rnorm(N)
W0 = dnorm(X0[,1],0,sigma^2/(1-alpha^2))*dnorm(y[i],0,beta*exp(X0[,i]/2))/dnorm(X0[,1])
w0 = W0/sum(W0)
x.estimate.ori[1] = sum(w0*X0[,1])
if (1/sum(w0^2) < N){
  idx = sample(1:N,N,prob = w0, replace = T)
  X0 = X0[idx,]
  rejuvs.ss = rejuvs.ss + 1
  W0 = rep(1/N,N)
}
for (i in 2:t) {
  X0[,i]=rnorm(N,X0[,i-1]+y[i])
  W0 = W0*dnorm(X0[,i],alpha*X0[,i-1],sigma^2)*dnorm(y[i],0,beta*exp(X0[,i]/2))/dnorm(X0[,i],X0[,i-1]+y[i])
  w0 = W0/sum(W0)
  x.estimate.ori[i] = sum(w0*X0[,i])
  if (1/sum(w0^2)<N){
    idx = sample(1:N,N,prob = w0, replace = T)
    X0 = X0[idx,]
    W0 = rep(1/n,n)
    rejuvs.ss = rejuvs.ss+1
  }
}

set.seed(NULL)
##using proposals N(xn-1+yn,1)
for (l in 1:m) {
  for (j in 1:length(threshold)){#change of series
    X[,1] = rnorm(n)
    W = dnorm(X[,1],0,sigma^2/(1-alpha^2))*dnorm(y[i],0,beta*exp(X[,i]/2))/dnorm(X[,1])
    w = W/sum(W)
    x.estimate.ss[l,1,j] = sum(w*X[,1])
    if (1/sum(w^2) < threshold[j]*n){
      idx = sample(1:n,n,prob = w, replace = T)
      X = X[idx,]
      rejuvs.ss = rejuvs.ss + 1
      W = rep(1/n,n)
    }
    for (i in 2:t) {
      X[,i]=rnorm(n,X[,i-1]+y[i])
      W = W*dnorm(X[,i],alpha*X[,i-1],sigma^2)*dnorm(y[i],0,beta*exp(X[,i]/2))/dnorm(X[,i],X[,i-1]+y[i])
      w = W/sum(W)
      x.estimate.ss[l,i,j] = sum(w*X[,i])
      if (1/sum(w^2)<threshold[j]*n){
        idx = sample(1:n,n,prob = w, replace = T)
        X = X[idx,]
        W = rep(1/n,n)
        rejuvs.ss = rejuvs.ss+1
      }
    }
    w = W/sum(W)
    Rej.ss[l,j] = rejuvs.ss
    rejuvs.ss = 0
  }
}

##Estimates using posterior mean
mse.ss <- matrix(rep(0,length(threshold)*t),nrow=length(threshold))
for (k in 1:length(threshold)){
  for (i in 1:m)
  {
    mse.ss[k,] = mse.ss[k,]+(x.estimate.ori-x.estimate.ss[i,,k])^2
  }
}
mse.ss = mse.ss/m
mse.sum.ss <- rowSums(mse.ss)

##Estimates
mse.ss <- matrix(rep(0,length(threshold)*t),nrow=length(threshold))
for (k in 1:length(threshold)){
  for (i in 1:m)
  {
    mse.ss[k,] = mse.ss[k,]+(x-x.estimate.ss[i,,k])^2
  }
}
mse.ss = mse.ss/m
mse.sum.ss <- rowSums(mse.ss)

mse <- data.frame(mse,mse.sum.ss)

plot(mse.sum.ss[-1])
for(i in 1:length(threshold)){
  boxplot(x.estimate.ss[,,i])
}


var.x <- c(var.x,var(x))
var.y <- c(var.y,var(y))

plot(x,ylim=c(-10,10),type="l")
par(new=T)
plot(x.estimate.ss[4,,11],col="red",ylim=c(-10,10),type="l")

## Plot
colnames(mse) <- c("a1","a2","b1","b2","c1","c2")
idx = 1:length(mse[,1])
mse <- data.frame(idx,mse)
mse <- mse[-1,]
library(tidyverse)
mse.new <- mse %>% pivot_longer(cols = "a1":"c2",names_to = "type", values_to = "value")
ggplot(mse.new)+
  geom_point(aes(x=idx,y=value,col=type))
ggplot(mse)+
  geom_point(aes(x=idx,y=c2))
var.y
