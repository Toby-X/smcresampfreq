t = 50##dimension
n = 500##number of particles
N = 1e6
m = 200##first parallel to yield mse
p = 192## parallel experiments
alpha=0.91
sigma=1.0
beta=0.5
W=rep(1,n)/n
threshold = c(0)

x = rep(0,t)
y = rep(0,t)
X = matrix(rep(0,n*t),nrow=n)
X0 = matrix(rep(0,N*t),nrow=N)
x.estimate.ss1 = rep(0,t)
x.estimate.ss2 = rep(0,t)

set.seed(2)
##Generate Samples
v = rnorm(t)
u = rnorm(t)
x[1] = rnorm(1,0,sigma^2/(1-alpha^2))
y[1] = beta*exp(x[1]/2)*u[1]
for (i in 2:t){
  x[i] = alpha*x[i-1]+sigma*v[i]
  y[i] = beta*exp(x[i]/2)*u[i]
}

X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
W = dnorm(y[1],0,beta*exp(X[,1]/2))
w = W/sum(W)
x.estimate.ss1[1] = sum(w*X[,1])
idx = sample(1:n,n,replace = T,prob = w)
X = X[idx,]
W = rep(1,n)/n
for (i in 2:t) {
  X[,i]=rnorm(n,alpha*X[,i-1],sigma)
  lW.tmp = log(dnorm(y[i],0,beta*exp(X[,i]/2)))
  W = W*exp(lW.tmp)
  w = W/sum(W)
  x.estimate.ss1[i] = sum(w*X[,i])
  idx = sample(1:n,n,replace = T,prob = w)
  X = X[idx,]
  W = rep(1,n)/n
}

X[,1] = rnorm(n,0,sqrt(sigma^2/(1-alpha^2)))
W = dnorm(y[1],0,beta*exp(X[,1]/2))
w = W/sum(W)
x.estimate.ss2[1] = sum(w*X[,1])
for (i in 2:t) {
  X[,i]=rnorm(n,alpha*X[,i-1],sigma)
  lW.tmp = log(dnorm(y[i],0,beta*exp(X[,i]/2)))
  W = W*exp(lW.tmp)
  w = W/sum(W)
  x.estimate.ss2[i] = sum(w*X[,i])
}

library(tidyverse)
data = data.frame(1:length(x),x,x.estimate.ss1,x.estimate.ss2)
colnames(data) = c("idx","ori","Resamp","NoResamp")
data.n = data %>% pivot_longer("ori":"NoResamp",names_to = "cat")


ggplot(data.n)+
  geom_line(aes(idx,value,col=cat),linewidth=1)
