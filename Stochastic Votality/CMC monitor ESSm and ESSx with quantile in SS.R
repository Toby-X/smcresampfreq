############################################################################
### EXAMPLE 3 CMC with rejuvenate for state space model
############################################################################
##CHOICE OF PARTITION: MAYBE QUANTILE, RANDOM, EQUAL DISTANCE(IN PAPER, RANDOM BETTER THAN EQUAL DIST)
##INITIAL VALUES
t = 50
n = 1e3
m = 1e3
partition = c(0.02,0.04,0.05,0.1,0.5,1.0)
x = rep(0,t)
y = rep(0,t)
alpha=0.91
sigma=1.0
beta=0.5
W=rep(1,n)/n
wm = W

threshold = seq(0,1,by=0.1)
X = matrix(rep(0,n*t),nrow=n)
rejuvs.ss = 0
x.estimate.ss.fixed = array(rep(0,t*m*length(threshold)*length(partition)),c(m,t,length(threshold),length(partition)))
ESSm <- NULL
Am = NULL
am = NULL


##FUNCTION
resampling <- function(w,n){
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

##Generate Samples
v = rnorm(t)
u = rnorm(t)
x[1] = rnorm(1,0,sigma^2/(1-alpha^2))
y[1] = beta*exp(x[1]/2)*u[1]
for (i in 2:t){
  x[i] = alpha*x[i-1]+sigma*v[i]
  y[i] = beta*exp(x[i]/2)*u[i]
}

##using proposals N(xn-1+yn,1)
##这里有一个问题，为什么和前面的不对应，但是bootstrap filter的时候基本是对应的，可能由随机性引起
for (iterM in length(partition)){#change series
  for (l in 1:m) {
    for (j in 1){#change series
      M = partition[iterM] * n
      X[,1] = rnorm(n)
      W = dnorm(X[,1],0,sigma^2/(1-alpha^2))/dnorm(X[,1])
      ##resampling step
      idx = order(X[,1])
      X = X[idx,]
      W = W[idx]
      Am = NULL
      for (p in 1:(n/M)) {
        Am[p] = sum(W[((p-1)*M+1):(p*M)])
      }
      am = Am/sum(Am)
      ESSx = 1/sum(am^2)
      if (ESSx < threshold[j]*M){
        idx <- resampling(am,n/M)
        for (p in 1:(n/M)) {
          X[((p-1)*M+1):(p*M),] = X[((idx[p]-1)*M+1):(idx[p]*M),]
          W[((p-1)*M+1):(p*M)] = W[((idx[p]-1)*M+1):(idx[p]*M)]
        }
      }
      idx = order(X[,1])
      X = X[idx,]
      W = W[idx]
      for (p in 1:(n/M)){
        wm[((p-1)*M+1):(p*M)] = W[((p-1)*M+1):(p*M)]/sum(W[((p-1)*M+1):(p*M)])
        ESSm[p] = 1/sum(wm[((p-1)*M+1):(p*M)]^2)
        if (ESSm[p] < threshold[j]*n/M){
          idx = resampling(wm[((p-1)*M+1):(p*M)],M)
          X[((p-1)*M+1):(p*M),] = X[idx+(p-1)*M,]
          W[((p-1)*M+1):(p*M)] = rep(sum(W[((p-1)*M+1):(p*M)])/M,M)
        }
      }
      ##这里还有一个问题，就是这是先后顺序的进行resampling，但是实际上这样出来的结果可能会不太ok
      
      
      for (i in 2:t) {
        X[,i]=rnorm(n,X[,i-1]+y[i])
        W = W*dnorm(X[,i],alpha*X[,i-1],sigma^2)/dnorm(X[,i],X[,i-1]+y[i])
        #resampling step
        idx = order(X[,i])
        X = X[idx,]
        W = W[idx]
        for (p in 1:(n/M)) {
          Am[p] = sum(W[((p-1)*M+1):(p*M)])
        }
        am = Am/sum(Am)
        ESSx = 1/sum(am^2)
        if (ESSx < threshold[j]*M){
          idx <- resampling(am,n/M)
          for (p in 1:(n/M)) {
            X[((p-1)*M+1):(p*M),] = X[((idx[p]-1)*M+1):(idx[p]*M),]
            W[((p-1)*M+1):(p*M)] = W[((idx[p]-1)*M+1):(idx[p]*M)]
          }
        }
        idx = order(X[,i])
        X = X[idx,]
        W = W[idx]
        for (p in 1:(n/M)){
          wm[((p-1)*M+1):(p*M)] = W[((p-1)*M+1):(p*M)]/sum(W[((p-1)*M+1):(p*M)])
          ESSm[p] = 1/sum(wm[((p-1)*M+1):(p*M)]^2)
          if (ESSm[p] < threshold[j]*n/M){
            idx = resampling(wm[((p-1)*M+1):(p*M)],M)
            X[((p-1)*M+1):(p*M),] = X[idx+(p-1)*M,]
            W[((p-1)*M+1):(p*M)] = rep(sum(W[((p-1)*M+1):(p*M)])/M,M)#这里好像还要商榷一下，在算总的weight的时候
          }
        }
      }
      #如果每段分析再加和确实会mitigate这个问题
      sm = matrix(rep(0,n/M*t),nrow=n/M)
      for (p in 1:(n/M)) {
        Am[p] = sum(W[((p-1)*M+1):(p*M)])
        wm[((p-1)*M+1):(p*M)] = W[((p-1)*M+1):(p*M)]/sum(W[((p-1)*M+1):(p*M)])
        sm[p,] = colSums(wm[((p-1)*M+1):(p*M)]*X[((p-1)*M+1):(p*M),])
      }
      am = Am/sum(Am)
      for (i in 1:t) {
        x.estimate.ss.fixed[l,i,j,iterM] = sum(am*sm[,i])
      }
    }
  }
}

##Estimates
mse.ss.fixed <- array(rep(0,length(threshold)*t*length(partition)),c(length(threshold),t,length(partition)))
for (j in length(partition)) {#change seires
  for (k in 1){#change series
    for (i in 1:m)
    {
      mse.ss.fixed[k,,j] = mse.ss.fixed[k,,j]+(x-x.estimate.ss.fixed[i,,k,j])^2
    }
  }
}
mse.ss.fixed = mse.ss.fixed/m
mse.sum.ss.fixed <- matrix(rep(0,length(threshold)*length(partition)),nrow=length(threshold))
for (i in 1){#change series
  for (j in length(partition)){#change series
    for (k in 1:t){
      mse.sum.ss.fixed[i,j] = mse.sum.ss.fixed[i,j]+mse.ss.fixed[i,k,j] 
    }
  }
}

for (i in 1:length(mse.sum.ss.fixed[,1])) {
  plot(partition,mse.sum.ss.fixed[i,],ylab="mse",xlab="partition")
}
plot(threshold,mse.sum.ss.fixed[,6])
plot(mse.sum.ss[-1])
for (j in length(partition)) {#change series
  boxplot(x.estimate.ss.fixed[,,1,j])
}
