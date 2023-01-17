##Resampling Scheme
#residual resampling
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

uniform_spacing = function(N){
  a = -log(runif(N+1))
  b = cumsum(a)
  c = sort(b,decreasing = F)
  return(c/b[length(b)])
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

multinomial = function(w){
  return(inv_cdf(uniform_spacing(length(w)),w))
}

# stratified resampling
stratified = function(w){
  m = length(w)
  su = (runif(m)+0:(m-1))/m
  return(inv_cdf(su,w))
}

# systematic resampling
systematic = function(w){
  m = length(w)
  su = (rep(runif(1),m)+0:(m-1))/m
  return(inv_cdf(su,w))
}