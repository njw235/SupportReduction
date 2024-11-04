library(momentLS)
load("MC_chains.Rdata")

# univariate MCSE estimation for chains in ch_blasso
avar_blasso <- c()
avar_blasso2 <- c()
error1 = c()
error2 = c()
for(j in 1:ncol(ch_blasso$x)){
  r <- autocov(ch_blasso$x[,j])
  dhat<- tune_delta(ch_blasso$x[,j], 5)$delta*0.8
  m <- SR1(r = r, delta = dhat)
  avar_blasso[j] <- avar(m)
  measure = julia_call("momentLSmod", r, dhat, 1e-3)
  
  support = vector(mode = "numeric", length = length(measure[[1]]) - 1)
  weight = vector(mode = "numeric", length = length(measure[[2]]) -1)
  
  for(i in 1:length(support)){
    support[i] = measure[[1]][[i]]
    weight[i] = measure[[2]][[i]]
  }
  avar_blasso2[j] = asympVariance(weight, support)
  error1[j] = L2diff_L2Moment(r, m$support, m$weights)
  error2[j] = L2diff_L2Moment(r, support, weight)
}
