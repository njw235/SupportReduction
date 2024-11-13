library(momentLS)
load("data/MC_chains.Rdata")

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
  measure = julia_call("momentLSmod", r, dhat, m$support, m$weights, 1e-8)
  
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

avar_mh = c()
avar_mh2 = c()
er_mh = c()
er_mh2 - c()

for(j in 1:ncol(ch_mh$x)){
  r = autocov(ch_mh$x[,j])
  dhat = tune_delta(ch_mh$x[,j],5)$delta*0.8
  
  m <- SR1(r = r, delta = dhat)
  avar_mh[j] <- avar(m)
  measure = julia_call("momentLSmod", r,  dhat, m$support, m$weights, 1e-8)
  
  support = vector(mode = "numeric", length = length(measure[[1]]) - 1)
  weight = vector(mode = "numeric", length = length(measure[[2]]) -1)
  
  for(i in 1:length(support)){
    support[i] = measure[[1]][[i]]
    weight[i] = measure[[2]][[i]]
  }
  avar_mh2[j] = asympVariance(weight, support)
  er_mh[j] = L2diff_L2Moment(r, m$support, m$weights)
  er_mh2[j] = L2diff_L2Moment(r, support, weight)
}

avar_pg = c()
avar_pg2 = c()
er_pg = c()
er_pg2 - c()
for(j in 1:ncol(ch_pg$x)){
  r = autocov(ch_pg$x[,j])
  dhat = tune_delta(ch_pg$x[,j],5)$delta*0.8
  
  m <- SR1(r = r, delta = dhat)
  avar_pg[j] <- avar(m)
  measure = julia_call("momentLSmod", r,  dhat, m$support, m$weights, 1e-8)
  
  support = vector(mode = "numeric", length = length(measure[[1]]) - 1)
  weight = vector(mode = "numeric", length = length(measure[[2]]) -1)
  
  for(i in 1:length(support)){
    support[i] = measure[[1]][[i]]
    weight[i] = measure[[2]][[i]]
  }
  avar_pg2[j] = asympVariance(weight, support)
  er_pg[j] = L2diff_L2Moment(r, m$support, m$weights)
  er_pg2[j] = L2diff_L2Moment(r, support, weight)
}



avar_var = c()
avar_var2 = c()
er_var = c()
er_var2 - c()
for(j in 1:ncol(ch_var$x)){
  r = autocov(ch_var$x[,j])
  dhat = tune_delta(ch_var$x[,j],5)$delta*0.8
  
  m <- SR1(r = r, delta = dhat)
  avar_var[j] <- avar(m)
  measure = julia_call("momentLSmod", r, dhat, m$support, m$weights, 1e-8)
  
  support = vector(mode = "numeric", length = length(measure[[1]]) - 1)
  weight = vector(mode = "numeric", length = length(measure[[2]]) -1)
  
  for(i in 1:length(support)){
    support[i] = measure[[1]][[i]]
    weight[i] = measure[[2]][[i]]
  }
  avar_var2[j] = asympVariance(weight, support)
  er_var[j] = L2diff_L2Moment(r, m$support, m$weights)
  er_var2[j] = L2diff_L2Moment(r, support, weight)
}

