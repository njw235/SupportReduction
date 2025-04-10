library(JuliaCall)

julia = julia_setup()

julia_source("julia/MomentLS.jl")

gridfreeSR <- function(r, delta_tilde, tol, type = "sub"){
  if(type == "sub"){
    measure = julia_call("momentLSmod", r, delta_tilde, 1e-7)
  }
  else{
    measure = julia_call("momentLS", -1+ delta_tilde, 1-delta_tilde, r, 1e-7)
  }

  support = vector(mode = "numeric", length = length(measure[[1]]) - 1)
  weight = vector(mode = "numeric", length = length(measure[[2]]) -1)
  
  for(i in 1:length(support)){
    support[i] = measure[[1]][[i]]
    weight[i] = measure[[2]][[i]]
  }


  return(c(support, weight))
}
