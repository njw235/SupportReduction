using Clarabel
using SumOfSquares
using DynamicPolynomials
using RCall
using Random
using Distributions
using LinearAlgebra
using LinearSolve
using Plots
using StatsBase

Empirical = function(x)
    
    n = length(x)
    supp = [0:1:findmax(x)[1];]
    count = proportions(x)
    
    return(Dict("support" => supp, "count" => count))

end


