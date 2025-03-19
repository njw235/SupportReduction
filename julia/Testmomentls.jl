include("/julia/momentLS.jl")
R"load('data/MC_chains.Rdata')"

errorsP = zeros(4,10)
errorsG = zeros(4,10)
AerrorsP = zeros(4,10)
AerrorsG = zeros(4,10)
R"set.seed(1234)"
for N in 2:5
    for j in 1:10
        @rput N
        R"x = generateChain(list(type  = 'AR', rho = 0.9, M = 6^N))$x"
        R"r = autocov(x)"
        R"dhat = tune_delta(x,5)$delta*0.8"
        @rget r
        @rget dhat
        m = momentLSmod(r, dhat, [0.0],[0.0], 3e-16)
        supp = m[1]
        weight = m[2]
        @rput supp
        @rput weight
        R"m = SR1(r,dhat)"
        R"errorp = L2diff_L2Moment(r,supp, weight)"
        R"errorg = L2diff_L2Moment(r, m$support, m$weights)"
        R"errG = (asympVariance(m$weights, m$support) - 100)^2"
        R"errP = (asympVariance(weight, supp) - 100)^2"
        @rget errorp
        @rget errorg
        @rget errP
        @rget errG
        errorsP[N-1, j] = errorp
        errorsG[N-1,j] = errorg
        AerrorsP[N-1,j] = errP
        AerrorsG[N-1,j] = errG
    end
end

println(errorsP)
println(errorsG)
println(AerrorsP)
println(AerrorsG)

