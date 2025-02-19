errorsP = zeros(4,5)
errorsG = zeros(4,5)
WerrorsP = zeros(4,5)
WerrorsG = zeros(4,5)
for N in 2:5
    for j in 1:5
        @rput N
        R"x = generateChain(list(type  = 'AR', rho = 0.5, M = 10^N))$x"
        R"r = autocov(x)"
        R"dhat = tune_delta(x,5)$delta*0.8"
        @rget r
        @rget dhat
        m = momentLSmod(r, dhat, [0.3],[0.0], 3e-16)
        supp = m[1]
        weight = m[2]
        @rput supp
        @rput weight
        R"m = SR1(r,dhat)"
        R"errorp = L2diff_L2Moment(r,supp, weight)"
        R"errorg = L2diff_L2Moment(r, m$support, m$weights)"
        R"errG = wasserstein1d(m$support, c(0.5), p = 1, m$weights, c(r[1]))"
        R"errP = wasserstein1d(supp, c(0.5), p = 1, weight, c(r[1]))"
        @rget errorp
        @rget errorg
        @rget errP
        @rget errG
        errorsP[N-1, j] = errorp
        errorsG[N-1,j] = errorg
        WerrorsP[N-1,j] = errP
        WerrorsG[N-1,j] = errG
    end
end
