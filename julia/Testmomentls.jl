R"load('data/MC_chains.Rdata')"

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
        m = momentLSmod(r, dhat, [0.0],[0.0], 3e-16)
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



errorbG = zeros(9,2)
errormhG = zeros(5,2)
errorpgG = zeros(6,2)
errorvarG = zeros(4,2)

errorbP = zeros(9,2)
errormhP = zeros(5,2)
errorpgP = zeros(6,2)
errorvarP = zeros(4,2)

for i in 1:9
    @rput i
    R"r = autocov(ch_blasso$x[,i])"
    R"dhat = tune_delta(ch_blasso$x[,i],5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0], [0.0], 3e-16)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight
    R"m = SR1(r, dhat)"
    R"err1 = L2diff_L2Moment(r, m$support, m$weights)"
    R"err2 = L2diff_L2Moment(r, supp, weight)"
    R"err3 = (asympVariance(weight, supp) - ch_blasso$varTruth[i,i])^2$"
    R"err4 = (asympVariance(m$weights, m$support) - ch_blasso$varTruth[i,i])^2$"
    @rget err1
    @rget err2
    @rget err3
    @rget err4
    errorbG[i,1] = err1
    errorbP[i,1] = err2
    errorbP[i,2] = err3
    errorbG[i,2] = err4
end

for i in 1:5
    @rput i
    R"r = autocov(ch_mh$x[,i])"
    R"dhat = tune_delta(ch_mh$x[,i],5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0], [0.0], 3e-16)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight
    R"m = SR1(r, dhat)"
    R"err1 = L2diff_L2Moment(r, m$support, m$weights)"
    R"err2 = L2diff_L2Moment(r, supp, weight)"
    R"err3 = (asympVariance(weight, supp) - ch_mh$varTruth[i,i])^2$"
    R"err4 = (asympVariance(m$weights, m$support) - ch_mh$varTruth[i,i])^2$"
    @rget err1
    @rget err2
    @rget err3
    @rget err4
    errormhG[i,1] = err1
    errormhP[i,1] = err2
    errormhP[i,2] = err3
    errormhG[i,2] = err4
end

for i in 1:6
    @rput i
    R"r = autocov(ch_pg$x[,i])"
    R"dhat = tune_delta(ch_pg$x[,i],5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0], [0.0], 3e-16)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight
    R"m = SR1(r, dhat)"
    R"err1 = L2diff_L2Moment(r, m$support, m$weights)"
    R"err2 = L2diff_L2Moment(r, supp, weight)"
    R"err3 = (asympVariance(weight, supp) - ch_pg$varTruth[i,i])^2$"
    R"err4 = (asympVariance(m$weights, m$support) - ch_pg$varTruth[i,i])^2$"
    @rget err1
    @rget err2
    @rget err3
    @rget err4
    errorpgG[i,1] = err1
    errorpgP[i,1] = err2
    errorpgP[i,2] = err3
    errorpgG[i,2] = err4
end

for i in 1:4
    @rput i
    R"r = autocov(ch_var$x[,i])"
    R"dhat = tune_delta(ch_var$x[,i],5)$delta*0.8"
    @rget r
    @rget dhat
    m = momentLSmod(r, dhat, [0.0], [0.0], 3e-16)
    supp = m[1]
    weight = m[2]
    @rput supp
    @rput weight
    R"m = SR1(r, dhat)"
    R"err1 = L2diff_L2Moment(r, m$support, m$weights)"
    R"err2 = L2diff_L2Moment(r, supp, weight)"
    R"err3 = (asympVariance(weight, supp) - ch_var$varTruth[i,i])^2$"
    R"err4 = (asympVariance(m$weights, m$support) - ch_var$varTruth[i,i])^2$"
    @rget err1
    @rget err2
    @rget err3
    @rget err4
    errorvarG[i,1] = err1
    errorvarP[i,1] = err2
    errorvarP[i,2] = err3
    errorvarG[i,2] = err4
end