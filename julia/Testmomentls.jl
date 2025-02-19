errors = zeros(4,5)
for N in 2:5
    for j in 1:5
        @rput N
        R"x = generateChain(list(type  = 'AR', rho = 0.5, M = N))$x"
        R"r = autocov(x)"
        R"dhat = tune_delta(x,5)$delta*0.8"
        @rget r
        @rget dhat
        m = momentLSmod(r, dhat, [0.3],[0.0], 1e-8)
        supp = m[1]
        weight = m[2]
        @rput supp
        @rput weight
        R"m = SR1(r,dhat)"
        R"errordiff = L2diff_L2Moment(r, m$support, m$weights) - L2diff_L2Moment(r,supp, weight)"
        @rget errordiff
        errors[N-1, j] = errordiff
    end
end

errors
