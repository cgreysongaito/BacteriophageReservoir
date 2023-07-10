#Eigenvalues integral balancing (gut responses to environmental variation)
#eigen integral bifurcation #C=1 #mid
function bifurcintegral_eigen1_mid(midrange, tsend)
    data = zeros(length(midrange), 3)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for midi in eachindex(midrange)
        par = BacPhageSineForcedPar(mid=midrange[midi])
        data[midi, 1] = calc_integral_eigen1(par)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000.0:1.0:tsend)
        data[midi, 2] = maximum(solseries)
        data[midi, 3] = minimum(solseries)
    end
    return data
end

# lowest C equilibrium analysis for mean transitory load
function meanTL_lowestC(bplusrrange, brratio, slow)
    data = zeros(length(bplusrrange), 2)
    @threads for bri in eachindex(bplusrrange)
        brvals = brratio_calc(brratio, bplusrrange[bri])
        par = BacPhagePar(b = brvals[1], r=brvals[2])
        data[bri, 1] = bplusrrange[bri]
        data[bri, 2] = interior_equil(slow, par)
    end
    return data
end

#Different periodicities and amplitudes
#test of bifurcmid with different periodicity (because model without bac is dependent on periodicity) - ANSWER - periodicity does not affect "bifurcation" point
function bifurcmid_per(midrange, perval, tsend)
    data = zeros(length(midrange), 3)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for midi in eachindex(midrange)
        par = BacPhageSineForcedPar(per=perval, mid=midrange[midi])
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000:1.0:tsend)
        data[midi, 1] = midrange[midi]
        data[midi, 2] = maximum(solseries)
        data[midi, 3] = minimum(solseries)
    end
    return data
end