#Sine wave
#"bifurcations" changing mid
function bifurcmid(midrange, tsend)
    data = zeros(length(midrange), 3)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for midi in eachindex(midrange)
        par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=midrange[midi])
        par.mid = midrange[midi]
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000:1.0:tsend)
        data[midi, 1] = midrange[midi]
        data[midi, 2] = maximum(solseries)
        data[midi, 3] = minimum(solseries)
    end
    return data
end

#Noise

function noise_mid_reps(μ, corr, reps, tend)
    maxC = zeros(reps)
    minC = zeros(reps)
    for i in 1:reps
        solseries = bacphage_pert_sol(0.001, 0.5, 1.0, μ, 0.05, corr, i, tend, tend-1000.0:1.0:tend)
        maxC[i] = maximum(solseries[1])
        minC[i] = minimum(solseries[1])
    end
    return [mean(maxC), mean(minC)]
end


#White noise

function bifurc_white_mid(midrange, reps, tend)
    data = zeros(length(midrange), 3)
    @threads for midi in eachindex(midrange)
        maxminmeans = noise_mid_reps(midrange[midi], 0.0, reps, tend)
        data[midi, 1] = midrange[midi]
        data[midi, 2] = maxminmeans[1]
        data[midi, 3] = maxminmeans[2]
    end
    return data
end

#Red noise

function bifurc_red_mid(corrrange, mid, reps, tend)
    data = zeros(length(corrrange), 3)
    @threads for corri in eachindex(corrrange)
        maxminmeans = noise_mid_reps(mid, corrrange[corri], reps, tend)
        data[corri, 1] = corrrange[corri]
        data[corri, 2] = maxminmeans[1]
        data[corri, 3] = maxminmeans[2]
    end
    return data
end