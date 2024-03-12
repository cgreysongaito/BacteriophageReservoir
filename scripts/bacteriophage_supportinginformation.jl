include("bacteriophage_commoncode.jl")

#Parameter ranges
include("bacteriophage_stability.jl")

function range_parameter_br(brratio, bplusrrange)
    brange = zeros(length(bplusrrange))
    rrange = zeros(length(bplusrrange))
    @threads for bri in eachindex(bplusrrange)
        br = brratio_calc(brratio, bplusrrange[bri])
        brange[bri] = br[1]
        rrange[bri] = br[2]
    end
    minb = minimum(brange)
    maxb = maximum(brange)
    minr = minimum(rrange)
    maxr = maximum(rrange)
    return hcat(minb, maxb, minr, maxr)
end

range_parameter_br(0.1, 0.00001:0.0001:0.001)
range_parameter_br(1.0, 0.00001:0.0001:0.001)
range_parameter_br(10.0, 0.00001:0.0001:0.001)

0.05*sin(π/2)-0.0005
-0.05*sin(π/2)-0.0005

function sel_noise_range(reps)
    mins = zeros(reps)
    maxs = zeros(reps)
    @threads for i in 1:reps
        noisedata = noise_creation(-0.0005, 0.005, 0.0, 1000, i)
        mins[i] = minimum(noisedata)
        maxs[i] = maximum(noisedata)
    end
    return [minimum(mins), maximum(maxs)]
end

sel_noise_range(100000)

#Gut responses to environmental variation
#Sine wave
function bifurcmid(midrange, tsend)
    data = zeros(length(midrange), 3)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for midi in eachindex(midrange)
        par = BacPhageSineForcedPar(mid=midrange[midi])
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000:1.0:tsend)
        data[midi, 1] = midrange[midi]
        data[midi, 2] = maximum(solseries)
        data[midi, 3] = minimum(solseries)
    end
    return data
end

#White Noise

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