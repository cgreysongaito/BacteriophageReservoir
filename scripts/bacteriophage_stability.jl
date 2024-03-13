#Tracking moving optimum (0 or 1 frequency)
function optimum(selectiondata)
    optimumdata = zeros(length(selectiondata))
    for si in eachindex(selectiondata)
        if selectiondata[si] > 0.0
            optimumdata[si] = 1.0
        elseif selectiondata[si] < 0.0
            optimumdata[si] = 0.0
        else
            optimumdata[si] = NaN
        end
    end
    return optimumdata
end


function setupbparam(brange)
    for bi in 1:length(brange)
        @eval $(Symbol(:b, bi)) = $brange[$bi]
    end
end

function trackoptimum(solutiondata, optimumdata)
    optimumdiff = zeros(length(solutiondata))
    for i in 1:length(solutiondata)
        optimumdiff[i] = abs(solutiondata[i] - optimumdata[i])
    end
    return optimumdiff
end

function sumsqddiff(data)
    meandata = mean(data)
    sqddiff = zeros(length(data))
    for i in eachindex(data)
        sqddiff[i] = (data[i]-meandata)^2
    end
    return sum(sqddiff)
end

function brratio_calc(brratio, bplusr)
    rval = bplusr/(1+brratio)
    bval = bplusr - rval
    return [bval, rval]
end

function sel_rel_conj(brratio, bplusrrange, avabssel)
    alphadata = zeros(length(bplusrrange))
    @threads for i in eachindex(bplusrrange)
        rval = brratio_calc(brratio, bplusrrange[i])[2]
        alphadata[i] = avabssel/rval
    end
    return hcat(bplusrrange, alphadata)
end

function bifurc_consrelsel(brratio, bplusrrange, αval)
    rsbdata=zeros(length(bplusrrange))
    @threads for i in eachindex(bplusrrange)
        brval = brratio_calc(brratio, bplusrrange[i])
        rsbdata[i] = brval[2]*(1+αval) + brval[1]
    end
    return hcat(bplusrrange, rsbdata)
end


function brconstrained_stabilitytracking_sine(bplusrrange, brratio, smid, freq, numsine)
    data = zeros(length(bplusrrange), 4)
    u0=[0.5]
    tsend = numsine*4*pi
    tspan=(0.0, tsend)
    @threads for bri in eachindex(bplusrrange)
        brvals = brratio_calc(brratio, bplusrrange[bri])
        par = BacPhageSineForcedPar(b = brvals[1], r=brvals[2], per=0.5, amp=0.05, mid=smid)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        # sol = solve(prob, RadauIIA5())
        sol = solve(prob, Rodas4P())
        solseries = sol(tsend-4*pi:freq:tsend)
        selectiondata = [sel_sine(par, t) for t in tsend-4*pi:freq:tsend]
        seriesoptimum = optimum(selectiondata)
        optimumdiff = trackoptimum(solseries[1,:], seriesoptimum)
        data[bri, 1] = bplusrrange[bri]
        data[bri, 2] = mean(optimumdiff)
        data[bri, 3] = sumsqddiff(optimumdiff)
        data[bri, 4] = sum(optimumdiff)
    end
    return data
end

function stepwise_indexes(optimumdata)
    optimum0 = []
    optimum1 = []
    for i in eachindex(optimumdata)
        if optimumdata[i] > 0.0
            append!(optimum1, i)
        else
            append!(optimum0, i)
        end
    end
    return [optimum0, optimum1]
end

function splitoptimumdiff(optimum, optimumdiffdata, optimumdata)
    if optimum == "1"
        indexes = stepwise_indexes(optimumdata)[2]
    elseif optimum =="0"
        indexes = stepwise_indexes(optimumdata)[1]
    else
        error("optimum should be either 1 or 0")
    end
    return [optimumdiffdata[i] for i in indexes]
end

function brconstrained_stabilitytracking_sine_splitoptimum(bplusrrange, brratio, smid, freq, numsine)
    data = zeros(length(bplusrrange), 5)
    u0=[0.5]
    tsend = numsine*4*pi
    tspan=(0.0, tsend)
    @threads for bri in eachindex(bplusrrange)
        brvals = brratio_calc(brratio, bplusrrange[bri])
        par = BacPhageSineForcedPar(b = brvals[1], r=brvals[2], per=0.5, amp=0.05, mid=smid)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        # sol = solve(prob, RadauIIA5())
        sol = solve(prob, Rodas4P())
        solseries = sol(tsend-4*pi:freq:tsend)
        selectiondata = [sel_sine(par, t) for t in tsend-4*pi:freq:tsend]
        seriesoptimum = optimum(selectiondata)
        optimumdiff = trackoptimum(solseries[1,:], seriesoptimum)
        optimumdiff0 = splitoptimumdiff("0", optimumdiff, seriesoptimum)
        optimumdiff1 = splitoptimumdiff("1", optimumdiff, seriesoptimum)
        data[bri, 1] = bplusrrange[bri]
        data[bri, 2] = mean(optimumdiff0)
        data[bri, 3] = sumsqddiff(optimumdiff0)
        data[bri, 4] = mean(optimumdiff1)
        data[bri, 5] = sumsqddiff(optimumdiff1)
    end
    return data
end

function noise_stabilityprep(bval, rval, smid, freq, tsend, reps)
    CVTLdata = zeros(length(reps))
    meanTLdata = zeros(length(reps))
    sumTLdata = zeros(length(reps))
    for i in 1:length(reps)
        solnoise = bacphage_pert_sol(bval, rval, 0.5, freq, smid, 0.005, 0.0, i, tsend, tsend-1000.0:1.0:tsend)
        seriesoptimum = optimum(solnoise[2])
        optimumdiff = trackoptimum(solnoise[1], seriesoptimum)
        meanTLdata[i] = mean(optimumdiff)
        # CVTLdata[i] = std(optimumdiff)/mean(optimumdiff)
        CVTLdata[i] = std(optimumdiff)
        sumTLdata[i] = sum(optimumdiff)
    end
    return [mean(meanTLdata), mean(CVTLdata), mean(sumTLdata)]
end

function brconstrained_stabilitytracking_noise(bplusrrange, brratio, smid, freq, tsend, reps)
    data = zeros(length(bplusrrange), 4)
    @threads for bri in eachindex(bplusrrange)
        brvals = brratio_calc(brratio, bplusrrange[bri])
        TLstats = noise_stabilityprep(brvals[1], brvals[2], smid, freq, tsend, reps)
        data[bri, 1] = bplusrrange[bri]
        data[bri, 2] = TLstats[1]
        data[bri, 3] = TLstats[2]
        data[bri, 4] = TLstats[3]
    end
    return data
end

#split optimum for noise
function noise_stabilityprep_splitoptimum(bval, rval, smid, freq, tsend, reps)
    CVTLdata0 = zeros(length(reps))
    CVTLdata1 = zeros(length(reps))
    meanTLdata0 = zeros(length(reps))
    meanTLdata1 = zeros(length(reps))
    for i in 1:length(reps)
        solnoise = bacphage_pert_sol(bval, rval, 0.5, freq, smid, 0.005, 0.0, i, tsend, tsend-1000.0:1.0:tsend)
        seriesoptimum = optimum(solnoise[2])
        optimumdiff = trackoptimum(solnoise[1], seriesoptimum)
        optimumdiff0 = splitoptimumdiff("0", optimumdiff, seriesoptimum)
        optimumdiff1 = splitoptimumdiff("1", optimumdiff, seriesoptimum)
        meanTLdata0[i] = mean(optimumdiff0)
        meanTLdata1[i] = mean(optimumdiff1)
        CVTLdata0[i] = std(optimumdiff0)
        CVTLdata1[i] = std(optimumdiff1)
    end
    return [mean(meanTLdata0), mean(CVTLdata0), mean(meanTLdata1), mean(CVTLdata1)]
end

function brconstrained_stabilitytracking_noise_splitoptimum(bplusrrange, brratio, smid, freq, tsend, reps)
    data = zeros(length(bplusrrange), 5)
    @threads for bri in eachindex(bplusrrange)
        brvals = brratio_calc(brratio, bplusrrange[bri])
        TLstats = noise_stabilityprep_splitoptimum(brvals[1], brvals[2], smid, freq, tsend, reps)
        data[bri, 1] = bplusrrange[bri]
        data[bri, 2] = TLstats[1]
        data[bri, 3] = TLstats[2]
        data[bri, 4] = TLstats[3]
        data[bri, 5] = TLstats[4]
    end
    return data
end

#time taken to reach new equilibrium after switch in selection for different values of b
function calc_time_selectionswitch(bval, rval, oldsel, newsel, tvals)
    solseries = selection_switch_bacphage(bval, rval, oldsel, newsel, tvals)
    newequil = stableequil(newsel, BacPhagePar(b=bval, r=rval))
    time = []
    for ti in eachindex(solseries.t)
        if isapprox(solseries.u[ti][1], newequil, atol=0.1)
            time = solseries.t[ti]
            break
        end
    end
    return time
end
    
function time_selectionswitch_b(brange, oldsel, newsel, tvals)
    timedata=zeros(length(brange))
    @threads for bi in eachindex(brange)
        timedata[bi] = calc_time_selectionswitch(brange[bi], 0.0001, oldsel, newsel, tvals)
    end
    return [brange, timedata]
end
    
function time_selectionswitch_r(rrange, oldsel, newsel, tvals)
    timedata=zeros(length(rrange))
    @threads for ri in eachindex(rrange)
        timedata[ri] = calc_time_selectionswitch(0.0001, rrange[ri], oldsel, newsel, tvals)
    end
    return [rrange, timedata]
end