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

function brconstrained_stabilitytracking_sine(bplusrrange, brratio, smid, freq, numsine)
    data = zeros(length(bplusrrange), 4)
    u0=[0.5]
    tsend = numsine*4*pi
    tspan=(0.0, tsend)
    @threads for bri in eachindex(bplusrrange)
        brvals = brratio_calc(brratio, bplusrrange[bri])
        par = BacPhageSineForcedPar(b = brvals[1], r=brvals[2], per=0.5, amp=0.4, mid=smid)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-4*pi:freq:tsend)
        selectiondata = [sel_sine(par, t) for t in tsend-8*pi:freq:tsend]
        seriesoptimum = optimum(selectiondata)
        optimumdiff = trackoptimum(solseries[1,:], seriesoptimum)
        data[bri, 1] = bplusrrange[bri]
        data[bri, 2] = mean(optimumdiff)
        data[bri, 3] = sumsqddiff(optimumdiff)
        data[bri, 4] = sum(optimumdiff)
    end
    return data
end

function stepwise_indexes(selectiondata)
    optimum0 = []
    optimum1 = []
    for i in eachindex(selectiondata)
        if selectiondata[i] >= 0.0
            append!(optimum1, i)
        else
            append!(optimum0, i)
        end
    end
    return [optimum0, optimum1]
end

function splitoptimumdiff(optimum, optimumdiffdata, selectiondata)
    if optimum == "1"
        indexes = stepwise_indexes(selectiondata)[2]
    elseif optimum =="0"
        indexes = stepwise_indexes(selectiondata)[1]
    else
        error("optimum should be either 1 or 0")
    end
    return [optimumdiffdata[i] for i in indexes]
end

function noise_stabilityprep(bval, rval, smid, freq, tsend, reps)
    CVTLdata = zeros(length(reps))
    meanTLdata = zeros(length(reps))
    sumTLdata = zeros(length(reps))
    for i in 1:length(reps)
        solnoise = bacphage_pert_sol(bval, rval, 0.5, freq, smid, 0.05, 0.0, i, tsend, tsend-1000.0:1.0:tsend)
        seriesoptimum = optimum(solnoise[2])
        optimumdiff = trackoptimum(solnoise[1], seriesoptimum)
        meanTLdata[i] = mean(optimumdiff)
        CVTLdata[i] = std(optimumdiff)/mean(optimumdiff)
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
        timedata[bi] = calc_time_selectionswitch(brange[bi], 0.001, oldsel, newsel, tvals)
    end
    return [brange, timedata]
end
    
function time_selectionswitch_r(rrange, oldsel, newsel, tvals)
    timedata=zeros(length(rrange))
    @threads for ri in eachindex(rrange)
        timedata[ri] = calc_time_selectionswitch(0.001, rrange[ri], oldsel, newsel, tvals)
    end
    return [rrange, timedata]
end