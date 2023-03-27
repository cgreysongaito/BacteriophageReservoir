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

#Tracking equilibrium curve
function attractordata(selectiondata, par)
    attractor = zeros(length(selectiondata))
    for i in 1:length(selectiondata)
        attractor[i] = stableequil(selectiondata[i], par)
    end
    return attractor
end 

function setupbparam(brange)
    for bi in 1:length(brange)
        @eval $(Symbol(:b, bi)) = $brange[$bi]
    end
end

function trackattractor_INT(solutiondata, attractordata)
    INTdiff = zeros(length(solutiondata))
    for i in 1:length(solutiondata)
        INTdiff[i] = solutiondata[i] - attractordata[i]
    end
    return log10(abs(sum(INTdiff)))
end

function trackattractor_INT2(solutiondata, attractordata)
    INTdiff = zeros(length(solutiondata))
    for i in 1:length(solutiondata)
        INTdiff[i] = abs(solutiondata[i] - attractordata[i])
    end
    return sum(INTdiff)
end

function trackattractor_INT3(solutiondata, attractordata)
    return log10(abs(sum(solutiondata)-sum(attractordata)))
end

function trackingcor_INT_b(brange, rval, freq, tsend)
    data = zeros(length(brange), 2)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for bi in eachindex(brange)
        par = BacPhageSineForcedPar(b = brange[bi], r=rval, per=0.5, amp=0.4, mid=-0.002)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000.0:freq:tsend)
        selection = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
        seriesattractor = attractordata(selection, par)
        INT = trackattractor_INT(solseries[1,:], seriesattractor)
        data[bi, 1] = brange[bi]
        data[bi, 2] = INT
    end
    return data
end

function brconstrained_tracking(rplusbrange, rbratio, smid, freq, tsend)
    data = zeros(length(rplusbrange), 2)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for rbi in eachindex(rplusbrange)
        if rbratio == Inf
            rval = rplusbrange[rbi]
            bval = 0.0
        else
            rval = rplusbrange[rbi]/(1+(1/rbratio))
            bval = rplusbrange[rbi] - rval
        end
        par = BacPhageSineForcedPar(b = bval, r=rval, per=0.5, amp=0.4, mid=smid)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000.0:freq:tsend)
        selection = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
        seriesattractor = attractordata(selection, par)
        INT = trackattractor_INT3(solseries[1,:], seriesattractor)
        data[rbi, 1] = rplusbrange[rbi]
        data[rbi, 2] = INT
    end
    return data
end