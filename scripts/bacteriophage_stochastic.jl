include("packages.jl")
include("bacteriophage_commoncode.jl")


#white noise should be exactly the same as sine (outcome depends on combination of r, b, and mid/mean)
#red noise effects will be dependent on "integral" of the red noise (above the bifurcation point) in combination with r, b, and mid/mean
#BUT unless you have really strong r and b - still highly unlikely to go to fixation


#Prediction for white and red noise for tracking attractor CV

#using PeriodicCallback

# let 
#     tsend = 10000.0
#     tvals = 0.0:1.0:tsend
#     data = bacphage_pert_sol(0.001, 0.5, 1.0, -0.007, 0.05, 0.9, 132, tsend, tvals)
#     test = figure()
#     plot(tvals, data[1])
#     plot(tvals, data[2])
#     return test
# end

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


#white noise

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

# let 
#     par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=0.0)
#     midrange = -0.01:0.0001:0.001
#     data = bifurc_white_mid(midrange, 6, 100000.0)
#     test = figure()
#     plot(data[:, 1], data[:, 2])
#     plot(data[:, 1], data[:, 3])
#     vlines(bifurc(par), 0.0, 1.0) #bifurc value or r+b=-mid
#     return test
# end


#Red noise maybe just a more probabilistic way - testing probability of going to fixation - at three mid points (above bifurc, at bifurc, and below bifurc)

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

# let 
#     par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=0.0)
#     corrrange = 0.0:0.1:0.9
#     dataabove = bifurc_red_mid(corrrange, 0.0, 6, 100000.0)
#     datacentred = bifurc_red_mid(corrrange, -0.002, 6, 100000.0)
#     databelow = bifurc_red_mid(corrrange, -0.004, 6, 100000.0)
#     test = figure()
#     subplot(3,1,1)
#     plot(dataabove[:, 1], dataabove[:, 2])
#     plot(dataabove[:, 1], dataabove[:, 3])
#     ylim(0,1.1)
#     subplot(3,1,2)
#     plot(datacentred[:, 1], datacentred[:, 2])
#     plot(datacentred[:, 1], datacentred[:, 3])
#     ylim(0,1.1)
#     subplot(3,1,3)
#     plot(databelow[:, 1], databelow[:, 2])
#     plot(databelow[:, 1], databelow[:, 3])
#     ylim(0,1.1)
#     tight_layout()
#     return test
# end

# let 
#     tsend = 100000.0
#     tvals = tsend-1000.0:1.0:tsend
#     data = bacphage_pert_sol(0.001, 0.5, 1.0, -0.004, 0.05, 0.8, 4, tsend, tvals)
#     println(maximum(data[1]))
#     println(minimum(data[1]))
#     test = figure()
#     plot(tvals, data[1])
#     plot(tvals, data[2])
#     return test
# end

#using SDE - where stochasticity is just on selection