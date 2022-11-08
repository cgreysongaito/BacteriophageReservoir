include("packages.jl")
include("bacteriophage_commoncode.jl")

#problem with model if selection long enough will just go to 1 (all carrier) and stay there regardless of selection changing (equilibrium)


#Trying to assess buffering capacity with and without lysogeny.
#CV
#changing with colour of noise
function CV_calc_colour(rrange)
    mnvals = zeros(length(rrange))
    stdevvals = zeros(length(rrange))
    CVvals = zeros(length(rrange))
    for (ri, rval) in enumerate(rrange)
        sol = bacphage_pert_sol(0.01, [0.5], 1.0, rval, 125, 500.0, 100.0:1.0:500.0)
        mnvals[ri] = mean(sol)
        stdevvals[ri] = std(sol)
        CVvals[ri] = stdevvals[ri]/mnvals[ri]
    end
    return [mnvals, stdevvals, abs.(CVvals)]
end

let
    rrange = 0.0:0.01:0.9
    data = CV_calc_colour(rrange)
    test = figure()
    plot(rrange, data[3])
    return test
end


#changing with periodicity of sine
function CV_calc_per(perrange)
    mnvals = zeros(length(perrange))
    stdevvals = zeros(length(perrange))
    CVvals = zeros(length(perrange))
    for (peri, perval) in enumerate(perrange)
        sol = bacphage_sine_sol(0.01, perval, 0.5, 500.0, 100.0:1.0:500.0)
        mnvals[peri] = mean(sol[1,:])
        stdevvals[peri] = std(sol[1,:])
        CVvals[peri] = stdevvals[peri]/mnvals[peri]
    end
    return [mnvals, stdevvals, abs.(CVvals)]
end

let
    perrange = 0.1:0.1:1.0
    data = CV_calc_per(perrange)
    test = figure()
    plot(perrange, data[3])
    return test
end

#changing with amplitude of sine
function CV_calc_amp(amprange)
    mnvals = zeros(length(amprange))
    stdevvals = zeros(length(amprange))
    CVvals = zeros(length(amprange))
    for (ampi, ampval) in enumerate(amprange)
        sol = bacphage_sine_sol(0.01, 0.2, ampval, 500.0, 100.0:1.0:500.0)
        mnvals[ampi] = mean(sol)
        stdevvals[ampi] = std(sol)
        CVvals[ampi] = stdevvals[ampi]/mnvals[ampi]
    end
    return [mnvals, stdevvals, abs.(CVvals)]
end

let
    amprange = 0.1:0.1:1.0
    data = CV_calc_amp(amprange)
    test = figure()
    plot(amprange, data[3])
    return test
end


#Changing with b

function CV_calc(sol)
    mn = mean(sol)
    stdev = std(sol)
    return [mn, stdev, stdev/mn]
end

#CV of C minus CV of environment  ####***** SOMETHING FEELS WRONG WHEN CALCULATING THE CV OF ENVIRONMENT
function StandCV_calc_per(perrange)
    solCVvals = zeros(length(perrange))
    envirCVvals = zeros(length(perrange))
    StandCVvals = zeros(length(perrange))
    for (peri, perval) in enumerate(perrange)
        sol = bacphage_sine_sol(0.01, perval, 0.5, 500.0, 100.0:1.0:500.0)
        solCVvals[peri] = abs(std(sol[1,:])/mean(sol[1,:]))
        envirCVvals[peri] = abs(std(sol[2,:])/mean(sol[2,:]))
        StandCVvals[peri] = solCVvals[peri]-envirCVvals[peri]
    end
    return [solCVvals,envirCVvals,StandCVvals]
end

StandCV_calc_per(0.1:0.1:1.0)

let
    perrange = 0.1:0.1:1.0
    data = StandCV_calc_per(perrange)
    test = figure()
    # plot(perrange, data[3], color = "blue")
    plot(perrange, data[1], color = "green")
    # plot(perrange, data[2], color = "red")
    return test
end

#changing with amplitude of sine
function StandCV_calc_amp(amprange)
    solCVvals = zeros(length(amprange))
    envirCVvals = zeros(length(amprange))
    StandCVvals = zeros(length(amprange))
    for (ampi, ampval) in enumerate(amprange)
        sol = bacphage_sine_sol(0.01, 0.2, ampval, 500.0, 100.0:1.0:500.0)
        solCVvals[ampi] = abs(std(sol[1,:])/mean(sol[1,:]))
        envirCVvals[ampi] = abs(std(sol[2,:])/mean(sol[2,:]))
        StandCVvals[ampi] = solCVvals[ampi]-envirCVvals[ampi]
    end
    return [solCVvals,envirCVvals,StandCVvals]
end

StandCV_calc_amp(0.1:0.1:1.0)

let
    perrange = 0.1:0.1:1.0
    data = StandCV_calc_amp(perrange)
    test = figure()
    plot(perrange, data[3], color = "blue")
    plot(perrange, data[1], color = "green")
    plot(perrange, data[2], color = "red")
    return test
end

#Speed of response after perturbation

# Tracking of equilibrium (integral area calculation)
#split if function - if below bifurc value return interior if above return 1


#sine
function tracking_sine_per(perrange)
    track = zeros(length(perrange))
    for (peri, perval) in enumerate(perrange)
        sol = bacphage_sine_sol(0.01, perval, 0.5, 500.0, 100.0:1.0:500.0)
        trackpoint = zeros(length(sol))
        for i in 1:length(sol)
            equilval = equil(sol[2,i])
            trackpoint[i] = abs(sol[1,i] - equilval)
        end
        track[peri] = mean(trackpoint)
    end
    return track
end

tracking_sine_per(0.1:0.1:0.9) #longer periods obviously mean better ability of the system to track the equilibrium

#temporal cross-correlation of equilibrium state with actual dynamics
#testing of crosscor
let
    lrange = 1:1:50
    soltest = bacphage_sine_sol(0.01, 0.2, 0.5, 500.0, 100.0:1.0:500.0)
    equilvals = zeros(length(soltest))
        for i in 1:length(soltest)
            equilvals[i] = equil(soltest[2,i])
        end
    cordata = crosscor(soltest[1,:], equilvals, lrange)
    # test = figure()
    # subplot(2,1,1)
    # plot(soltest.t, equilvals)
    # plot(soltest.t, soltest[1,:])
    # subplot(2,1,2)
    # plot(lrange, cordata)
    # tight_layout()
    return findmax(cordata)[2]
end

soltest = bacphage_sine_sol(0.01, 0.2, 0.5, 500.0, 100.0:1.0:500.0)
soltest[1,:]
#maybe best thing is to maximise crosscor and what is this lag - make sure lrange is larger than period
function calc_equil(sel)
    len = length(sel)
        equilvals = zeros(len)
        for i in 1:len
            equilvals[i] = equil(sel[i])
        end
        return equilvals
end

function trackingcor_sine_per(perrange, lrange)
    trackcor = zeros(length(perrange))
    for (peri, perval) in enumerate(perrange)
        sol = bacphage_sine_sol(0.01, perval, 0.5, 500.0, 100.0:1.0:500.0)
        splitsolsel = vector_prod(sol)
        equilvals = calc_equil(splitsolsel[2])
        cordata = crosscor(splitsolsel[1], equilvals, lrange)
        trackcor[peri] = findmax(cordata)[2]
    end
    return trackcor
end

trackingcor_sine_per(0.1:0.1:1.8, 1:1:100)
#need to test this on one per value to check how to summarize over multiple per values

function trackingcor_noise_r(rrange, lrange)
    trackcor = zeros(length(rrange))
    for (ri, rval) in enumerate(rrange)
        sol = bacphage_pert_sol(0.01, [0.5], 1.0, rval, 125, 500.0, 100.0:1.0:500.0)
        splitsolsel = vector_prod(sol)
        equilvals = calc_equil(splitsolsel[2])
        cordata = crosscor(splitsolsel[1], equilvals, lrange)
        trackcor[ri] = findmax(cordata)[2]
    end
    return trackcor
end

trackingcor_noise_r(0.1:0.1:0.9, 1:1:100)