include("packages.jl")
include("bacteriophage_commoncode.jl")

#Proof 1 bacteriophage does most of the work (but depends)

#Examining HGT across C with r, b, s(mid) stacked on top of each other (geometric intuition)

function fillbetween_setup(data_low, data_high)
    final_data = zeros(length(data_low))
    for i in 1:length(data_low)
        final_data[i] = data_low[i] + data_high[i]
    end
    return final_data
end

let 
    par = BacPhagePar(s = -0.002)
    Crange = 0.0:0.01:1.0
    conj = [conjugation(C, par.r) for C in Crange]
    bac = [lysogeny(C, par.b) for C in Crange]
    sel = [selection(C, par.s) for C in Crange]
    conj_bac = fillbetween_setup(conj, bac)
    conj_bac_sel = fillbetween_setup(conj_bac, sel)
    test = figure()
    # plot(Crange, sel)
    fill_between(Crange, sel, color="#73D055FF")
    fill_between(Crange, conj, color="#440154FF")
    fill_between(Crange, conj, conj_bac, color="#404788FF")
    title("s=-0.002")
    # fill_between(Crange, conj_bac, conj_bac_sel, color="#73D055FF")
    return test
end

let 
    par = BacPhagePar(s = -0.004)
    Crange = 0.0:0.01:1.0
    conj = [conjugation(C, par.r) for C in Crange]
    bac = [lysogeny(C, par.b) for C in Crange]
    sel = [selection(C, par.s) for C in Crange]
    conj_bac = fillbetween_setup(conj, bac)
    conj_bac_sel = fillbetween_setup(conj_bac, sel)
    test = figure()
    # plot(Crange, sel)
    fill_between(Crange, sel, color="#73D055FF")
    fill_between(Crange, conj, color="#440154FF")
    fill_between(Crange, conj, conj_bac, color="#404788FF")
    # fill_between(Crange, conj_bac, conj_bac_sel, color="#73D055FF")
    title("s=-0.004")
    return test
end

let 
    par = BacPhagePar(s = -0.001)
    Crange = 0.0:0.01:1.0
    conj = [conjugation(C, par.r) for C in Crange]
    bac = [lysogeny(C, par.b) for C in Crange]
    sel = [selection(C, par.s) for C in Crange]
    conj_bac = fillbetween_setup(conj, bac)
    conj_bac_sel = fillbetween_setup(conj_bac, sel)
    test = figure()
    # plot(Crange, sel)
    fill_between(Crange, sel, color="#73D055FF")
    fill_between(Crange, conj, color="#440154FF")
    fill_between(Crange, conj, conj_bac, color="#404788FF")
    # fill_between(Crange, conj_bac, conj_bac_sel, color="#73D055FF")
    title("s=-0.001")
    return test
end

##exploring contributions of conjugation versus selection versus lysogeny after a change to positive selection
#May be able to get same intuition from the equations - but time series still useful


let 
    u0 = [0.5]
    tend = 5000.0
    tspan=(0.0, tend)
    par = BacPhagePar(b = 0.001, r=0.001, s=-0.004)
    condition(u,t,integrator) = t > 2000.0 
    function changesel!(integrator)
        integrator.p.s=-0.001
    end

    cb = DiscreteCallback(condition, changesel!)
    prob = ODEProblem(bacphage!, u0, tspan, par)
    sol = solve(prob,RadauIIA5(), callback=cb)
    shortsol = sol(0.0:0.5:tend)
    # test = figure()
    # plot(sol.t, sol.u)
    conjdata = [conjugation(C, par.r) for C in shortsol[1, 1:end]]
    selecdata = [selection(C, par.s) for C in shortsol[1,1:end]]
    lysodata = [lysogeny(C, par.b) for C in shortsol[1,1:end]]

    test = figure()
    subplot(2,1,1)
    plot(shortsol.t, shortsol.u)
    subplot(2,1,2)
    plot(shortsol.t, conjdata, color="blue", label="Conjugation")
    plot(shortsol.t, selecdata, color="orange", label="Selection")
    plot(shortsol.t, lysodata, color="red", label="Lysogeny")
    xlabel("Time")
    ylabel("d(comp)/dt")
    legend()
    tight_layout()
    return test
end

#cumulative HGT graph as s changes from one value to larger value and from larger value to smaller value


#NOTE that selection is not HGT!!!*** but r and b are

#graphing reducing b but keeping s width the same
let 
    u0=[0.5]
    tsend=10000.0
    tspan=(0.0, tsend)
    freq = .1
    par = BacPhageSineForcedPar(b = 0.00001, per=0.5, amp=0.1, mid=-0.002)
    prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
    sol = solve(prob, RadauIIA5())
    solseries = sol(tsend-50.0:freq:tsend)
    selection = [sel_sine(par, t) for t in tsend-50.0:freq:tsend]
    seriesattractor = attractordata(selection, par)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, seriesattractor)
    return test
end

#what is the speed of "recovery" with and without lysogeny
#so problem - depending on selection if too negative then should lose gene BUT because continuous model never reach exactly 0.00000
#but we could use that to our advantage?

let 
    u0 = [0.5]
    tspan=(0.0, 900.0)

    condition(u,t,integrator) = t > 200.0 
    function changesel!(integrator)
        integrator.p.s=0.01
    end

    cb = DiscreteCallback(condition, changesel!)
    prob = ODEProblem(bacphage!, u0, tspan, BacPhagePar(s=-0.25, b=0.0)) ###NEED TO CHECK IF WORKS WHEN PUT b=0.0 instead of wobac function
    sol = solve(prob, Rodas5(), callback=cb)
    shortsol = sol(0.0:0.5:900.0)
    # test = figure()
    # plot(sol.t, sol.u)
    # return test
    conjdata = [conjugation(C, 0.1) for C in shortsol[1, 1:end]]
    selecdata = [selection(C, 0.01) for C in shortsol[1,1:end]]
    # lysodata = [lysogeny(C, 0.01) for C in shortsol[1,1:end]]

    test = figure()
    plot(shortsol.t, conjdata, color="blue", label="Conjugation")
    plot(shortsol.t, selecdata, color="orange", label="Selection")
    xlabel("Time")
    ylabel("d(comp)/dt")
    legend()
    # plot(shortsol.t, lysodata, color="red")
    return test
end

#all change in proportion of carrier is dominated by conjugation
#but with lysogeny shift in proportion much faster than without lysogeny

#REVISIT PROOF cause r and b dictate where equilbria are
#proof - #not sure I have taken into account initial starting conditions of conjugation and bacteriophage

@vars x

g(x) = -2*(x^3) + 4 * (x^2) -2*x

SymPy.solve(g(x), x)

eq2= 0.1:0.1:1.0

test2 = rand(Float64, (2, 2))
test2[1,2]
test2[2,1]

function conj_proof(eq1,eq2)
    return (2 * eq2 - eq2^2 - 2 * eq1 + eq1^2) / (3 * eq2^2 - 2 * eq2^3 - 3 * eq1^2 + 2 * eq1^3)
end

conj_proof(0.5,1.0)
conj_proof(0.6,0.5)

#seems to be symmetrical! what does this mean?

#rows are eq2
#columns are eq1
let
    eq1 = eq2 = 0.1:0.1:1.0
    test = zeros(10,10)
    for i in 1:10
        for j in 1:10
            test[i,j] = (2 * eq2[i] - eq2[i]^2 - 2 * eq1[j] + eq1[j]^2) / (3 * eq2[i]^2 - 2 * eq2[i]^3 - 3 * eq1[j]^2 + 2 * eq1[j]^3)
        end
    end
    return test
end


function dcdt(C,par)
    @unpack r, b, s = par
    return r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + b * ( 1 - C )
end

dcdt(0.4499999, BacPhagePar(s = 0.01))

function dHGTdt(C,dcdtval, par)
    @unpack r, b = par
    return 100 * (r * dcdtval - 2 * r * C * dcdtval - b * dcdtval)
end

dHGTdt(0.449999, dcdt(0.4499999, BacPhagePar(s = 0.01)), BacPhagePar(s = 0.01))

#time taken to reach gene fixation after switch in selection for different values of b

function calc_time_fixation(b, oldsel, newsel, tvals)
    solseries = selection_switch_bacphage(b, oldsel, newsel, tvals)
    time = []
    for ti in eachindex(solseries.t)
        if isapprox(solseries.u[ti][1], 1.0, atol=0.1)
            time = solseries.t[ti]-50.0
            break
        end
    end
    return time
end

let 
    data = selection_switch_bacphage(0.0001, -0.05, 0.01, 0.0:1.0:10000.0)
    test = figure()
    plot(data.t, data.u)
    return test
end

function time_fixation_b(brange, oldsel, newsel, tvals)
    timedata=zeros(length(brange))
    @threads for bi in eachindex(brange)
        timedata[bi] = calc_time_fixation(brange[bi], oldsel, newsel, tvals)
    end
    return [brange, timedata]
end

let 
    data = time_fixation_b(0.0001:0.0001:0.01, -0.05, 0.01, 0.0:1.0:10000.0)
    delayfigure = figure()
    plot(data[1], data[2])
    xlabel("b")
    ylabel("Time")
    # return delayfigure
    savefig(joinpath(abpath(), "figs/delay_selectionswitch.png"))
end

## PROOF 2 BUT INCREASING BACTERIOPHAGE REDUCES THE DELAY BETWEEN EQUILIBRIUM AND SYSTEM STATE. ALSO REDUCED CV
#I should be able to show that b/r never affects delay due to the eigenvalue function just shifting but not flattening.  #ACTUALLY b does affect delay - just very small
#So yes 1 fixed point - dynamics move off 1 earlier when b is smaller but overall trough in dynamics I think is the smaller (but should check this with higher resolution)
#look at the speed of down and up (with same sine wave) 
#draw equilibria on graph


# let #Rodas5 - CAN LIKELY REMOVE
#     par1 = BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.5, mid=0.0)
#     par2 = BacPhageSineForcedPar(b = 0.004, per=0.2, amp=0.5, mid=0.0)
#     u0 = [0.5]
#     tspan=(0.0, 500.0)
#     prob1 = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par1)
#     sol1 = solve(prob1, Rodas5())
#     solseries1 = sol1(50.0:0.05:100.0)
#     prob2 = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par2)
#     sol2 = solve(prob2, Rodas5())
#     solseries2 = sol2(50.0:0.05:100.0)
#     sel = [sel_sine(par1, t) for t in solseries1.t]
#     equil1 = [stableequil_ver1(s, par1) for s in sel]
#     equil2 = [stableequil_ver1(s, par2) for s in sel]
#     test = figure()
#     plot(solseries1.t, solseries1.u, color="green", label = "b=0.001")
#     plot(solseries2.t, solseries2.u, color="red", label = "b=0.004")
#     plot(solseries1.t, equil1, color="green", linestyle="dashed")
#     plot(solseries2.t, equil2, color="red", linestyle="dashed")
#     ylim(-0.1,1.1)
#     legend()
#     return test
#     # savefig(joinpath(abpath(), "figs/delay_bvalues.png"))
# end

#eigenvalues
let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data1 = [interior_equil_ver2(s, BacPhagePar()) for s in srange]
    data2 = [interior_equil_ver2(s, BacPhagePar(b=0.02)) for s in srange]
    data3 = [interior_equil_ver2(s, BacPhagePar(b=0.09)) for s in srange]
    bifurcplot = figure(figsize=(5,4))
    plot(srange, data1, color = "blue", label = "b = 0.01")
    plot(srange, data2, color = "green", label = "b = 0.02")
    plot(srange, data3, color = "red", label = "b = 0.09")
    ylabel("Ĉ")
    xlabel("s")
    xlim(-1.00, 1.0)
    ylim(0, 1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    legend()
    title("Linear Selection Function")
    # return bifurcplot
    savefig(joinpath(abpath(), "figs/bifurcver2_changingb.png"))
end

#Proof of 2 is
#Delay: geometric examination of equilibrium as s is changed. increasing b shifts the equilibrium "in" and the intersection of the system state with the equilbrium is the minimum which is always "in"
#CV: larger interior equilibrium for increased b. plus smaller push from 1 and larger pull from 1 for increased b

#Other ideas
#We know that CV will be reduced (minimum C) because b shifts the interior equilibrium up
# but increasing b decreases the pushing force away from 1 (lower eigenvalue)
#Maybe proof of decreasing delay is geometric - the interior equilibrium line has to cross with the system state line
# how to proove that intersection of equilibrium and system state is the minimum of the system state - not sure need to proove - just logical that when equilibria above system state then has to increase

## PROOF 3 DEPENDING ON THE PARAMETERS OF r, b, per, and amp, we get different qualitative behaviours

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

#eigen integral bifurcation #C=1 #mid
function bifurcintegral_eigen1_mid(midrange, tsend)
    data = zeros(length(midrange), 3)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for midi in eachindex(midrange)
        par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=midrange[midi])
        data[midi, 1] = calc_integral_eigen1(par)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000.0:1.0:tsend)
        data[midi, 2] = maximum(solseries)
        data[midi, 3] = minimum(solseries)
    end
    return data
end


#"bifurcations" changing b (keeping in case need to show b and r do the same thing as changing mid)
function bifurcb(brange)
    par = BacPhageSineForcedPar(r = 0.001, per=0.5, amp=0.4, mid=-0.01)
    maxC = zeros(length(brange))
    minC = zeros(length(brange))
    u0=[0.5]
    tspan=(0.0, 10000.0)
    for bi in eachindex(brange)
        par.b = brange[bi]
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(9000.0:1.0:10000.0)
        maxC[bi] = maximum(solseries)
        minC[bi] = minimum(solseries)
    end
    return [maxC, minC]
end

#"bifurcations" changing r (keeping in case need to show b and r do the same thing as changing mid)
function bifurcr(rrange)
    par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=-0.01)
    maxC = zeros(length(rrange))
    minC = zeros(length(rrange))
    u0=[0.5]
    tspan=(0.0, 10000.0)
    for ri in eachindex(rrange)
        par.r = rrange[ri]
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(9000.0:1.0:10000.0)
        maxC[ri] = maximum(solseries)
        minC[ri] = minimum(solseries)
    end
    return [maxC, minC]
end


function attractordata(selectiondata, par)
    attractor = zeros(length(selectiondata))
    for i in 1:length(selectiondata)
        attractor[i] = stableequil(selectiondata[i], par)
    end
    return attractor
end  #USED in tracking analysis

#using CV to think about tracking attractor
function trackattractor_CV(solutiondata, attractordata)
    solutionCV = std(solutiondata)/mean(solutiondata)
    attractorCV = std(attractordata)/mean(solutiondata)
    return [solutionCV, attractorCV]
end

function trackingcor_CV_b(brange, rval, freq, tsend)
    data = zeros(length(brange), 4)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for bi in eachindex(brange)
        par = BacPhageSineForcedPar(b = brange[bi], r=rval, per=0.5, amp=0.4, mid=-0.002)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000.0:freq:tsend)
        selection = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
        seriesattractor = attractordata(selection, par)
        CV = trackattractor_CV(solseries, seriesattractor)
        data[bi, 1] = brange[bi]
        data[bi, 2] = CV[1]
        data[bi, 3] = CV[2]
        data[bi, 4] = CV[1]/CV[2]
    end
    return data
end

let 
    data_lowr = trackingcor_CV_b(0.00001:0.00001:0.003, 0.0, 0.1, 10000.0)
    data_medr = trackingcor_CV_b(0.00001:0.00001:0.003, 0.0005, 0.1, 10000.0)
    data_highr = trackingcor_CV_b(0.00001:0.00001:0.003, 0.001, 0.1, 10000.0)
    stability_b = figure()
    subplot(1,3,1)
    plot(data_lowr[:,1], data_lowr[:,4], color="blue")
    ylabel("CV(Solution)/CV(Attractor)")
    xlabel("bacteriophage (b)")
    subplot(1,3,2)
    plot(data_medr[:,1], data_medr[:,4], color="blue")
    ylabel("CV(Solution)/CV(Attractor)")
    xlabel("bacteriophage (b)")
    subplot(1,3,3)
    plot(data_highr[:,1], data_highr[:,4], color="blue")
    ylabel("CV(Solution)/CV(Attractor)")
    xlabel("bacteriophage (b)")
    tight_layout()
    return stability_b
    # savefig(joinpath(abpath(), "figs/conjugation_stability.pdf"))
end

function trackingcor_CV_r(rrange, bval, freq, tsend)
    data = zeros(length(rrange), 4)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for ri in eachindex(rrange)
        par = BacPhageSineForcedPar(r = rrange[ri], b = bval, per=0.5, amp=0.4, mid=-0.002)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000.0:freq:tsend)
        selection = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
        seriesattractor = attractordata(selection, par)
        CV = trackattractor_CV(solseries, seriesattractor)
        data[ri, 1] = rrange[ri]
        data[ri, 2] = CV[1]
        data[ri, 3] = CV[2]
        data[ri, 4] = CV[1]/CV[2]
    end
    return data
end

let 
    data_lowb = trackingcor_CV_r(0.00001:0.00001:0.003, 0.0, 0.1, 10000.0)
    data_medb = trackingcor_CV_r(0.00001:0.00001:0.003, 0.0005, 0.1, 10000.0)
    data_highb = trackingcor_CV_r(0.00001:0.00001:0.003, 0.001, 0.1, 10000.0)
    stability_r = figure()
    subplot(1,3,1)
    plot(data_lowb[:,1], data_lowb[:,4], color="blue")
    ylabel("CV(Solution)/CV(Attractor)")
    xlabel("conjugation (r)")
    subplot(1,3,2)
    plot(data_medb[:,1], data_medb[:,4], color="blue")
    ylabel("CV(Solution)/CV(Attractor)")
    xlabel("conjugation (r)")
    subplot(1,3,3)
    plot(data_highb[:,1], data_highb[:,4], color="blue")
    ylabel("CV(Solution)/CV(Attractor)")
    xlabel("conjugation (r)")
    tight_layout()
    return stability_r
    # savefig(joinpath(abpath(), "figs/conjugation_stability.pdf"))
end

let
    u0=[0.5]
    tsend = 10000.0
    freq=0.1
    tspan=(0.0, tsend)
    par = BacPhageSineForcedPar(r = 0.00001, b = 0.001, per=0.5, amp=0.4, mid=-0.002)
    prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
    sol = solve(prob, RadauIIA5())
    solseries = sol(tsend-1000.0:freq:tsend)
    selection = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
    seriesattractor = attractordata(selection, par)
    test = figure()
    plot(solseries.t,solseries.u)
    plot(solseries.t, seriesattractor)
    return test
end

let
    u0=[0.5]
    tsend = 10000.0
    freq=0.1
    tspan=(0.0, tsend)
    par = BacPhageSineForcedPar(r = 0.001, b = 0.00001, per=0.5, amp=0.4, mid=-0.002)
    prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
    sol = solve(prob, RadauIIA5())
    solseries = sol(tsend-1000.0:freq:tsend)
    selection = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
    seriesattractor = attractordata(selection, par)
    test = figure()
    plot(solseries.t,solseries.u)
    plot(solseries.t, seriesattractor)
    return test
end

bifurc1 = [interior_equil(s, BacPhageSineForcedPar(r = 0.001, b = 0.00001, per=0.5, amp=0.4, mid=-0.002)) for s in -0.01:0.001:0.01]

let 
    srange = -0.01:0.00001:-0.00101
    bifurc1 = [interior_equil(s, BacPhageSineForcedPar(r = 0.001, b = 0.00001, per=0.5, amp=0.4, mid=-0.002)) for s in -0.01:0.00000001:-0.00101]
    bifurc2 = [interior_equil(s, BacPhageSineForcedPar(r = 0.00001, b = 0.001, per=0.5, amp=0.4, mid=-0.002)) for s in -0.01:0.0001:-0.00101]
    test = figure()
    plot(-0.01:0.00000001:-0.00101, bifurc1, color="blue")
    plot(-0.01:0.0001:-0.00101, bifurc2, color="red")
    ylim(0,1)
    xlim(-0.01,0.01)
    return test
end

let 
    srange = -0.01:0.00001:-0.00101
    bifurc1 = [interior_equil(s, BacPhageSineForcedPar(r = 0.001, b = 0.001, per=0.5, amp=0.4, mid=-0.002)) for s in -0.01:0.00000001:-0.00101]
    bifurc2 = [interior_equil(s, BacPhageSineForcedPar(r = 0.00001, b = 0.001, per=0.5, amp=0.4, mid=-0.002)) for s in -0.01:0.0001:-0.00101]
    test = figure()
    plot(-0.01:0.00000001:-0.00101, bifurc1, color="blue")
    plot(-0.01:0.0001:-0.00101, bifurc2, color="red")
    ylim(0,1)
    xlim(-0.01,0.01)
    return test
end


function trackingcor_CV_br(brange, rrange, freq, tsend)
    data = zeros(length(brange)+1, length(rrange)+1)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for bi in eachindex(brange)
        for ri in eachindex(rrange)
        par = BacPhageSineForcedPar(b = brange[bi], r=rrange[ri], per=0.5, amp=0.4, mid=-0.002)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000.0:freq:tsend)
        selection = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
        seriesattractor = attractordata(selection, par)
        CV = trackattractor_CV(solseries, seriesattractor)
        data[bi+1,1] = brange[bi]
        data[1,ri+1] = rrange[ri]
        data[bi+1, ri+1] = CV[1]/CV[2]
        end
    end
    return data
end

test = trackingcor_CV_br(0.00:0.001:0.003, 0.00001:0.001:0.003, 0.1, 10000)


function splitting_trackingdata(trackingdata, r_rb)
    if r_rb == "r"
        return hcat(trackingdata[1,2:end],trackingdata[2,2:end])
    elseif r_rb == "rb"
        rbdata = zeros((size(trackingdata,1)-2)*(size(trackingdata,2)-1), 2)
        for i in 3:size(trackingdata,1)
            for j in 2:size(trackingdata, 2)
                rbdata[i+j-4,1] = trackingdata[i,1]+trackingdata[1,j]
                rbdata[i+j-4,2] = trackingdata[i,j]
            end
        end
        return rbdata
    else
        error("r_rb should either be r or rb")
    end
end
#need to fix overwriting data problem
splitting_trackingdata(test, "r")
splitting_trackingdata(test, "rb")

#using range to think about tracking attractor
function trackattractor_range(solutiondata, attractordata)
    solutionrange = maximum(solutiondata) - minimum(solutiondata)
    attractorrange = maximum(attractordata) - minimum(attractordata)
    return [solutionrange, attractorrange]
end

function trackingcor_range_b(brange, freq, tsend)
    data = zeros(length(brange), 4)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for bi in eachindex(brange)
        par = BacPhageSineForcedPar(b = brange[bi], per=0.5, amp=0.4, mid=-0.002)
        prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
        sol = solve(prob, RadauIIA5())
        solseries = sol(tsend-1000.0:freq:tsend)
        selection = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
        seriesattractor = attractordata(selection, par)
        range = trackattractor_range(solseries, seriesattractor)
        data[bi, 1] = brange[bi]
        data[bi, 2] = range[1]
        data[bi, 3] = range[2]
        data[bi, 4] = range[1]/range[2]
    end
    return data
end

let 
    data = trackingcor_range_b(0.0001:0.0001:0.003, 0.1, 10000.0)
    test = figure()
    plot(data[:,1], data[:,4], color="blue")
    # plot(data[:,1], data[:,3], color="red")
    return test
end