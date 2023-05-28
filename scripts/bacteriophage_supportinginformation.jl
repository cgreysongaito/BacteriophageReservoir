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


#times series decomposition of conjugation, bacteriophage, selection
let 
    u0=[0.5]
    tsend = 10000.0
    freq = 0.1
    tspan=(0.0, tsend)
    par = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.003)
    prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
    sol = solve(prob, RadauIIA5())
    solseries = sol(tsend-100.0:freq:tsend)
    conjworkdata = [conjugation(C, par.r) for C in solseries[1, :]]
    lysoworkdata = [lysogeny(C, par.b) for C in solseries[1,:]]
    selectiondata = [sel_sine(par, t) for t in tsend-100.0:freq:tsend]
    selecworkdata = [selection(C, s) for (C,s) in zip(solseries[1,:],selectiondata)]
    # seriesattractor = attractordata(selection, par)
    test = figure()
    plot(solseries.t, conjworkdata, color="blue", label="Conjugation")
    plot(solseries.t, lysoworkdata, color="red", label="Bacteriophage")
    # plot(solseries.t, selecworkdata, color="green", label="Selection")
    # # plot(solseries.t, seriesattractor)
    return test
end

#definitely selection does most of the movement work! but the balance for sure dictates the actual dynamics
#anyway to proove the order of "work" of selection, then bacteriophage, then conjugation

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



#Different periodicities and amplitudes
#test of bifurcmid with different periodicity (because model without bac is dependent on periodicity) - ANSWER - periodicity does not affect "bifurcation" point
function bifurcmid_per(midrange, perval, tsend)
    data = zeros(length(midrange), 3)
    u0=[0.5]
    tspan=(0.0, tsend)
    @threads for midi in eachindex(midrange)
        par = BacPhageSineForcedPar(b = 0.001, per=perval, amp=0.4, mid=midrange[midi])
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

let 
    midrange = -0.01:0.0001:0.001
    bifurcmid_data = bifurcmid_per(midrange, 0.2, 100000.0)
    test = figure()
    plot(bifurcmid_data[:, 1], bifurcmid_data[:, 2], color="black")
    plot(bifurcmid_data[:, 1], bifurcmid_data[:, 3], color="black")
    return test
end

#Eigenvalues integral balancing (gut responses to environmental variation)
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