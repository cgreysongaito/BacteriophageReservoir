include("packages.jl")
include("bacteriophage_commoncode.jl")

#Proof 1 YES CONJUGATION DOES MOST OF THE WORK
##exploring contributions of conjugation versus selection versus lysogeny after a change to positive selection
#May be able to get same intuition from the equations - but time series still useful
function conjugation(C, r)
    return r * C * (1 - C)
end

function selection(C, s)
    return ( s / ( 1 + s * C ) ) * C * ( 1 - C )
end

function lysogeny(C, b)
    return  b * ( 1 - C )
end

let 
    u0 = [0.5]
    tspan=(0.0, 900.0)

    condition(u,t,integrator) = t > 200.0 
    function changesel!(integrator)
        integrator.p.s=0.01
    end

    cb = DiscreteCallback(condition, changesel!)
    prob = ODEProblem(bacphage!, u0, tspan, BacPhagePar(b = 0.001, s=-0.25))
    sol = solve(prob,Rodas5(), callback=cb)
    shortsol = sol(0.0:0.5:900.0)
    # test = figure()
    # plot(sol.t, sol.u)
    conjdata = [conjugation(C, 0.1) for C in shortsol[1, 1:end]]
    selecdata = [selection(C, 0.01) for C in shortsol[1,1:end]]
    lysodata = [lysogeny(C, 0.01) for C in shortsol[1,1:end]]

    test = figure()
    plot(shortsol.t, conjdata, color="blue", label="Conjugation")
    plot(shortsol.t, selecdata, color="orange", label="Selection")
    plot(shortsol.t, lysodata, color="red", label="Lysogeny")
    xlabel("Time")
    ylabel("d(comp)/dt")
    legend()
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
    prob = ODEProblem(bacphage_wobac!, u0, tspan, BacPhagePar(s=-0.25))
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


## PROOF 2 BUT INCREASING BACTERIOPHAGE REDUCES THE DELAY BETWEEN EQUILIBRIUM AND SYSTEM STATE. ALSO REDUCED CV
#I should be able to show that b/r never affects delay due to the eigenvalue function just shifting but not flattening.  #ACTUALLY b does affect delay - just very small
#So yes 1 fixed point - dynamics move off 1 earlier when b is smaller but overall trough in dynamics I think is the smaller (but should check this with higher resolution)
#look at the speed of down and up (with same sine wave) 
#draw equilibria on graph
let #Rodas5
    par1 = BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.5, mid=0.0)
    par2 = BacPhageSineForcedPar(b = 0.004, per=0.2, amp=0.5, mid=0.0)
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob1 = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par1)
    sol1 = solve(prob1, Rodas5())
    solseries1 = sol1(50.0:0.05:100.0)
    prob2 = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par2)
    sol2 = solve(prob2, Rodas5())
    solseries2 = sol2(50.0:0.05:100.0)
    sel = [sel_sine(par1, t) for t in solseries1.t]
    equil1 = [stableequil_ver1(s, par1) for s in sel]
    equil2 = [stableequil_ver1(s, par2) for s in sel]
    test = figure()
    plot(solseries1.t, solseries1.u, color="green", label = "b=0.001")
    plot(solseries2.t, solseries2.u, color="red", label = "b=0.004")
    plot(solseries1.t, equil1, color="green", linestyle="dashed")
    plot(solseries2.t, equil2, color="red", linestyle="dashed")
    ylim(-0.1,1.1)
    legend()
    return test
    # savefig(joinpath(abpath(), "figs/delay_bvalues.png"))
end

#eigenvalues
#version 1


let 
    srange = -0.2:0.01:0.1
    data_ver1a = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.0)) for s in srange]
    data_ver1b = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.001)) for s in srange]
    data_ver1c = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.01)) for s in srange]    
    test = figure()
    plot(srange, data_ver1a, color = "blue", label = "b = 0.00")
    plot(srange, data_ver1b, color = "green", label = "b = 0.001")
    plot(srange, data_ver1c, color = "purple", label = "b = 0.01")
    hlines(0.0, -0.2, 0.1, linestyles = "dashed")
    xlabel("s")
    ylabel("λ")
    legend()
    title("Eigenvalue for Ĉ = 1")
    return test
    # savefig(joinpath(abpath(), "figs/selectionfunction_eigenvalues1b.png"))
end
#b just shifts where the eigenvalue line sits but does not flatten the graph


let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data1 = [interior_equil_ver1(s, BacPhagePar(b=0.001)) for s in srange]
    data2 = [interior_equil_ver1(s, BacPhagePar(b=0.005)) for s in srange]
    data3 = [interior_equil_ver1(s, BacPhagePar(b=0.01)) for s in srange]
    bifurcplot = figure(figsize=(5,4))
    plot(srange, data1, color = "blue", label = "b = 0.001")
    plot(srange, data2, color = "green", label = "b = 0.005")
    plot(srange, data3, color = "red", label = "b = 0.01")
    ylabel("Ĉ")
    xlabel("s")
    xlim(-1.00, 1.0)
    ylim(0, 1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    legend()
    title("Non-Linear Selection Function")
    return bifurcplot
    # savefig(joinpath(abpath(), "figs/bifurcver1_changingb.png"))
end

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
##next step figure out way to find bifurcation values of different qualitative behaviours for version 1 and version 2
let #Rodas5
    par1 = BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.3, mid=-0.2)
    u0 = [0.5]
    tspan=(0.0, 800.0)
    prob1 = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par1)
    sol1 = solve(prob1, Rodas5())
    solseries1 = sol1(0.0:0.05:800.0)

    test = figure()
    plot(solseries1.t, solseries1.u)
    return test
end

#looks to be oscilattory attractor always (unless pushed to 1) but amplitude of oscilattory attractor dependent on parameters
#so periodicity doesn't really affect qualitative pattern (except for the number of oscillations). the eventual end point is the same
#amplitude and mid point do affect the qualitative patter though (keeping b and r the same)


#Best QUESTION - what makes go to fixation versus oscilating attractor?
let #Rodas5
    par1 = BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.3, mid=0.1)
    u0 = [0.5]
    tspan=(0.0, 600.0)
    prob1 = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par1)
    sol1 = solve(prob1, Rodas5())
    solseries1 = sol1(500.0:0.05:600.0)
    sel = [sel_sine(par1, t) for t in solseries1.t]
    equil1 = [stableequil_ver1(s, par1) for s in sel]
    test = figure()
    plot(solseries1.t, solseries1.u)
    # plot(solseries1.t, equil1, color="green", linestyle="dashed")
    return test
end

let #Rodas5
    par1 = BacPhageSineForcedPar(b = 0.001, per=0.1, amp=0.3, mid=-0.09999999)
    u0 = [0.5]
    tspan=(0.0, 2000.0)
    prob1 = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par1)
    sol1 = solve(prob1, Rodas5())
    solseries1 = sol1(0.0:0.05:2000.0)
    sel = [sel_sine(par1, t) for t in solseries1.t]
    equil1 = [stableequil_ver1(s, par1) for s in sel]
    test = figure()
    plot(solseries1.t, solseries1.u)
    plot(solseries1.t, equil1, color="green", linestyle="dashed")
    return test
end
#one difference is that the time spent at 1 is longer for higher mid. (shorter time inbetween 1)
#wonder whether we can get at speed travelling versus space must travel versus time available?

#do the same analysis with version 2 (below) then examine geometric mean fitness
bifurc_ver2(BacPhageSineForcedPar(b = 0.001, per=0.1, amp=0.3, mid=-0.09))


let 
    sel1 = [sel_sine(BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.3, mid=-0.01), t) for t in 0.0:0.1:100.0]
    maxminsel1 = [maximum(sel1), minimum(sel1)]
    sel2 = [sel_sine(BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.3, mid=-0.06), t) for t in 0.0:0.1:100.0]
    maxminsel2 = [maximum(sel2), minimum(sel2)]
    srange = -0.5:0.01:0.5
    data_ver1a = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.001)) for s in srange]   
    test = figure()
    plot(srange, data_ver1a, color = "blue", label = "b = 0.001")
    vlines(maxminsel1[1], -1, 1, color = "green", linestyles = "dashed")
    vlines(maxminsel1[2], -1, 1, color = "green", linestyles = "dashed")
    vlines(maxminsel2[1], -1, 1, color = "red", linestyles = "dashed")
    vlines(maxminsel2[2], -1, 1, color = "red", linestyles = "dashed")
    hlines(0.0, -0.5, 0.5, linestyles = "dashed")
    xlabel("s")
    ylabel("λ")
    legend()
    title("Eigenvalue for Ĉ = 1")
    return test
    # savefig(joinpath(abpath(), "figs/selectionfunction_eigenvalues1b.png"))
end

#maybe when integral of eigenvalues of balance we get balanced oscillation at 1. heavy on 1 side we get fixation. heavy on other side we get oscillator below 1
#doesn't seem to be balance of integral of eigenvalues because we can fixation when negative integral

calc_integral_eigen(BacPhageSineForcedPar(per = 0.2, amp=0.3, mid=-0.065))

bifurc_ver1(BacPhageSineForcedPar(per = 0.2, amp=0.3, mid=-0.065))

function balance_eigen_mid(midrange, par)
    parnew = par
    integral = zeros(length(midrange))
    for i in eachindex(midrange)
        parnew.mid = midrange[i]
        sel = [sel_sine(parnew, t) for t in 0.0:0.1:100.0]
        maxminsel = [maximum(sel), minimum(sel)]
        integral[i], err = quadgk(s -> eigen1_ver1(s, parnew), maxminsel[2], maxminsel[1]) 
    end
    for j in eachindex(integral)
        if isapprox(0.0, integral[j], atol=0.0001)
            return [midrange[j],integral[j]]
        end
    end
end

balance_eigen_mid(-0.2:0.00001:0.2, BacPhageSineForcedPar(per = 0.2, amp=0.3))


quadgk(s -> eigen1_ver1(s, BacPhageSineForcedPar()), -0.3, 0.2)


#i wonder if I could calculate the potential function for this model and see how s changes this potential function

function integral_calc_bifurc(data, bifurcvalue)
    integral_sum = 0.0
    for i in 1:length(data)
        if data[i] > bifurcvalue
            integral_sum += data[i]-bifurcvalue
        end
    end
    return integral_sum
end
#graph of integral above bifurc value for per and amp
function integral_surface_per(perrange, amprange, pardefault)
    bifurcvalue = bifurc_ver1(pardefault)
    trange = 0.0:0.1:100.0
    surface = zeros(length(perrange),length(amprange))
    for i in eachindex(perrange)
        for j in eachindex(amprange)
            parnew = pardefault
            parnew.per = perrange[i]
            parnew.amp = amprange[j]
            seldata = [sel_sine(parnew, t) for t in trange]
            surface[i,j] = integral_calc_bifurc(seldata, bifurcvalue)
        end
    end
    return surface
end

function integral_line_amp(amprange, pardefault)
    bifurcvalue = bifurc_ver1(pardefault)
    trange = 0.0:0.1:100.0
    intdata = zeros(length(amprange))
    for i in eachindex(amprange)
        parnew = pardefault
        parnew.amp = amprange[i]
        seldata = [sel_sine(parnew, t) for t in trange]
        intdata[i] = integral_calc_bifurc(seldata, bifurcvalue)
    end
    return intdata
end

let 
    n = 100
    perrange = range(0.01, stop=0.9, length=n)
    amprange = range(0.01, stop=1.0, length=n)
    xgrid = repeat(perrange',n,1)
    ygrid = repeat(amprange,1,n)
    data = integral_surface_per(perrange, amprange, BacPhageSineForcedPar())
    surface_figure = figure()
    plot_surface(perrange, amprange, data, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
    return surface_figure
end

#line
let 
    n = 100
    amprange = range(0.01, stop=1.0, length=n)
    data1 = integral_line_amp(amprange, BacPhageSineForcedPar(per=0.1))
    data2 = integral_line_amp(amprange, BacPhageSineForcedPar(per=0.4))
    data3 = integral_line_amp(amprange, BacPhageSineForcedPar(per=1.0))
    test = figure()
    plot(amprange, data1, label = "per = 0.1")
    plot(amprange, data2, label = "per = 0.4")
    plot(amprange, data3, label = "per = 1.0")
    legend()
    return test
end


function integral_surface_mid(midrange, amprange, par)
    bifurcvalue = bifurc_ver1(par)
    trange = 0.0:0.1:100.0
    surface = zeros(length(midrange),length(amprange))
    for i in 1:length(midrange)
        for j in 1:length(amprange)
            par.per = midrange[i]
            par.amp = amprange[j]
            seldata = [sel_sine(par, t) for t in trange]
            surface[i,j] = integral_calc_bifurc(seldata, bifurcvalue)
        end
    end
    return surface
end


function integral_line_mid(midrange, pardefault)
    bifurcvalue = bifurc_ver1(pardefault)
    trange = 0.0:0.1:100.0
    intdata = zeros(length(midrange))
    for i in eachindex(midrange)
        parnew = pardefault
        parnew.mid = midrange[i]
        seldata = [sel_sine(parnew, t) for t in trange]
        intdata[i] = integral_calc_bifurc(seldata, bifurcvalue)
    end
    return intdata
end

let 
    n = 100
    midrange = range(-0.5, stop=0.5, length=n)
    amprange = range(0.01, stop=1.0, length=n)
    xgrid = repeat(midrange',n,1)
    ygrid = repeat(amprange,1,n)
    data = integral_surface_per(midrange, amprange, BacPhageSineForcedPar())
    surface_figure = figure()
    plot_surface(midrange, amprange, data, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
    return surface_figure
end

#line
let 
    n = 100
    midrange = range(-0.5, stop=0.5, length=n)
    data1 = integral_line_mid(midrange, BacPhageSineForcedPar(amp=0.1))
    data2 = integral_line_mid(midrange, BacPhageSineForcedPar(amp=0.4))
    data3 = integral_line_mid(midrange, BacPhageSineForcedPar(amp=1.0))
    test = figure()
    plot(midrange, data1, label = "amp = 0.1")
    plot(midrange, data2, label = "amp = 0.4")
    plot(midrange, data3, label = "amp = 1.0")
    legend()
    return test
end

#graph of integral above bifurc value for per and amp

#trying geometric mean of selection (for bifurcation)
#fixation pattern
let #Rodas5
    par1 = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.3, mid=-0.000000001)
    u0 = [0.5]
    tspan=(0.0, 800.0)
    prob1 = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par1)
    sol1 = solve(prob1, Rodas5())
    solseries1 = sol1(0.0:0.05:800.0)

    test = figure()
    plot(solseries1.t, solseries1.u)
    return test
end

function geomean(par, a, b)
    @unpack per, amp, mid = par
    integral, err = quadgk(t -> log(Complex(sel_sine(par, t))), a, b)
    return exp((1 / (b-a)) * integral )
end

function sine_return(per)
    return 2 * pi / per
end

geomean(BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.3, mid=-0.000000001), 0.0, 2*sine_return(0.2))
geomean(BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.3, mid=0.0), 0.0, sine_return(0.2))
geomean(BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.3, mid=0.000000001), 0.0, 2*sine_return(0.2))
#oscillation pattern
let #Rodas5
    par1 = BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.3, mid=-0.102)
    u0 = [0.5]
    tspan=(0.0, 800.0)
    prob1 = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par1)
    sol1 = solve(prob1, Rodas5())
    solseries1 = sol1(0.0:0.05:800.0)

    test = figure()
    plot(solseries1.t, solseries1.u)
    return test
end

geomean(BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.3, mid=-0.102), 0.0, sine_return(0.2))
geomean(BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.3, mid=0.0), 0.0, sine_return(0.2))
geomean(BacPhageSineForcedPar(b = 0.001, per=0.2, amp=0.3, mid=0.102), 0.0, sine_return(0.2))


#eigenvalue changes for b and r
#C = 1
let 
    srange = -0.2:0.01:0.2
    data1 = [eigen1_ver2(s, BacPhageSineForcedPar(b = 0.001)) for s in srange]
    data2 = [eigen1_ver2(s, BacPhageSineForcedPar(b = 0.01)) for s in srange]
    data3 = [eigen1_ver2(s, BacPhageSineForcedPar(b = 0.1)) for s in srange]
    eigenb = figure()
    plot(srange, data1, color="blue", label = "b=0.001")
    plot(srange, data2, color="green", label = "b=0.01")
    plot(srange, data3, color="red", label = "b=0.1")
    xlabel("s")
    ylabel("λ")
    legend()
    return eigenb
end

let 
    srange = -0.2:0.01:0.2
    data1 = [eigen1_ver2(s, BacPhageSineForcedPar(r = 0.001)) for s in srange]
    data2 = [eigen1_ver2(s, BacPhageSineForcedPar(r = 0.01)) for s in srange]
    data3 = [eigen1_ver2(s, BacPhageSineForcedPar(r = 0.1)) for s in srange]
    eigenr = figure()
    plot(srange, data1, color="blue", label = "r=0.001")
    plot(srange, data2, color="green", label = "r=0.01")
    plot(srange, data3, color="red", label = "r=0.1")
    xlabel("s")
    ylabel("λ")
    legend()
    return eigenr
end

#Interior equilibrium
let 
    srange = -0.2:0.01:0.2
    data1 = [eigenint_ver2(s, BacPhageSineForcedPar(b = 0.001)) for s in srange]
    data2 = [eigenint_ver2(s, BacPhageSineForcedPar(b = 0.01)) for s in srange]
    data3 = [eigenint_ver2(s, BacPhageSineForcedPar(b = 0.1)) for s in srange]
    eigenintb = figure()
    plot(srange, data1, color="blue", label = "b=0.001")
    plot(srange, data2, color="green", label = "b=0.01")
    plot(srange, data3, color="red", label = "b=0.1")
    xlabel("s")
    ylabel("λ")
    legend()
    return eigenintb
end

let 
    srange = -0.2:0.001:0.2
    data1 = [eigenint_ver2(s, BacPhageSineForcedPar(r = 0.001)) for s in srange]
    data2 = [eigenint_ver2(s, BacPhageSineForcedPar(r = 0.01)) for s in srange]
    data3 = [eigenint_ver2(s, BacPhageSineForcedPar(r = 0.1)) for s in srange]
    eigenintr = figure()
    plot(srange, data1, color="blue", label = "r=0.001")
    plot(srange, data2, color="green", label = "r=0.01")
    plot(srange, data3, color="red", label = "r=0.1")
    xlabel("s")
    ylabel("λ")
    legend()
    return eigenintr
end #why is there a gap in the eigen lines  - when r=-S

#"bifurcations" changing b
function bifurcb(brange)
    par = BacPhageSineForcedPar(r = 0.001, per=0.5, amp=0.4, mid=-0.01)
    maxC = zeros(length(brange))
    minC = zeros(length(brange))
    u0=[0.5]
    tspan=(0.0, 10000.0)
    for bi in eachindex(brange)
        par.b = brange[bi]
        prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
        sol = solve(prob, Rodas5())
        solseries = sol(9000.0:1.0:10000.0)
        maxC[bi] = maximum(solseries)
        minC[bi] = minimum(solseries)
    end
    return [maxC, minC]
end

let 
    brange = 0.00001:0.00001:0.0088
    data = bifurcb(brange)
    test = figure()
    plot(brange, data[1])
    plot(brange, data[2])
    return test
end #need to deal with instabilities

#"bifurcations" changing r
function bifurcr(rrange)
    par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=-0.01)
    maxC = zeros(length(rrange))
    minC = zeros(length(rrange))
    u0=[0.5]
    tspan=(0.0, 10000.0)
    for ri in eachindex(rrange)
        par.r = rrange[ri]
        prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
        sol = solve(prob, Rodas5())
        solseries = sol(9000.0:1.0:10000.0)
        maxC[ri] = maximum(solseries)
        minC[ri] = minimum(solseries)
    end
    return [maxC, minC]
end

let 
    rrange = 0.00001:0.00001:0.0091
    data = bifurcr(rrange)
    test = figure()
    plot(rrange, data[1])
    plot(rrange, data[2])
    return test
end
#"bifurcations" changing mid
function bifurcmid(midrange)
    par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=-0.01)
    maxC = zeros(length(midrange))
    minC = zeros(length(midrange))
    u0=[0.5]
    tspan=(0.0, 10000.0)
    for midi in eachindex(midrange)
        par.mid = midrange[midi]
        prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
        sol = solve(prob, Rodas5())
        solseries = sol(9000.0:1.0:10000.0)
        maxC[midi] = maximum(solseries)
        minC[midi] = minimum(solseries)
    end
    return [maxC, minC]
end

let 
    midrange = -0.01:0.0001:0.001
    data = bifurcmid(midrange)
    test = figure()
    plot(midrange, data[1])
    plot(midrange, data[2])
    return test
end

#eigen integral bifurcation
#C=1
#mid
function bifurcintegral_eigen1_mid(midrange)
    data = zeros(length(midrange), 3)
    u0=[0.5]
    tspan=(0.0, 10000.0)
    @threads for midi in eachindex(midrange)
        par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=midrange[midi])
        data[midi, 1] = calc_integral_eigen1_ver2(par)
        prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
        sol = solve(prob, Rodas5())
        solseries = sol(9000.0:1.0:10000.0)
        data[midi, 2] = maximum(solseries)
        data[midi, 3] = minimum(solseries)
    end
    return data
end

let 
    midrange = -0.01:0.0001:0.001
    data = bifurcintegral_eigen1_mid(midrange)
    test = figure()
    plot(data[:, 1], data[:, 2])
    plot(data[:, 1], data[:, 3])
    vlines(0.0, 0.0, 1.0, linestyles="dashed", color="black")
    return test
end

#int C
#mid

function bifurcintegral_eigenint_mid(midrange)
    data = zeros(length(midrange), 3)
    u0=[0.5]
    tspan=(0.0, 10000.0)
    @threads for midi in eachindex(midrange)
        par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=midrange[midi])
        data[midi, 1] = calc_integral_eigenint_ver2(par)
        prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
        sol = solve(prob, Rodas5())
        solseries = sol(9000.0:1.0:10000.0)
        data[midi, 2] = maximum(solseries)
        data[midi, 3] = minimum(solseries)
    end
    return data
end

let 
    midrange = -0.01:0.0001:0.001
    data = bifurcintegral_eigenint_mid(midrange)
    test = figure()
    plot(data[:, 1], data[:, 2])
    plot(data[:, 1], data[:, 3])
    return test
end

let #eigenvalues are just "conjugates" of each other so adding or subtracting eigenvalue integrals does nothing
    srange = -0.2:0.01:0.2
    data1 = [eigen1_ver2(s, BacPhageSineForcedPar(b = 0.001)) for s in srange]
    dataint = [eigenint_ver2(s, BacPhageSineForcedPar(b = 0.001)) for s in srange]
    eigenb = figure()
    plot(srange, data1, color="blue")
    plot(srange, dataint, color="green")
    xlabel("s")
    ylabel("λ")
    legend()
    return eigenb
end


#relative time spent on 1 versus interior in a cycle - "bifurcation"
function time_equil1(par)
    bifurcval=bifurc_ver2(par)
    if (bifurcval-par.mid)/par.amp > 1.0 || (bifurcval-par.mid)/par.amp < -1
        return 0.0
    end
    t = asin((bifurcval-par.mid)/par.amp)/par.per
    if par.mid > bifurcval
        return (abs(t)+((pi/par.per)+abs(t)))/(2*pi/par.per)
    elseif par.mid < bifurcval
        return (((pi/par.per)-abs(t))-abs(t))/(2*pi/par.per)
    else
        return 0.5
    end
end

let 
    bifurc = -0.2
    xrange = 0.0:0.01:2*pi
    sinedata1 = [sin(x)+0.0 for x in xrange]
    sinedata2 = [sin(x)-0.1 for x in xrange]
    sinedata3 = [sin(x)-0.3 for x in xrange]
    test = figure()
    subplot(3,1,1)
    plot(xrange, sinedata1)
    hlines(bifurc, 0.0, 2*pi, linestyles = "dashed")
    hlines(0.0, 0.0, 2*pi)
    println(bifurc-0.0)
    subplot(3,1,2)
    plot(xrange, sinedata2)
    hlines(bifurc, 0.0, 2*pi, linestyles = "dashed")
    hlines(-0.1, 0.0, 2*pi)
    println(bifurc-(-0.1))
    subplot(3,1,3)
    plot(xrange, sinedata3)
    hlines(bifurc, 0.0, 2*pi, linestyles = "dashed")
    hlines(-0.3, 0.0, 2*pi)
    println(bifurc-(-0.3))
    return test
end


###NOT SURE IF THIS IS WORKING
let 
    par = BacPhageSineForcedPar(b = 0.001)
    trange = 0.0:0.01:2*pi/par.per
    data = [sel_sine(par, t) for t in trange]
    test = figure()
    plot(trange, data)
    hlines(bifurc_ver2(par), 0.0, 2*pi/par.per)
    return test
end


#"bifurcation" integral of sine wave above 0 and above bifurcation (shift sine wave by bifurcation value)

#gemometric mean - "bifucation"

#CHECK IF INITIAL VALUE PROBLEM
#PROOF 4 GENERALIZE FOR white and red noise


### Numerical Proof
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
    xlabel()
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
    xlabel("Periodicity")
    ylabel("CV")
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

let 
    sol = bacphage_sine_sol(0.01,0.2, 0.5, 500.0, 100.0:1.0:500.0)
    return std(sol[2,:])
end

let 
    sol = bacphage_sine_sol(0.01,1.0, 0.5, 500.0, 100.0:1.0:500.0)
    test=figure()
    plot(sol.t, sol.u)
    return test
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

let
    perrange = 0.1:0.1:1.0
    data = StandCV_calc_per(perrange)
    test = figure()
    plot(perrange, data[3], color = "blue")
    plot(perrange, data[1], color = "green")
    plot(perrange, data[2], color = "red")
    xlabel("Periodicity")
    ylabel("CV")
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

let
    data = tracking_sine_per(0.1:0.1:0.9) #longer periods obviously mean better ability of the system to track the equilibrium
    test = figure()
    plot(0.1:0.1:0.9, data)
    xlabel("Periodicity")
    ylabel("Mean difference between solution and equilibrium")
    return test
end

#temporal cross-correlation of equilibrium state with actual dynamics
#testing of crosscor

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

let 
    data = trackingcor_sine_per(0.1:0.1:0.9, 1:1:100)
    test = figure()
    plot(0.1:0.1:0.9, data)
    xlabel("Periodicity")
    ylabel("Delay (from crosscorrelation)")
    return test
end

function trackingcor_sine(bval)
    sol = bacphage_sine_sol(bval, 0.05, 1.0, 500.0, 100.0:0.01:500.0)
    return vector_prod(sol)
    # equilvals = calc_equil(splitsolsel[2])
    # cordata = crosscor(splitsolsel[1], equilvals, lrange)
    # return findmax(cordata)[2]
end


let 
    par = BacPhagePar(b=0.01, s=-0.2)
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage!, u0, tspan, par)
    return sol = solve(prob)
    # test = figure()
    # plot(sol.t, sol.u)
    # return test
end

let 
    testdata = trackingcor_sine(0.01)
    test = figure()
    plot(100.0:0.01:500.0, testdata[1])
    plot(100.0:0.01:500.0, testdata[2])
    return test
end

#Note - increasing b increases the lower limit of the equilibrium value

#tracking with changing b (still sine wave)
function trackingcor_sine_b(brange, lrange)
    trackcor = zeros(length(brange))
    for (bi, bval) in enumerate(brange)
        
    end
    return trackcor
end




let 
    data = trackingcor_sine_b(0.0:0.001:0.06, 1:1:100)
    test = figure()
    plot(0.0:0.001:0.06, data)
    xlabel("b")
    ylabel("Delay (from crosscorrelation)")
    return test
end


#DO I NEED TO STANDARDIZE BY UNDERLYING TIME SCALE OF SELECTION WHEN CHANGING PERIOD

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



#check eigenvalues of 1 fixed point as change selection for version 1 and version 2
#version 1
function eigen1_ver1(s, par)
    @unpack b, r = par
    return (-b*s - b - r*s - r - s)/(s+1)
end

function eigen1_ver2(s, par)
    @unpack b, r = par
    return - b - r - s
end

let 
    srange = -0.5:0.01:1.00
    data_ver1a = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.0)) for s in srange]
    data_ver2 = [eigen1_ver2(s, BacPhageSineForcedPar()) for s in srange]
    test = figure()
    plot(srange, data_ver1a, color = "blue")

    plot(srange, data_ver2, color = "red")
    hlines(0.0, -0.5, 1.00)
    xlabel("s")
    ylabel("λ")
    return test
end
#b just shifts where the eigenvalue line sits but does not flatten the graph


let 
    srange = -0.5:0.01:1.00
    data_ver1a = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.01)) for s in srange]
    data_ver2 = [eigen1_ver2(s, BacPhageSineForcedPar()) for s in srange]
    test = figure()
    plot(srange, data_ver1a, color = "blue", label="Nonlinear")
    plot(srange, data_ver2, color = "red", label="Linear")
    hlines(0.0, -0.5, 1.00, linestyles = "dashed")
    xlabel("s")
    ylabel("λ")
    legend()
    title("Eigenvalue for Ĉ = 1")
    # return test
    savefig(joinpath(abpath(), "figs/selectionfunction_eigenvalues1.png"))
end

let 
    srange = -0.5:0.01:1.00
    data_ver1a = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.0)) for s in srange]
    data_ver1b = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.04)) for s in srange]
    data_ver1c = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.08)) for s in srange]    
    test = figure()
    plot(srange, data_ver1a, color = "blue", label = "b = 0.00")
    plot(srange, data_ver1b, color = "green", label = "b = 0.04")
    plot(srange, data_ver1c, color = "purple", label = "b = 0.08")
    hlines(0.0, -0.5, 1.00, linestyles = "dashed")
    xlabel("s")
    ylabel("λ")
    legend()
    title("Eigenvalue for Ĉ = 1")
    # return test
    savefig(joinpath(abpath(), "figs/selectionfunction_eigenvalues1b.png"))
end
#r also doesn't flatten the eigenvalue line, just shifts the line


#why does version 2 go to fixation or 0 (can it ever find oscilating "attractor")?
#version 1 can go to fixation just different point where this happens - so 
#next step figure out way to find bifurcation values of different qualitative behaviours for version 1 and version 2
#maybe "bifrucation" is dependent of delay of system state matching equilibrium (i.e if we make the delay longer then quicker to fixation OR if we flatten the equilibria line then also quicker to fixation)
#I should be able to show that b/r never affects delay due to the eigenvalue function just shifting but not flattening.
let #Rodas5
    par = BacPhageSineForcedPar()
    par.b = 0.02
    par.per = 0.2
    par.amp= 0.5
    par.mid = 0.00
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par)
    sol = solve(prob, Rodas5())
    solseries = sol(50.0:0.05:100.0)
    sel = [sel_sine(par, t) for t in solseries.t]
    equil = [equil_ver1(s, par.b) for s in sel]
    equil2 = [equil_ver1(s, 0.01) for s in sel]
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, sel)
    plot(solseries.t, equil, color="green")
    plot(solseries.t, equil2, color="red")
    hlines(bifurc_ver1(par), 50.0, 100.0, colors="green" )
    ylim(-1,1.6)
    return test
end


let #Rodas5
    par = BacPhageSineForcedPar()
    par.b = 0.01
    par.per = 0.2
    par.amp= 0.5
    par.mid = 0.0
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
    sol = solve(prob, Rodas5())
    solseries = sol(50.0:1.0:500.0)
    return solseries
end


