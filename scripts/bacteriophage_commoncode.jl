#Helper functions
function abpath()
    replace(@__DIR__, "scripts" => "")
end

#Basic model
@with_kw mutable struct BacPhagePar
    r = 0.001 #to match niehus 10^-4 top value
    s = 0.1
    b = 0.001
end

function bacphage!(du, u, p, t,)
    @unpack r, b, s = p
    du[1] = r * u[1] * (1 - u[1]) +  s * u[1] * ( 1 - u[1] ) + b * ( 1 - u[1] )
    return
end

function selection_switch_bacphage(b, oldsel, newsel, tvals)
    par = BacPhagePar(b = b, s = oldsel)
    u0 = [stableequil(oldsel, par)]
    tspan=(0.0, 900.0)

    condition(u,t,integrator) = t > 50.0 
    function changesel!(integrator)
        integrator.p.s=newsel
    end

    cb = DiscreteCallback(condition, changesel!)
    prob = ODEProblem(bacphage!, u0, tspan, par)
    sol = solve(prob,RadauIIA5(), callback=cb)
    solseries = sol(tvals)
    return solseries
end

function conjugation(C, r)
    return r * C * (1 - C)
end

function selection(C, s)
    return s * C * ( 1 - C )
end

function lysogeny(C, b)
    return  b * ( 1 - C )
end

function parse_r_b_s(solseries, par)
    @unpack r, s, b = par
    conjdata = [conjugation(C, r) for C in solseries[1, :]]
    selecdata = [selection(C, s) for C in solseries[1,:]]
    lysodata = [lysogeny(C, b) for C in solseries[1,:]]
    return [conjdata, selecdata, lysodata]
end

#Equilibrium and bifurcation functions
function stableequil(s, par)
    @unpack b, r = par
    if s > -b - r
        return 1
    else
        return -b / (r+s)
    end
end

function interior_equil(s, par)
    @unpack b, r = par
    return -b / (r+s)
end


function bifurc(par)
    @unpack b, r = par
    return -b - r
end

#Eigenvalues for equilibrium 1
function eigen1(s, par)
    @unpack b, r = par
    return - b - r - s
end

#Eigenvalues for interior equilibrium
function eigenint(s, par)
    @unpack b, r = par
    C = interior_equil(s, par)
    return -2⋅C⋅r - 2⋅C⋅s - b + r + s
end

function calc_integral_eigen1(par)
    sel = [sel_sine(par, t) for t in 0.0:0.0001:1000.0]
    maxminsel = [maximum(sel), minimum(sel)]
    integral, err = quadgk(s -> eigen1(s, par), maxminsel[2], maxminsel[1])
    return integral
end

function calc_integral_eigenint(par)
    sel = [sel_sine(par, t) for t in 0.0:0.0001:1000.0]
    maxminsel = [maximum(sel), minimum(sel)]
    integral, err = quadgk(s -> eigenint(s, par), maxminsel[2], maxminsel[1])
    return integral
end

#With sine environmental selection
@with_kw mutable struct BacPhageSineForcedPar
    r = 0.001 #changed to neihus parameter
    b = 0.001
    per = 0.5
    amp = 1.0
    mid = 0.0
    selec::Function = sel_sine
end

function sel_sine(p, t)
    @unpack amp, per, mid = p
    return amp * sin(per * t) + mid
end

function bacphage_sine_forced!(du, u, p, t,)
    @unpack r, b, per, amp = p
    s = p.selec(p,t)
    du[1] = r * u[1] * (1 - u[1]) + (s * u[1] * ( 1 - u[1] )) + b * ( 1 - u[1] )
    return
end


#Noise environmental selection

function noise_creation(μ, σ, corr, len, seed)
    Random.seed!(seed)
    white = rand(Normal(0.0, σ), Int64(len))
    intnoise = [white[1]]
    for i in 2:Int64(len)
        intnoise = append!(intnoise, corr * intnoise[i-1] + white[i] )
    end
    c = std(white)/std(intnoise)
    meanintnoise = mean(intnoise)
    scalednoise = zeros(Int64(len))
    for i in 1:Int64(len)
        scalednoise[i] = c * (intnoise[i] - meanintnoise)
    end
    recentrednoise = zeros(Int64(len))
    for i in 1:Int64(len)
        recentrednoise[i] = scalednoise[i]+μ
    end
    return recentrednoise
end #produces noise with a certain autocorrelation. variance of the noise is scaled using method in Wichmann et al. 2005

function bacphage_pert_sol(b, u0, freq, μ, σ, corr, seed, tsend, tvals)
    par = BacPhagePar()
    par.b = b
    par.s = μ
    noise = noise_creation(μ, σ, corr, tsend / freq, seed)
    count = 1
    tspan = (0.0, tsend)

    function pert_cb2(integrator)
        count += 1
        integrator.p.s = noise[count] #https://discourse.julialang.org/t/change-parameters-in-a-model-in-specific-time-with-differentialequations/36930
    end

    cb = PeriodicCallback(pert_cb2, freq, initial_affect = false)
    prob = ODEProblem(bacphage!, [u0], tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, alg=RadauIIA5()) #reltol = 1e-8
    solend = sol(tvals)
    return [solend[1,:], append!([μ], noise[Int64((minimum(tvals) / freq) + 1.0):Int64(tsend / freq)])]
end

