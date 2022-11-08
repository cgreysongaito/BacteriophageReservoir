function abpath()
    replace(@__DIR__, "scripts" => "")
end

@with_kw mutable struct BacPhagePar
    r = 0.1
    s = 0.1
    b = 0.01
end

@with_kw mutable struct BacPhageSinePar
    r = 0.1
    b = 0.01
    per = 0.5
    amp = 1.0
    # selec::Function = sine
end

function bacphage!(du, u, p, t,)
    @unpack r, b, s = p
    du[1] = r * u[1] * (1 - u[1]) + ( s / ( 1 + s * u[1] ) ) * u[1] * ( 1 - u[1] ) + b * ( 1 - u[1] )
    return
end

function bacphage_wobac!(du, u, p, t,)
    @unpack r, b, s = p
    du[1] = r * u[1] * (1 - u[1]) + ( s / ( 1 + s * u[1] ) ) * u[1] * ( 1 - u[1] )
    return
end #can probably remove this and just change the b parameter to zero

function bacphage_sine!(du, u, p, t,)
    @unpack r, b, per, amp = p
    C , S  = u
    du[1] = r * C * (1 - C) + ( S / ( 1 + S * C ) ) * C * ( 1 - C ) + b * ( 1 - C )
    du[2] = amp * per * cos(per * t)
    return
end

function bacphage_sine_sol(b, per, amp, tsend, tvals)
    par = BacPhageSinePar()
    par.b = b
    par.per = per
    par.amp = amp
    u0 = [0.5, 0.0]
    tspan=(0.0, tsend)
    prob = ODEProblem(bacphage_sine!, u0, tspan, par)
    sol = solve(prob)
    return solend = sol(tvals)
end


function noise_creation(r, len)
    white = rand(Normal(0.0, 0.5), Int64(len))
    intnoise = [white[1]]
    for i in 2:Int64(len)
        intnoise = append!(intnoise, r * intnoise[i-1] + white[i] )
    end
    c = std(white)/std(intnoise)
    meanintnoise = mean(intnoise)
    scalednoise = zeros(Int64(len))
    for i in 1:Int64(len)
        scalednoise[i] = c * (intnoise[i] - meanintnoise)
    end
    return scalednoise
end #produces noise with a certain autocorrelation. variance of the noise is scaled using method in Wichmann et al. 2005

function bacphage_pert_sol(b, u0, freq, r, seed, tsend, tvals)
    Random.seed!(seed)
    par = BacPhagePar()
    par.b = b
    noise = noise_creation(r, tsend / freq)
    count = 1
    tspan = (0.0, tsend)

    function pert_cb2(integrator)
        count += 1
        integrator.p.s = noise[count] #https://discourse.julialang.org/t/change-parameters-in-a-model-in-specific-time-with-differentialequations/36930
    end

    cb = PeriodicCallback(pert_cb2, freq, initial_affect = false)
    prob = ODEProblem(bacphage!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    solend = sol(tvals)
    return [[i,j] for (i,j) in zip(solend[1,:], noise[100:500])]
end

function equil(s, b::Float64=0.01)
    r = BacPhagePar().r
    if s > (-(b+r))/(b+r+1)
        return 1
    else
        return (-(b*s + r + s) - sqrt(b^2*s^2 - 2*b*r*s + 2*b*s^2 + r^2 + 2*r*s + s^2))/(2*r*s)
    end
end


function vector_prod(sol)
    soldata=zeros(length(sol))
    seldata=zeros(length(sol))
    for i in 1:length(sol)
        soldata[i] = sol[i][1]
        seldata[i] = sol[i][2]
    end
    return [soldata, seldata]
end
