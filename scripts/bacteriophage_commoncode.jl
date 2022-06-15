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
    # selec::Function = sine
end

function bacphage!(du, u, p, t,)
    @unpack r, b, s = p
    du = r * u * (1 - u) + ( s / ( 1 + s * u ) ) * u * ( 1 - u ) + b * ( 1 - u )
    return
end

function bacphage_sine!(du, u, p, t,)
    @unpack r, b, per = p
    C , S  = u
    du[1] = r * C * (1 - C) + ( S / ( 1 + S * C ) ) * C * ( 1 - C ) + b * ( 1 - C )
    du[2] = per * cos(per * t)
    return
end

function noise_creation(r, len)
    white = rand(Normal(0.0, 0.01), Int64(len))
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

function bacphage_pert(b, freq, r, seed, tsend, tvals)
    Random.seed!(seed)
    par = BacPhagePar()
    par.b = b
    noise = noise_creation(r, tsend / freq)
    count = 1
    u0 = [eq_II(par)[1], eq_II(par)[2] + noise[1]] #need to fix
    tspan = (0.0, tsend)

    function pert_cb2(integrator)
        count += 1
        if isapprox(integrator.u[2], 0.00000000; atol = 1e-8)
            integrator.u[2] = 0.00000000
        else
            integrator.u[2] = maximum([integrator.u[2] + noise[count], 0.0]) #https://discourse.julialang.org/t/change-parameters-in-a-model-in-specific-time-with-differentialequations/36930
        end
    end

    cb = PeriodicCallback(pert_cb2, freq, initial_affect = false)
    prob = ODEProblem(roz_mac_II!, u0, tspan, par)
    sol = DifferentialEquations.solve(prob, callback = cb, reltol = 1e-8)
    return solend = sol(tvals)
end
