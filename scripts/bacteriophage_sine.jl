include("packages.jl")
include("bacteriophage_commoncode.jl")


let 
    u0 = [0.5, 0.0]
    tspan=(0.0, 300.0)
    prob = ODEProblem(bacphage_sine!, u0, tspan, BacPhageSinePar(per = 0.11))
    sol = solve(prob)

    test = figure()
    plot(sol.t, sol.u)
    return test
end

#problem with model if selection long enough will just go to 1 (all carrier) and stay there regardless of selection changing (equilibrium)


let 
    u0 = [0.5, 0.0]
    tspan=(0.0, 300.0)
    prob = ODEProblem(bacphage!, [0.5], (0.0, 300.0),  BacPhagePar())
    sol = solve(prob)

    test = figure()
    plot(sol.t, sol.u)
    return test
end

prob = ODEProblem(bacphage!, [0.5], (0.0, 300.0),  BacPhagePar())
integrator = init(prob)