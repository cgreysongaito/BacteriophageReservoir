include("packages.jl")
include("bacteriophage_commoncode.jl")

# Background for problem 1
# Here we are using version 1 of the model (with the denominator of the selection part as 1+s(t)C(t))

let
    par = BacPhageSineForcedPar()
    par.b = 0.02
    par.per = 0.1
    par.amp= 1.0
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par)
    sol = solve(prob)
    solseries = sol(0.0:1.0:50.0)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    ylim(-1,1.6)
    xlabel("Time (t)")
    ylabel("Ĉ & s")
    return test
end
#Blue line is C in time. Orange line is s in time.

#First check whether above 1 is due to degeneracy at s=-1

let 
    par = BacPhageSineForcedPar()
    par.b = 0.02
    par.per = 0.1
    par.amp= 1.0
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par)
    sol = solve(prob)
    solseries = sol(40.0:1.0:50.0)
    return solseries
end #WE ARE GETTING VALUES ABOVE 1


let 
    par = BacPhageSineForcedPar()
    par.b = 0.02
    par.per = 0.1
    par.amp= 1.0
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par)
    sol = solve(prob)
    solseries = sol(40.0:1.0:50.0)
    return [solseries[1,3],sel_sine(par, solseries.t[3])]
end #JUMP UP OF C is not to do with degeneracy at -1


#Why is C going above 1

#Next check change the solver
let  #default - auto so picks Rosenbrock23 but this has tolerances of 1e-2
    par = BacPhageSineForcedPar()
    par.b = 0.02
    par.per = 0.1
    par.amp= 1.0
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par)
    sol = solve(prob)
    return sol.alg
end 


let  #Rodas5
    par = BacPhageSineForcedPar()
    par.b = 0.02
    par.per = 0.1
    par.amp= 1.0
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par)
    sol = solve(prob, Rodas5())
    solseries = sol(0.0:1.0:300.0)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    ylim(-1,1.6)
    xlabel("Time (t)")
    ylabel("Ĉ & s")
    return test
end 
#Rodas5 looks to be the best (for the moment)

let
    par = BacPhageSinePar()
    par.b = 0.01
    par.per = 0.2
    par.amp= 0.5
    u0 = [0.5]
    tspan=(0.0, 500.0)
    # condition(u,t,integrator) = integrator.u[1] > 1.0
    # function returnC!(integrator)
    #     integrator.u[1] = 0.999999999999999
    # end

    # cb = DiscreteCallback(condition, returnC!)

    prob = ODEProblem(bacphage_sine_ver1!, u0, tspan, par)
    sol = solve(prob)
    solseries = sol(0.0:1.0:300.0)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    ylim(-1,1.6)
    return test
end




#Version 2 problems

#why is s C (1-C) settling at 1.0 even though s is cycling

let #default
    par = BacPhageSineForcedPar()
    par.b = 0.01
    par.per = 0.2
    par.amp= 0.5
    u0 = [0.5]
    tspan=(0.0, 500.0)
    # condition(u,t,integrator) = integrator.u[1] > 1.0
    # function returnC!(integrator)
    #     integrator.u[1] = 0.999999999999999
    # end

    # cb = DiscreteCallback(condition, returnC!)

    prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
    sol = solve(prob)
    solseries = sol(0.0:1.0:300.0)
    println(sol.alg)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    ylim(-1,1.6)
    return test
end

let #default
    par = BacPhageSineForcedPar()
    par.b = 0.01
    par.per = 0.2
    par.amp= 0.5
    u0 = [0.5]
    tspan=(0.0, 500.0)
    # condition(u,t,integrator) = integrator.u[1] > 1.0
    # function returnC!(integrator)
    #     integrator.u[1] = 0.999999999999999
    # end

    # cb = DiscreteCallback(condition, returnC!)

    prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
    sol = solve(prob, Rodas4P())
    solseries = sol(0.0:1.0:300.0)
    println(sol.alg)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    ylim(-1,1.6)
    return test
end



#why are we getting C values above 1
#could we fix the C values above 1 by using discretecallback to return to 1 if above 1 - (maybe 0.99999999)
