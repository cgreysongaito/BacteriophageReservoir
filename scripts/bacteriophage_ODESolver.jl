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
    solseries = sol(0.0:1.0:300.0)
    println(sol.alg)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    ylim(-1,1.6)
    return test
end

let #Rodas5
    par = BacPhageSineForcedPar()
    par.b = 0.01
    par.per = 0.2
    par.amp= 0.5
    par.mid = +(0.11-0.099099099)
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par)
    sol = solve(prob, Rodas5())
    solseries = sol(0.0:1.0:300.0)
    println(sol.alg)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    ylim(-1,1.6)
    return test
end


let #Rodas5
    par = BacPhageSineForcedPar()
    par.b = 0.01
    par.per = 0.1
    par.amp= 0.5
    par.mid = -0.1
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver2!, u0, tspan, par)
    sol = solve(prob, Rodas5())
    solseries = sol(0.0:1.0:300.0)
    println(sol.alg)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    ylim(-1,1.6)
    return test
end

#why does version one fluctuate consistently whereas version two does not fluctuate


bifurc_ver1(BacPhageSineForcedPar())
bifurc_ver2(BacPhageSineForcedPar())

#Comparing whether c goes above 1 in each version
let #Rodas5
    par = BacPhageSineForcedPar()
    par.b = 0.01
    par.per = 0.2
    par.amp= 0.5
    par.mid = +(0.11-0.099099099)
    u0 = [0.5]
    tspan=(0.0, 500.0)
    prob = ODEProblem(bacphage_sine_forced_ver1!, u0, tspan, par)
    sol = solve(prob, Rodas5())
    solseries = sol(40.0:1.0:50.0)
    return solseries
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
    solseries = sol(50.0:1.0:100.0)
    return solseries
end
# version 2 does not go above 1

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
    data_ver1b = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.04)) for s in srange]
    data_ver1c = [eigen1_ver1(s, BacPhageSineForcedPar(b=0.08)) for s in srange]
    data_ver2 = [eigen1_ver2(s, BacPhageSineForcedPar()) for s in srange]
    test = figure()
    plot(srange, data_ver1a, color = "blue")
    plot(srange, data_ver1b, color = "green")
    plot(srange, data_ver1c, color = "purple")
    plot(srange, data_ver2, color = "red")
    xlabel("s")
    ylabel("λ")
    return test
end
#b just shifts where the eigenvalue line sits but does not flatten the graph

let 
    srange = -0.5:0.01:1.00
    data_ver1a = [eigen1_ver1(s, BacPhageSineForcedPar()) for s in srange]
    data_ver1b = [eigen1_ver1(s, BacPhageSineForcedPar(r=0.2)) for s in srange]
    data_ver1c = [eigen1_ver1(s, BacPhageSineForcedPar(r=0.3)) for s in srange]
    data_ver2 = [eigen1_ver2(s, BacPhageSineForcedPar()) for s in srange]
    test = figure()
    plot(srange, data_ver1a, color = "blue")
    plot(srange, data_ver1b, color = "green")
    plot(srange, data_ver1c, color = "purple")
    plot(srange, data_ver2, color = "red")
    xlabel("s")
    ylabel("λ")
    return test
end
#r also doesn't flatten the eigenvalue line, just shifts the line


#why does version 2 go to fixation or 0 (can it ever find oscilating "attractor")?
#version 1 can go to fixation just different point where this happens - so 
#next step figure out way to find bifurcation values of different qualitative behaviours for version 1 and version 2
#maybe "bifrucation" is dependent of delay of system state matching equilibrium (i.e if we make the delay longer then quicker to fixation OR if we flatten the equilibria line then also quicker to fixation)
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
    solseries = sol(50.0:1.0:500.0)
    test = figure()
    plot(solseries.t, solseries.u)
    plot(solseries.t, [sel_sine(par, t) for t in solseries.t])
    hlines(bifurc_ver1(par), 50.0, 500.0, colors="green" )
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
