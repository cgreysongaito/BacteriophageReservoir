include("packages.jl")
include("bacteriophage_commoncode.jl")



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
    tspan=(0.0, 300.0)

    condition(u,t,integrator) = t > 100.0 
    function changesel!(integrator)
        integrator.p.s=0.01
    end

    cb = DiscreteCallback(condition, changesel!)
    prob = ODEProblem(bacphage!, u0, tspan, BacPhagePar(s=-0.25))
    sol = solve(prob, callback=cb)
    shortsol = sol(100.0:0.5:200.0)
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
    tspan=(0.0, 300.0)

    condition(u,t,integrator) = t > 100.0 
    function changesel!(integrator)
        integrator.p.s=0.01
    end

    cb = DiscreteCallback(condition, changesel!)
    prob = ODEProblem(bacphage_wobac!, u0, tspan, BacPhagePar(s=-0.25))
    sol = solve(prob, callback=cb)
    shortsol = sol(100.0:0.5:200.0)
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