include("packages.jl")
include("bacteriophage_commoncode.jl")
include("bacteriophage_gutresponses.jl")
include("bacteriophage_stability.jl")

#Figure 2 (Stability Schematic)
#Set up of equilibrium and solution curves
let
    u0=[0.5]
    tsend = 10000.0
    freq = 0.1
    tspan=(0.0, tsend)
    par = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.003)
    prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
    sol = solve(prob, RadauIIA5())
    solseries = sol(tsend-16.0:freq:tsend)
    sel = [sel_sine(par, t) for t in solseries.t]
    opt = optimum(sel)
    trackoptimumdata = trackoptimum(solseries[1,:], opt)
    println(mean(trackoptimumdata))
    println(std(trackoptimumdata))
    stabilitysetup = figure(figsize=(4,8))
    subplot(4,1,1)
    plot(solseries.t, sel, color="black", linewidth="3")
    ylabel("s(t)", fontsize=15)
    xlabel("Time", fontsize=15)
    hlines(0.0, minimum(solseries.t), maximum(solseries.t), color="black")
    ylim(-0.5,0.5)
    yticks([])
    xticks([])
    subplot(4,1,2)
    plot(solseries.t, opt, color="black", linestyle="dashed", linewidth="3")
    ylabel("Optimum C", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(-0.1, 1.1)
    yticks([0.0, 1.0], fontsize=12)
    xticks([])
    subplot(4,1,3)
    plot(solseries.t, opt, color="black", linestyle="dashed", linewidth="3")
    plot(solseries.t, solseries.u, color="black", linewidth="3")
    ylabel("Optimum C &\nC Solution", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(-0.1, 1.1)
    yticks([0.0, 1.0], fontsize=12)
    xticks([])
    subplot(4,1,4)
    plot(solseries.t, trackoptimumdata, color="black", linewidth="3")
    ylabel("Transitory Load", fontsize=15)
    xlabel("Time", fontsize=15)
    yticks([0.3, 0.5, 0.7], fontsize=12)
    xticks([])
    tight_layout()
    return stabilitysetup
    # savefig(joinpath(abpath(), "figs/stabilitysetup.pdf"))
end

#Figure - showing the different patterns of model
let 
    u0=[0.5]
    tsend = 10000.0
    freq = 0.1
    tspan=(0.0, tsend)
    par1 = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.03)
    prob1 = ODEProblem(bacphage_sine_forced!, u0, tspan, par1)
    sol1 = solve(prob1, RadauIIA5())
    solseries1 = sol1(tsend-20.0:freq:tsend)
    par2 = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.003)
    prob2 = ODEProblem(bacphage_sine_forced!, u0, tspan, par2)
    sol2 = solve(prob2, RadauIIA5())
    solseries2 = sol2(tsend-20.0:freq:tsend)
    par3 = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.001)
    prob3 = ODEProblem(bacphage_sine_forced!, u0, tspan, par3)
    sol3 = solve(prob3, RadauIIA5())
    solseries3 = sol3(tsend-20.0:freq:tsend)
    patternsfigure = figure(figsize = (7,2))
    subplot(1,3,1)
    plot(solseries1.t, solseries1.u, color="black")
    ylabel("C(t)", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    subplot(1,3,2)
    plot(solseries2.t, solseries2.u, color="black")
    ylabel("C(t)", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    subplot(1,3,3)
    plot(solseries3.t, solseries3.u, color="black")
    ylabel("C(t)", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    tight_layout()
    return patternsfigure
    # savefig(joinpath(abpath(), "figs/patternsfigure.pdf"))
end


#Figure

bifurcmid_data = bifurcmid(-0.01:0.0001:0.001, 100000.0)
white_noise_mid_data = bifurc_white_mid(-0.01:0.0001:0.001, 6, 100000.0)

let 
    par  = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=0.0)
    bifurcfigure = figure(figsize=(8,3))
    subplot(1,2,1)
    plot(bifurcmid_data[:, 1], bifurcmid_data[:, 2], color="black")
    plot(bifurcmid_data[:, 1], bifurcmid_data[:, 3], color="black")
    vlines(bifurc(par), 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Average selection", fontsize = 15)
    ylabel("C(t) min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.010, -0.002, 0.0])
    title("Sine Wave", fontsize = 15)
    subplot(1,2,2)
    plot(white_noise_mid_data[:, 1], white_noise_mid_data[:, 2], color="black")
    plot(white_noise_mid_data[:, 1], white_noise_mid_data[:, 3], color="black")
    vlines(bifurc(par), 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Average selection", fontsize = 15)
    ylabel("Mean C(t) min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.010, -0.002, 0.0])
    title("White Noise", fontsize = 15)
    tight_layout()
    return bifurcfigure
    # savefig(joinpath(abpath(), "figs/conjbacsel_balance_sinewhite.pdf"))
end



##Stability of Gut Function

#Final figure (tracking)

let 
    data001_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 0.1, -0.002, 0.01, 10000.0)
    data1_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 1, -0.002, 0.001, 10000.0)
    data10_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 10, -0.002, 0.001, 10000.0)
    data001_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 0.1, -0.002, 1.0, 10000.0, 100)
    data1_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 1, -0.002, 1.0, 10000.0, 100)
    data10_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 10, -0.002, 1.0, 10000.0, 100)
    rplusbconstrainedtrackingfigure = figure(figsize = (5,7))
    subplot(2,1,1)
    plot(data001_sine[:,1], data001_sine[:,2], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,2], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,2], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("CV of Transitory Load", fontsize = 15)
    xticks([0.0, 0.001, 0.002, 0.003, 0.004], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Sine Wave", fontsize = 15)
    subplot(2,1,2)
    plot(data001_noise[:,1], data001_noise[:,2], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,2], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,2], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("CV of Transitory Load", fontsize = 15)
    xticks([0.0, 0.001, 0.002, 0.003, 0.004], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("White Noise", fontsize = 15)
    tight_layout()
    return rplusbconstrainedtrackingfigure
    # savefig(joinpath(abpath(), "figs/rplusbconstrainedtrackingfigure.pdf"))
end

#Figure
let 
    data_increaseb = time_selectionswitch_b(0.0001:0.0001:0.01, -0.05, 0.01, 0.0:1.0:10000.0)
    data_decreaseb = time_selectionswitch_b(0.0001:0.0001:0.01, 0.01, -0.05, 0.0:1.0:10000.0)
    data_increaser = time_selectionswitch_r(0.0001:0.0001:0.01, -0.05, 0.01, 0.0:1.0:10000.0)
    data_decreaser = time_selectionswitch_r(0.0001:0.0001:0.01, 0.01, -0.05, 0.0:1.0:10000.0)
    delayfigure = figure(figsize=(9,7))
    plot(data_increaseb[1], data_increaseb[2], color="#73D055FF", linewidth = 3, label="Positive (b)")
    plot(data_decreaseb[1], data_decreaseb[2], color="#440154FF", linewidth = 3, label="Negative (b)")
    plot(data_increaser[1], data_increaser[2], color="#73D055FF", linewidth = 3, label="Positive (r)", linestyle="dashed")
    plot(data_decreaser[1], data_decreaser[2], color="#440154FF", linewidth = 3, label="Negative (r)", linestyle="dashed")
    xticks(fontsize=12)
    yticks(fontsize=12)
    xlabel("Bacteriophage level (\$b\$)\nConjugation level (\$r\$)", fontsize = 15)
    ylabel("Return Time", fontsize = 15)
    ylim(0.0, 620)
    legend(title = "Selection Switch", title_fontsize = 15, fontsize = 12)
    return delayfigure
    # savefig(joinpath(abpath(), "figs/delay_selectionswitch.pdf"))
end


#Supporting Information
include("bacteriophage_supportinginformation.jl")

#Mean transitory load
let 
    data001_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 0.1, -0.002, 0.1, 10000.0)
    data1_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 1, -0.002, 0.1, 10000.0)
    data10_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 10, -0.002, 0.1, 10000.0)
    data001_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 0.1, -0.002, 1.0, 10000.0, 100)
    data1_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 1, -0.002, 1.0, 10000.0, 100)
    data10_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 10, -0.002, 1.0, 10000.0, 100)
    rplusbconstrainedtrackingfigure = figure(figsize = (5,7))
    subplot(2,1,1)
    plot(data001_sine[:,1], data001_sine[:,3], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,3], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,3], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Mean Transitory Load", fontsize = 15)
    xticks([0.0, 0.001, 0.002, 0.003, 0.004], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Sine Wave", fontsize = 15)
    subplot(2,1,2)
    plot(data001_noise[:,1], data001_noise[:,3], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,3], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,3], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Mean Transitory Load", fontsize = 15)
    xticks([0.0, 0.001, 0.002, 0.003, 0.004], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("White Noise", fontsize = 15)
    tight_layout()
    # return rplusbconstrainedtrackingfigure
    savefig(joinpath(abpath(), "figs/rplusbconstrainedtrackingfigure_meantransitoryload.pdf"))
end

#understanding why mean transitory load is between 0.4990 and 0.5010
function brconstrained_stabilitytracking_sine_test(bplusr, brratio, smid, freq, tsend)
    u0=[0.5]
    tspan=(0.0, tsend)
    rval = bplusr/(1+brratio)
    bval = bplusr - rval
    par = BacPhageSineForcedPar(b = bval, r=rval, per=0.5, amp=0.4, mid=smid)
    prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
    sol = solve(prob, RadauIIA5())
    solseries = sol(tsend-1000.0:freq:tsend)
    selectiondata = [sel_sine(par, t) for t in tsend-1000.0:freq:tsend]
    seriesoptimum = optimum(selectiondata)
    optimumdiff = trackoptimum(solseries[1,:], seriesoptimum)
    return optimumdiff
end

test = brconstrained_stabilitytracking_sine_test(0.003, 0.1, -0.002, 0.1, 10000.0)
length(test)

let 
    data1= brconstrained_stabilitytracking_sine_test(0.003, 0.001, -0.002, 0.1, 10000.0)
    println(mean(data1))
    data2= brconstrained_stabilitytracking_sine_test(0.0005, 0.001, -0.002, 0.1, 10000.0)
    println(mean(data2))
    test = figure()
    plot(1.0:1.0:Float64(length(data1)), data1)
    plot(1.0:1.0:Float64(length(data2)), data2)
    xlim(9800, 10000)
    return test 
end

#Total transitory load
let 
    data001_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 0.1, -0.002, 0.1, 10000.0)
    data1_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 1, -0.002, 0.1, 10000.0)
    data10_sine = brconstrained_stabilitytracking_sine(0.00001:0.0001:0.004, 10, -0.002, 0.1, 10000.0)
    data001_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 0.1, -0.002, 1.0, 10000.0, 100)
    data1_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 1, -0.002, 1.0, 10000.0, 100)
    data10_noise = brconstrained_stabilitytracking_noise(0.00001:0.0001:0.004, 10, -0.002, 1.0, 10000.0, 100)
    rplusbconstrainedtrackingfigure = figure(figsize = (5,7))
    subplot(2,1,1)
    plot(data001_sine[:,1], data001_sine[:,4], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,4], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,4], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Total Transitory Load", fontsize = 15)
    xticks([0.0, 0.001, 0.002, 0.003, 0.004], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Sine Wave", fontsize = 15)
    subplot(2,1,2)
    plot(data001_noise[:,1], data001_noise[:,4], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,4], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,4], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Total Transitory Load", fontsize = 15)
    xticks([0.0, 0.001, 0.002, 0.003, 0.004], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("White Noise", fontsize = 15)
    tight_layout()
    # return rplusbconstrainedtrackingfigure
    savefig(joinpath(abpath(), "figs/rplusbconstrainedtrackingfigure.pdf"))
end

#Bifurcation analysis
let 
    bifurcval = bifurc(BacPhagePar())
    st = -0.3
    en = 0.3 
    srange1 = st:0.0001:bifurcval
    srange2 = bifurcval:0.0001:en
    data1 = [interior_equil(s, BacPhagePar()) for s in srange1]
    data2 = [interior_equil(s, BacPhagePar()) for s in srange2]
    bifurcwlrplot = figure(figsize=(6,5))
    plot(srange1, data1, color = "black")
    plot(srange2, data2, linestyle= "dashed",color = "black")
    ylabel("Ĉ", fontsize = 15)
    xlabel("s", fontsize = 15)
    xlim(-0.3, 0.3)
    ylim(-0.1, 1.1)
    xticks(fontsize = 12)
    yticks(fontsize = 12)
    hlines(1.0, st, bifurcval, linestyle= "dashed", colors= "black")
    hlines(1.0, bifurcval, en, colors= "black")
    # return bifurcwlrplot
    savefig(joinpath(abpath(), "figs/SIbifurcfigure.pdf"))
end

#Eigen integral balancing
bifurcintegral_data = bifurcintegral_eigen1_mid(-0.01:0.0001:0.001, 100000.0)

let 
    srange = -0.02:0.01:0.02
    par  = BacPhageSineForcedPar(b = 0.001)
    eigenb001 = [eigen1(s, par) for s in srange]
    bifurc_integralfigure = figure(figsize=(8,8))
    subplot(1,2,1)
    plot(srange, eigenb001, color="blue", linewidth=3)
    hlines(0.0, -0.02, 0.02, linewidth=0.5)
    vlines(bifurc(par), -0.02, 0.02, linewidth=0.5)
    xlabel("s")
    ylabel("λ (Ĉ=1)")
    subplot(1,2,2)
    plot(bifurcintegral_data[:, 1], bifurcintegral_data[:, 2], color="black")
    plot(bifurcintegral_data[:, 1], bifurcintegral_data[:, 3], color="black")
    vlines(0.0, 0.0, 1.0, linestyles="dashed", color="black")
    xlabel("λ Integral")
    ylabel("C(t) min & max")
    tight_layout()
    return bifurc_integralfigure
    # savefig(joinpath(abpath(), "figs/conjbacsel_balance.pdf"))
end

#Periodicity
bifurcmid_data_lowper = bifurcmid_per(-0.01:0.0001:0.001, 0.3, 100000.0)
bifurcmid_data_normper = bifurcmid_per(-0.01:0.0001:0.001, 0.5, 100000.0)
bifurcmid_data_highper = bifurcmid_per(-0.01:0.0001:0.001, 0.7, 100000.0)
let 
    periodicityplot = figure(figsize=(5,8))
    subplot(3,1,1)
    plot(bifurcmid_data_lowper[:, 1], bifurcmid_data_lowper[:, 2], color="black")
    plot(bifurcmid_data_lowper[:, 1], bifurcmid_data_lowper[:, 3], color="black")
    vlines(-0.002, 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Average selection", fontsize = 15)
    ylabel("C(t) min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.010, -0.002, 0.0], fontsize = 12)
    subplot(3,1,2)
    plot(bifurcmid_data_normper[:, 1], bifurcmid_data_normper[:, 2], color="black")
    plot(bifurcmid_data_normper[:, 1], bifurcmid_data_normper[:, 3], color="black")
    vlines(-0.002, 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Average selection", fontsize = 15)
    ylabel("C(t) min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.010, -0.002, 0.0], fontsize = 12)
    subplot(3,1,3)
    plot(bifurcmid_data_highper[:, 1], bifurcmid_data_highper[:, 2], color="black")
    plot(bifurcmid_data_highper[:, 1], bifurcmid_data_highper[:, 3], color="black")
    vlines(-0.002, 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Average selection", fontsize = 15)
    ylabel("C(t) min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.010, -0.002, 0.0], fontsize = 12)
    tight_layout()
    # return periodicityplot
    savefig(joinpath(abpath(), "figs/SIperiodicityplot.pdf"))
end
