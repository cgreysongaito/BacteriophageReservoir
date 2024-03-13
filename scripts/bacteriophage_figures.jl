include("packages.jl")
include("bacteriophage_commoncode.jl")
include("bacteriophage_stability.jl")

#Figure 2 (Stability Schematic)
#Set up of equilibrium and solution curves
let
    u0=[0.5]
    tsend = 10000.0
    freq = 0.1
    tspan=(0.0, tsend)
    par = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.003) #keeping amp=0.4 because conceptual figure and easier to see
    prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
    sol = solve(prob, Rodas4P())
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

##Stability of Gut Function

#Figure 3 Mean transitory load
let 
    data001_sine = brconstrained_stabilitytracking_sine(0.00001:0.00005:0.001, 0.1, -0.0005, 0.001, 2000.0)
    data1_sine = brconstrained_stabilitytracking_sine(0.00001:0.00005:0.001, 1, -0.0005, 0.001, 2000.0)
    data10_sine = brconstrained_stabilitytracking_sine(0.00001:0.00005:0.001, 10, -0.0005, 0.001, 2000.0)
    data001_noise = brconstrained_stabilitytracking_noise(0.00001:0.00005:0.001, 0.1, -0.0005, 1.0, 10000.0, 100)
    data1_noise = brconstrained_stabilitytracking_noise(0.00001:0.00005:0.001, 1, -0.0005, 1.0, 10000.0, 100)
    data10_noise = brconstrained_stabilitytracking_noise(0.00001:0.00005:0.001, 10, -0.0005, 1.0, 10000.0, 100)
    rplusbconstrainedtrackingfigure = figure(figsize = (5,7))
    subplot(2,1,1)
    plot(data001_sine[:,1], data001_sine[:,2], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,2], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,2], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Mean Transitory Load", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Sine Wave", fontsize = 15)
    subplot(2,1,2)
    plot(data001_noise[:,1], data001_noise[:,2], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,2], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,2], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Mean Transitory Load", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("White Noise", fontsize = 15)
    tight_layout()
    # return rplusbconstrainedtrackingfigure
    savefig(joinpath(abpath(), "figs/rplusbconstrainedtrackingfigure_meantransitoryload.pdf"))
end

#Figure 4 Fluctuations of transitory load
let 
    data001_sine = brconstrained_stabilitytracking_sine(0.00001:0.00001:0.001, 0.1, -0.0005, 0.001, 500.0)
    data1_sine = brconstrained_stabilitytracking_sine(0.00001:0.00001:0.001, 1, -0.0005, 0.001, 500.0)
    data10_sine = brconstrained_stabilitytracking_sine(0.00001:0.00001:0.001, 10, -0.0005, 0.001, 500.0)
    data001_noise = brconstrained_stabilitytracking_noise(0.00001:0.00001:0.001, 0.1, -0.0005, 1.0, 10000.0, 100)
    data1_noise = brconstrained_stabilitytracking_noise(0.00001:0.00001:0.001, 1, -0.0005, 1.0, 10000.0, 100)
    data10_noise = brconstrained_stabilitytracking_noise(0.00001:0.00001:0.001, 10, -0.0005, 1.0, 10000.0, 100)
    rplusbconstrainedtrackingfigure = figure(figsize = (5,7))
    subplot(2,1,1)
    plot(data001_sine[:,1], data001_sine[:,3], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,3], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,3], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Transitory Load Fluctuation", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Sine Wave", fontsize = 15)
    subplot(2,1,2)
    plot(data001_noise[:,1], data001_noise[:,3], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,3], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,3], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Transitory Load Fluctuation", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("White Noise", fontsize = 15)
    tight_layout()
    # return rplusbconstrainedtrackingfigure
    savefig(joinpath(abpath(), "figs/rplusbconstrainedtrackingfigure_TLfluc.pdf"))
end

#Figure 5
let 
    data_increaseb = time_selectionswitch_b(0.00001:0.000005:0.001, -0.05, 0.01, 0.0:1.0:10000.0)
    data_decreaseb = time_selectionswitch_b(0.00001:0.000005:0.001, 0.01, -0.05, 0.0:1.0:10000.0)
    data_increaser = time_selectionswitch_r(0.00001:0.000005:0.001, -0.05, 0.01, 0.0:1.0:10000.0)
    data_decreaser = time_selectionswitch_r(0.00001:0.000005:0.001, 0.01, -0.05, 0.0:1.0:10000.0)
    TLmintime = figure(figsize=(8.5,7))
    plot(data_increaseb[1], data_increaseb[2], color="#73D055FF", linewidth = 3, label="Positive (b)")
    plot(data_decreaseb[1], data_decreaseb[2], color="#440154FF", linewidth = 3, label="Negative (b)")
    plot(data_increaser[1], data_increaser[2], color="#73D055FF", linewidth = 3, label="Positive (r)", linestyle="dashed")
    plot(data_decreaser[1], data_decreaser[2], color="#440154FF", linewidth = 3, label="Negative (r)", linestyle="dashed")
    xticks(fontsize=12)
    yticks(fontsize=12)
    xlabel("Bacteriophage rate (\$b\$)\nConjugation rate (\$r\$)", fontsize = 15)
    ylabel("Transitory load minimization time", fontsize = 15)
    ylim(0.0, 1000)
    legend(title = "Optimum C Switch", title_fontsize = 15, fontsize = 12)
    # return TLmintime
    savefig(joinpath(abpath(), "figs/TLmintime_selectionswitch.pdf"))
end


#Supporting Information
include("bacteriophage_supportinginformation.jl")

range_parameter_br(0.1, 0.00001:0.0001:0.001)
range_parameter_br(1.0, 0.00001:0.0001:0.001)
range_parameter_br(10.0, 0.00001:0.0001:0.001)

0.1*sin(π/2)-0.0005
-0.1*sin(π/2)-0.0005

sel_noise_range(-0.0005, 0.015, 100000) #setting noise standard deviation to 0.015

let 
    bplusrrange = 0.00001:0.00001:0.001
    avabssel = -0.0005
    alpha001data = sel_rel_conj(0.1, bplusrrange, avabssel)
    alpha1data = sel_rel_conj(1.0, bplusrrange, avabssel)
    alpha10data = sel_rel_conj(10.0, bplusrrange, avabssel)
    alpha_relsel_figure = figure()
    plot(alpha001data[:,1],alpha001data[:,2], color="#FDE725FF", label="b/r=0.1")
    plot(alpha1data[:,1],alpha1data[:,2],  color="#29AF7FFF", label="b/r=1")
    plot(alpha10data[:,1],alpha10data[:,2], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    ylabel("α", fontsize = 15)
    xticks(fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    # return alpha_relsel_figure
    savefig(joinpath(abpath(), "figs/alpha_relsel_figure.pdf"))
end

#Balance of selection, conjugation, and bacteriophage transduction

#Bifurcation analysis of varying s
let 
    bifurcval = bifurc(BacPhagePar())
    st = -0.1
    en = 0.1
    srange1 = st:0.0001:bifurcval
    srange2 = bifurcval:0.0001:en
    data1 = [interior_equil(s, BacPhagePar()) for s in srange1]
    data2 = [interior_equil(s, BacPhagePar()) for s in srange2]
    bifurcwlrplot = figure(figsize=(6,5))
    plot(srange1, data1, color = "black")
    plot(srange2, data2, linestyle= "dashed",color = "black")
    ylabel("Ĉ", fontsize = 15)
    xlabel("s", fontsize = 15)
    xlim(-0.1, 0.1)
    ylim(-0.1, 1.1)
    xticks([-0.1,-0.05,0.0,0.05,0.1],fontsize = 12)
    yticks(fontsize = 12)
    hlines(1.0, st, bifurcval, linestyle= "dashed", colors= "black")
    hlines(1.0, bifurcval, en, colors= "black")
    # return bifurcwlrplot
    savefig(joinpath(abpath(), "figs/SIbifurcfigure.pdf"))
end

#Gut responses to environmental variation
let 
    u0=[0.5]
    tsend = 10000.0
    freq = 0.1
    tspan=(0.0, tsend)
    par1 = BacPhageSineForcedPar(b = 0.0001, r=0.0001, per=0.5, amp=0.1, mid=-0.003)
    prob1 = ODEProblem(bacphage_sine_forced!, u0, tspan, par1)
    sol1 = solve(prob1, Rodas4P())
    solseries1 = sol1(tsend-20.0:freq:tsend)
    par2 = BacPhageSineForcedPar(b = 0.0001, r=0.0001, per=0.5, amp=0.1, mid=-0.0003)
    prob2 = ODEProblem(bacphage_sine_forced!, u0, tspan, par2)
    sol2 = solve(prob2, Rodas4P())
    solseries2 = sol2(tsend-20.0:freq:tsend)
    par3 = BacPhageSineForcedPar(b = 0.0001, r=0.0001, per=0.5, amp=0.1, mid=0.01)
    prob3 = ODEProblem(bacphage_sine_forced!, u0, tspan, par3)
    sol3 = solve(prob3, Rodas4P())
    solseries3 = sol3(tsend-20.0:freq:tsend)
    patternsfigure = figure(figsize = (7,2))
    subplot(1,3,1)
    plot(solseries1.t, solseries1.u, color="black")
    ylabel("C(t)", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    xticks([])
    subplot(1,3,2)
    plot(solseries2.t, solseries2.u, color="black")
    ylabel("C(t)", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    xticks([])
    subplot(1,3,3)
    plot(solseries3.t, solseries3.u, color="black")
    ylabel("C(t)", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    xticks([])
    tight_layout()
    # return patternsfigure
    savefig(joinpath(abpath(), "figs/patternsfigure.pdf"))
end

bifurcmid_data = bifurcmid(-0.001:0.00001:0.0001, 100000.0)
white_noise_mid_data = bifurc_white_mid(-0.001:0.00001:0.0001, 0.015, 6, 100000.0)

let 
    par  = BacPhageSineForcedPar(b = 0.0001, r=0.0001, per=0.5, amp=0.1, mid=0.0)
    bifurcfigure = figure(figsize=(8,3))
    subplot(1,2,1)
    plot(bifurcmid_data[:, 1], bifurcmid_data[:, 2], color="black")
    plot(bifurcmid_data[:, 1], bifurcmid_data[:, 3], color="black")
    vlines(bifurc(par), 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Average selection", fontsize = 15)
    ylabel("C(t) min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.001, -0.0002, 0.0])
    title("Sine Wave", fontsize = 15)
    subplot(1,2,2)
    plot(white_noise_mid_data[:, 1], white_noise_mid_data[:, 2], color="black")
    plot(white_noise_mid_data[:, 1], white_noise_mid_data[:, 3], color="black")
    vlines(bifurc(par), 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Average selection", fontsize = 15)
    ylabel("Mean C(t) min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.001, -0.0002, 0.0])
    title("White Noise", fontsize = 15)
    tight_layout()
    # return bifurcfigure
    savefig(joinpath(abpath(), "figs/conjbacsel_balance_sinewhite.pdf"))
end

#SI Figure 4 - Mean Transitory Load and Transitory Load split between 0 and 1 Optimum C (Sine)
let 
    data001_sine = brconstrained_stabilitytracking_sine_splitoptimum(0.00001:0.00001:0.001, 0.1, -0.0005, 0.001, 500.0)
    data1_sine = brconstrained_stabilitytracking_sine_splitoptimum(0.00001:0.00001:0.001, 1, -0.0005, 0.001, 500.0)
    data10_sine = brconstrained_stabilitytracking_sine_splitoptimum(0.00001:0.00001:0.001, 10, -0.0005, 0.001, 500.0)
    stability_splitoptimum_sine = figure(figsize = (8,7))
    subplot(2,2,1)
    plot(data001_sine[:,1], data001_sine[:,2], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,2], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,2], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Mean Transitory Load", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Optimum C = 0", fontsize = 15)
    subplot(2,2,2)
    plot(data001_sine[:,1], data001_sine[:,3], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,3], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,3], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Transitory Load Fluctuation", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Optimum C = 0", fontsize = 15)
    subplot(2,2,3)
    plot(data001_sine[:,1], data001_sine[:,4], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,4], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,4], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Mean Transitory Load", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Optimum C = 1", fontsize = 15)
    subplot(2,2,4)
    plot(data001_sine[:,1], data001_sine[:,5], color="#FDE725FF", label="b/r=0.1")
    plot(data1_sine[:,1], data1_sine[:,5], color="#29AF7FFF", label="b/r=1")
    plot(data10_sine[:,1], data10_sine[:,5], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Transitory Load Fluctuation", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Optimum C = 1", fontsize = 15)
    tight_layout()
    # return stability_splitoptimum_sine
    savefig(joinpath(abpath(), "figs/stability_splitoptimum_sine.pdf"))
end

#SI Figure 5 - Mean Transitory Load and Transitory Load split between 0 and 1 Optimum C (Noise)
let 
    data001_noise = brconstrained_stabilitytracking_noise_splitoptimum(0.00001:0.00001:0.001, 0.1, -0.0005, 1.0, 10000.0, 100)
    data1_noise = brconstrained_stabilitytracking_noise_splitoptimum(0.00001:0.00001:0.001, 1.0, -0.0005, 1.0, 10000.0, 100)
    data10_noise = brconstrained_stabilitytracking_noise_splitoptimum(0.00001:0.00001:0.001, 10.0, -0.0005, 1.0, 10000.0, 100)
    stability_splitoptimum_noise = figure(figsize = (8,7))
    subplot(2,2,1)
    plot(data001_noise[:,1], data001_noise[:,2], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,2], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,2], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Mean Transitory Load", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Optimum C = 0", fontsize = 15)
    subplot(2,2,2)
    plot(data001_noise[:,1], data001_noise[:,3], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,3], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,3], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Transitory Load Fluctuation", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Optimum C = 0", fontsize = 15)
    subplot(2,2,3)
    plot(data001_noise[:,1], data001_noise[:,4], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,4], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,4], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Mean Transitory Load", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Optimum C = 1", fontsize = 15)
    subplot(2,2,4)
    plot(data001_noise[:,1], data001_noise[:,5], color="#FDE725FF", label="b/r=0.1")
    plot(data1_noise[:,1], data1_noise[:,5], color="#29AF7FFF", label="b/r=1")
    plot(data10_noise[:,1], data10_noise[:,5], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Transitory Load Fluctuation", fontsize = 15)
    xticks([0.0, 0.0005, 0.001], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    title("Optimum C = 1", fontsize = 15)
    tight_layout()
    # return stability_splitoptimum_noise
    savefig(joinpath(abpath(), "figs/stability_splitoptimum_noise.pdf"))
end

#Lowest C analysis for Mean Transitory Load
let
    slow = 0.1*sin(1.5*pi)-.0005
    bplusrrange = 0.000001:0.000001:0.001
    data01 = meanTL_lowestC(bplusrrange, 0.1, slow)
    data1 = meanTL_lowestC(bplusrrange, 1, slow)
    data10 = meanTL_lowestC(bplusrrange, 10, slow)
    lowestCplot = figure(figsize=(7,5))
    plot(data01[:,1], data01[:,2], color="#FDE725FF", label="b/r=0.1")
    plot(data1[:,1], data1[:,2], color="#29AF7FFF", label="b/r=1")
    plot(data01[:,1], data10[:,2], color="#39568CFF", label="b/r=10")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Lowest Ĉ", fontsize = 15)
    xticks([0.0,0.0005,0.001],fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    # return lowestCplot
    savefig(joinpath(abpath(), "figs/lowestCplot.pdf"))
end

