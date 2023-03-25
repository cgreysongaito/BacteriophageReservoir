include("packages.jl")
include("bacteriophage_commoncode.jl")
include("bacteriophage_buffer.jl")
include("bacteriophage_stochastic.jl")


#Figure 2 (Stability Schematic)
#Set up of equilibrium and solution curves
let
    u0=[0.5]
    tsend = 10000.0
    freq = 0.1
    tspan=(0.0, tsend)
    par1 = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.04)
    prob1 = ODEProblem(bacphage_sine_forced!, u0, tspan, par1)
    sol1 = solve(prob1, RadauIIA5())
    solseries1 = sol1(tsend-15.0:freq:tsend)
    sel1 = [sel_sine(par1, t) for t in solseries1.t]
    equil1 = [stableequil(s, par1) for s in sel1]
    par2 = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.0025)
    prob2 = ODEProblem(bacphage_sine_forced!, u0, tspan, par2)
    sol2 = solve(prob2, RadauIIA5())
    solseries2 = sol2(tsend-15.0:freq:tsend)
    sel2 = [sel_sine(par2, t) for t in solseries2.t]
    equil2 = [stableequil(s, par2) for s in sel2]
    par3 = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.001)
    prob3 = ODEProblem(bacphage_sine_forced!, u0, tspan, par3)
    sol3 = solve(prob3, RadauIIA5())
    solseries3 = sol3(tsend-15.0:freq:tsend)
    sel3 = [sel_sine(par3, t) for t in solseries3.t]
    equil3 = [stableequil(s, par3) for s in sel3]
    trackingsetup = figure(figsize=(8,2.5))
    subplot(1,3,1)
    plot(solseries1.t, solseries1.u, color="black")
    plot(solseries1.t, equil1, color="black", linestyle="dashed")
    ylabel("C", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(-0.05,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    subplot(1,3,2)
    plot(solseries2.t, solseries2.u, color="black")
    plot(solseries2.t, equil2, color="black", linestyle="dashed")
    ylabel("C", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(-0.05,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    subplot(1,3,3)
    plot(solseries3.t, solseries3.u, color="black")
    plot(solseries3.t, equil3, color="black", linestyle="dashed")
    ylabel("C", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(-0.05,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    tight_layout()
    # return trackingsetup
    savefig(joinpath(abpath(), "figs/trackingsetup.pdf"))
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
    ylabel("C", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    subplot(1,3,2)
    plot(solseries2.t, solseries2.u, color="black")
    ylabel("C", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    subplot(1,3,3)
    plot(solseries3.t, solseries3.u, color="black")
    ylabel("C", fontsize=15)
    xlabel("Time", fontsize=15)
    ylim(0.0,1.05)
    yticks([0.0,0.5,1.0])
    tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    tight_layout()
    # return patternsfigure
    savefig(joinpath(abpath(), "figs/patternsfigure.pdf"))
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
    ylabel("C min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.010, -0.002, 0.0])
    title("Sine Wave", fontsize = 15)
    subplot(1,2,2)
    plot(white_noise_mid_data[:, 1], white_noise_mid_data[:, 2], color="black")
    plot(white_noise_mid_data[:, 1], white_noise_mid_data[:, 3], color="black")
    vlines(bifurc(par), 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Average selection", fontsize = 15)
    ylabel("Mean C min & max", fontsize = 15)
    yticks([0.0, 0.5,1.0], fontsize = 12)
    xticks([-0.010, -0.002, 0.0])
    title("White Noise", fontsize = 15)
    tight_layout()
    # return bifurcfigure
    savefig(joinpath(abpath(), "figs/conjbacsel_balance_sinewhite.pdf"))
end

red_noise_midbelow = bifurc_red_mid(0.0:0.1:0.9, -0.004, 6, 100000.0)
red_noise_midcentred = bifurc_red_mid(0.0:0.1:0.9, -0.002, 6, 100000.0)
red_noise_midabove = bifurc_red_mid(0.0:0.1:0.9, 0.0, 6, 100000.0)

let 
    par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=0.0)
    corrrange = 0.0:0.1:0.9
    red_noise = figure(figsize=(9,3))
    subplot(1,3,1)
    plot(red_noise_midbelow[:, 1], red_noise_midbelow[:, 2], color="black")
    plot(red_noise_midbelow[:, 1], red_noise_midbelow[:, 3], color="black")
    xlabel("Noise correlation", fontsize = 15)
    ylabel("Mean C min & max", fontsize = 15)
    ylim(0,1.05)
    yticks([0.0,0.5,1.0], fontsize = 12)
    xticks(fontsize = 12)
    subplot(1,3,2)
    plot(red_noise_midcentred[:, 1], red_noise_midcentred[:, 2], color="black")
    plot(red_noise_midcentred[:, 1], red_noise_midcentred[:, 3], color="black")
    xlabel("Noise correlation", fontsize = 15)
    ylabel("Mean C min & max", fontsize = 15)
    ylim(0,1.05)
    yticks([0.0,0.5,1.0], fontsize = 12)
    xticks(fontsize = 12)
    subplot(1,3,3)
    plot(red_noise_midabove[:, 1], red_noise_midabove[:, 2], color="black")
    plot(red_noise_midabove[:, 1], red_noise_midabove[:, 3], color="black")
    xlabel("Noise correlation", fontsize = 15)
    ylabel("Mean C min & max", fontsize = 15)
    ylim(0,1.05)
    yticks([0.0,0.5,1.0], fontsize = 12)
    xticks(fontsize = 12)
    tight_layout()
    # return red_noise
    savefig(joinpath(abpath(), "figs/rednoise_mid_figure.pdf"))
end

##Stability of Gut Function

#Figure
let 
    data_increaseb = time_fixation_b(0.0001:0.0001:0.01, -0.05, 0.01, 0.0:1.0:10000.0)
    data_decreaseb = time_fixation_b(0.0001:0.0001:0.01, 0.01, -0.05, 0.0:1.0:10000.0)
    data_increaser = time_fixation_r(0.0001:0.0001:0.01, -0.05, 0.01, 0.0:1.0:10000.0)
    data_decreaser = time_fixation_r(0.0001:0.0001:0.01, 0.01, -0.05, 0.0:1.0:10000.0)
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
    # return delayfigure
    savefig(joinpath(abpath(), "figs/delay_selectionswitch.pdf"))
end

#Final figure (tracking)

let 
    data05 = brconstrained_tracking(0.00001:0.0001:0.004, 0.01, -0.002, 0.1, 10000.0)
    data1 = brconstrained_tracking(0.00001:0.0001:0.004, 1, -0.002, 0.1, 10000.0)
    data2 = brconstrained_tracking(0.00001:0.0001:0.004, 3, -0.002, 0.1, 10000.0)
    rplusbconstrainedtrackingfigure = figure(figsize = (7,5))
    plot(data05[:,1], data05[:,2], color="#FDE725FF", label="r/b=0.01")
    plot(data1[:,1], data1[:,2], color="#29AF7FFF", label="r/b=1")
    plot(data2[:,1], data2[:,2], color="#39568CFF", label="r/b=3")
    xlabel("Horizontal Gene Transfer (\$b\$ + \$r\$)", fontsize = 15)
    ylabel("Log 10 absolute difference \nbetween solution and equilibrium", fontsize = 15)
    xticks([0.0, 0.001, 0.002, 0.003, 0.004], fontsize=12)
    yticks(fontsize=12)
    legend(fontsize = 12)
    return rplusbconstrainedtrackingfigure
    # savefig(joinpath(abpath(), "figs/rplusbconstrainedtrackingfigure.pdf"))
end


#Supporting Information

#Figure - Decomposition of conjugation, bacteriophage, and selection "work"
let 
    u0=[0.5]
    tsend = 10000.0
    freq = 0.1
    tspan=(0.0, tsend)
    par = BacPhageSineForcedPar(b = 0.001, r=0.001, per=0.5, amp=0.4, mid=-0.002)
    prob = ODEProblem(bacphage_sine_forced!, u0, tspan, par)
    sol = solve(prob, RadauIIA5())
    solseries = sol(tsend-100.0:freq:tsend)
    conjworkdata = [conjugation(C, par.r) for C in solseries[1, :]]
    lysoworkdata = [lysogeny(C, par.b) for C in solseries[1,:]]
    selectiondata = [sel_sine(par, t) for t in tsend-100.0:freq:tsend]
    selecworkdata = [selection(C, s) for (C,s) in zip(solseries[1,:],selectiondata)]
    workdecompositionfigure = figure()
    subplot(2,1,1)
    plot(solseries.t, lysoworkdata, color="red", label="Bacteriophage")
    plot(solseries.t, conjworkdata, color="blue", label="Conjugation")
    ylabel("Work (d/dt)")
    xlabel("Time")
    legend()
    subplot(2,1,2)
    plot(solseries.t, selecworkdata, color="green", label="Selection")
    ylabel("Work (d/dt)")
    xlabel("Time")
    legend()
    tight_layout()
    # return workdecompositionfigure
    savefig(joinpath(abpath(), "figs/workdecompositionfigure.pdf"))
end


#Figure - HGT and selection work
let 
    par = BacPhagePar(s = -0.002)
    Crange = 0.0:0.01:1.0
    conj = [conjugation(C, par.r) for C in Crange]
    bac = [lysogeny(C, par.b) for C in Crange]
    sel = [selection(C, par.s) for C in Crange]
    conj_bac = fillbetween_setup(conj, bac)
    conj_bac_sel = fillbetween_setup(conj_bac, sel)
    HGTsel_work = figure(figsize=(8,7))
    # plot(Crange, sel)
    fill_between(Crange, sel, color="#73D055FF")
    fill_between(Crange, conj, color="#440154FF")
    fill_between(Crange, conj, conj_bac, color="#404788FF")
    ylabel("Total HGT work")
    xlabel("C")
    title("s=-0.002")
    # fill_between(Crange, conj_bac, conj_bac_sel, color="#73D055FF")
    # return HGTsel_work
    savefig(joinpath(abpath(), "figs/HGTsel_work_002.pdf"))
end


#Eigen integral balancing
bifurcintegral_data = bifurcintegral_eigen1_mid(-0.01:0.0001:0.001, 100000.0)

let 
    srange = -0.02:0.01:0.02
    par  = BacPhageSineForcedPar(b = 0.001)
    eigenb001 = [eigen1(s, par) for s in srange]
    eigenb01 = [eigen1(s, BacPhageSineForcedPar(b = 0.01)) for s in srange]
    eigenb1 = [eigen1(s, BacPhageSineForcedPar(b = 0.1)) for s in srange]
    eigenr001 = [eigen1(s, BacPhageSineForcedPar(r = 0.001)) for s in srange]
    eigenr01 = [eigen1(s, BacPhageSineForcedPar(r = 0.01)) for s in srange]
    eigenr1 = [eigen1(s, BacPhageSineForcedPar(r = 0.1)) for s in srange]
    midrange = -0.01:0.0001:0.001
    # bifurcmid_data = bifurcmid(midrange, 100000.0)
    # bifurcintegral_data = bifurcintegral_eigen1_mid(midrange, 100000.0)
    bifurcfigure = figure(figsize=(8,8))
    subplot(3,2,1)
    plot(srange, eigenb001, color="blue", label = "b=0.001", linewidth=3)
    plot(srange, eigenb01, color="green", label = "b=0.01", linewidth=3)
    plot(srange, eigenb1, color="red", label = "b=0.1", linewidth=3)
    xlabel("s")
    ylabel("λ (Ĉ=1)")
    legend()
    subplot(3,2,2)
    plot(srange, eigenr001, color="blue", label = "r=0.001", linewidth=3)
    plot(srange, eigenr01, color="green", label = "r=0.01", linewidth=3)
    plot(srange, eigenr1, color="red", label = "r=0.1", linewidth=3)
    xlabel("s")
    ylabel("λ (Ĉ=1)")
    legend()
    subplot(3,2,3)
    plot(srange, eigenb001, color="blue", linewidth=3)
    hlines(0.0, -0.02, 0.02, linewidth=0.5)
    vlines(bifurc(par), -0.02, 0.02, linewidth=0.5)
    xlabel("s")
    ylabel("λ (Ĉ=1)")
    subplot(3,2,5)
    plot(bifurcmid_data[:, 1], bifurcmid_data[:, 2], color="black")
    plot(bifurcmid_data[:, 1], bifurcmid_data[:, 3], color="black")
    vlines(bifurc(par), 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Selection midpoint")
    ylabel("C min/max")
    subplot(3,2,6)
    plot(bifurcintegral_data[:, 1], bifurcintegral_data[:, 2], color="black")
    plot(bifurcintegral_data[:, 1], bifurcintegral_data[:, 3], color="black")
    vlines(0.0, 0.0, 1.0, linestyles="dashed", color="black")
    xlabel("λ Integral")
    ylabel("C min/max")
    tight_layout()
    # return bifurcfigure
    savefig(joinpath(abpath(), "figs/conjbacsel_balance.pdf"))
end
