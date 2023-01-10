include("packages.jl")
include("bacteriophage_commoncode.jl")
include("bacteriophage_buffer.jl")
include("bacteriophage_stochastic.jl")

#Figure
let 
    data = time_fixation_b(0.0001:0.0001:0.01, -0.05, 0.01, 0.0:1.0:10000.0)
    delayfigure = figure()
    plot(data[1], data[2])
    xlabel("b")
    ylabel("Time")
    # return delayfigure
    savefig(joinpath(abpath(), "figs/delay_selectionswitch.png"))
end

#Figure

bifurcmid_data = bifurcmid(-0.01:0.0001:0.001, 100000.0)
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


#Figure
white_noise_mid_data = bifurc_white_mid(-0.01:0.0001:0.001, 6, 100000.0)

let 
    par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=0.0)
    midrange = -0.01:0.0001:0.00
    white_noise = figure()
    plot(white_noise_mid_data[:, 1], white_noise_mid_data[:, 2], color="black")
    plot(white_noise_mid_data[:, 1], white_noise_mid_data[:, 3], color="black")
    vlines(bifurc(par), 0.0, 1.0, linestyles="dashed", color="black") #bifurc value or r+b=-mid
    xlabel("Selection midpoint")
    ylabel("Mean C min/max")
    # return white_noise
    savefig(joinpath(abpath(), "figs/whitenoise_mid_figure.pdf"))
end

red_noise_midabove = bifurc_red_mid(0.0:0.1:0.9, 0.0, 6, 100000.0)
red_noise_midcentred = bifurc_red_mid(0.0:0.1:0.9, -0.002, 6, 100000.0)
red_noise_midbelow = bifurc_red_mid(0.0:0.1:0.9, -0.004, 6, 100000.0)

let 
    par = BacPhageSineForcedPar(b = 0.001, per=0.5, amp=0.4, mid=0.0)
    corrrange = 0.0:0.1:0.9
    red_noise = figure(figsize=(3.5,5))
    subplot(3,1,1)
    plot(red_noise_midabove[:, 1], red_noise_midabove[:, 2], color="black")
    plot(red_noise_midabove[:, 1], red_noise_midabove[:, 3], color="black")
    xlabel("Noise correlation")
    ylabel("Mean C min/max")
    ylim(0,1.1)
    subplot(3,1,2)
    plot(red_noise_midcentred[:, 1], red_noise_midcentred[:, 2], color="black")
    plot(red_noise_midcentred[:, 1], red_noise_midcentred[:, 3], color="black")
    xlabel("Noise correlation")
    ylabel("Mean C min/max")
    ylim(0,1.1)
    subplot(3,1,3)
    plot(red_noise_midbelow[:, 1], red_noise_midbelow[:, 2], color="black")
    plot(red_noise_midbelow[:, 1], red_noise_midbelow[:, 3], color="black")
    xlabel("Noise correlation")
    ylabel("Mean C min/max")
    ylim(0,1.1)
    tight_layout()
    # return red_noise
    savefig(joinpath(abpath(), "figs/rednoise_mid_figure.pdf"))
end