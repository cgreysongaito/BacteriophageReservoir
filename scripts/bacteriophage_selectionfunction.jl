include("packages.jl")
include("bacteriophage_commoncode.jl")

function selection_surface(srange, Crange)
    surface = zeros(length(srange),length(Crange))
    for i in 1:length(srange)
        for j in 1:length(Crange)
            surface[i,j] = srange[i] / (1 + srange[i] * Crange[j])
        end
    end
    return surface
end

let 
    n = 100
    srange = range(-0.9, stop=0.9, length=n)
    Crange = range(0.01, stop=1.0, length=n)
    xgrid = repeat(srange',n,1)
    ygrid = repeat(Crange,1,n)
    data = selection_surface(srange, Crange)
    surface_figure = figure()
    plot_surface(srange, Crange, data, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
    return surface_figure
end


function selection_line(srange, C)
    line = zeros(length(srange))
    for i in 1:length(srange)
        line[i] = srange[i] / (1 + srange[i] * C)
    end
    return line
end

let 
    n = 100
    srange = range(-0.9, stop=0.9, length=n)
    data1 = selection_line(srange, 0.01)
    data2 = selection_line(srange, 0.25)
    data3 = selection_line(srange, 0.5)
    data4 = selection_line(srange, 0.75)
    data5 = selection_line(srange, 1.0)
    lineplots = figure(figsize=(6,10))
    subplot(5,1,1)
    plot(srange, data1)
    xlabel("S")
    ylabel("Sel")
    ylim(-8, 1.0)
    title("C = 0.01")
    subplot(5,1,2)
    plot(srange, data2)
    xlabel("S")
    ylabel("Sel")
    ylim(-8, 1.0)
    title("C = 0.25")
    subplot(5,1,3)
    plot(srange, data3)
    xlabel("S")
    ylabel("Sel")
    ylim(-8, 1.0)
    title("C = 0.5")
    subplot(5,1,4)
    plot(srange, data4)
    xlabel("S")
    ylabel("Sel")
    ylim(-8, 1.0)
    title("C = 0.75")
    subplot(5,1,5)
    plot(srange, data5)
    xlabel("S")
    ylabel("Sel")
    ylim(-8, 1.0)
    title("C = 1.0")
    tight_layout()
    return lineplots
end

let 
    n = 100
    srange = range(-0.25, stop=0.25, length=n)
    data1 = selection_line(srange, 0.01)
    data2 = selection_line(srange, 0.25)
    data3 = selection_line(srange, 0.5)
    data4 = selection_line(srange, 0.75)
    data5 = selection_line(srange, 1.0)
    lineplots = figure(figsize=(6,10))
    subplot(5,1,1)
    plot(srange, data1)
    xlabel("S")
    ylabel("Sel")
    ylim(-0.4, 0.25)
    title("C = 0.01")
    subplot(5,1,2)
    plot(srange, data2)
    xlabel("S")
    ylabel("Sel")
    ylim(-0.4, 0.25)
    title("C = 0.25")
    subplot(5,1,3)
    plot(srange, data3)
    xlabel("S")
    ylabel("Sel")
    ylim(-0.4, 0.25)
    title("C = 0.5")
    subplot(5,1,4)
    plot(srange, data4)
    xlabel("S")
    ylabel("Sel")
    ylim(-0.4, 0.25)
    title("C = 0.75")
    subplot(5,1,5)
    plot(srange, data5)
    xlabel("S")
    ylabel("Sel")
    ylim(-0.4, 0.25)
    title("C = 1.0")
    tight_layout()
    return lineplots
end

#As there is more C, the effects of positive s is weaker and the effects of negative s is stronger
#As there is more C, moves away from linear function.