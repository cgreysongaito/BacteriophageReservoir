include("packages.jl")
include("bacteriophage_commoncode.jl")

## Symbolic analysis
@vars C

@vars r s b γ

# Without lysogeny
f(C) = r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C )

SymPy.solve(f(C), C)

SymPy.simplify(diff(f(C),C))

SymPy.simplify((( 0 * s^2 * (0 - 1) + r * (1 - 2 * 0) * (0 * s + 1)^2 + s * (1 - 2 * 0) * (0 * s + 1) ) ) / (0 * s + 1)^2)
# 0 equilibrium stable when s <-r
SymPy.simplify((( 1 * s^2 * (1 - 1) + r * (1 - 2 * 1) * (1 * s + 1)^2 + s * (1 - 2 * 1) * (1 * s + 1) ) ) / (1 * s + 1)^2)
# 1 equilibrium stable when s> -r / (r+1)
SymPy.simplify((( (-(r + s)/(r*s)) * s^2 * ((-(r + s)/(r*s)) - 1) + r * (1 - 2 * (-(r + s)/(r*s))) * ((-(r + s)/(r*s)) * s + 1)^2 + s * (1 - 2 * (-(r + s)/(r*s))) * ((-(r + s)/(r*s)) * s + 1) ) ) / ((-(r + s)/(r*s)) * s + 1)^2)
SymPy.solve(SymPy.simplify((( (-(r + s)/(r*s)) * s^2 * ((-(r + s)/(r*s)) - 1) + r * (1 - 2 * (-(r + s)/(r*s))) * ((-(r + s)/(r*s)) * s + 1)^2 + s * (1 - 2 * (-(r + s)/(r*s))) * ((-(r + s)/(r*s)) * s + 1) ) ) / ((-(r + s)/(r*s)) * s + 1)^2),s)


(( C * s^2 * (C - 1) + r * (1 - 2 * C) * (C * s + 1)^2 + s * (1 - 2 * C) * (C * s + 1) ) ) / (C * s + 1)^2

function modelwol(C,p)
    @unpack s, r = p
    return r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C )
end

let
    Crange = 0.0:0.0001:1.0
    data1 = [modelwol(C, BacPhagePar(s = -0.12)) for C in Crange]
    data2 = [modelwol(C, BacPhagePar(s = -0.0991)) for C in Crange]
    data3 = [modelwol(C, BacPhagePar(s = -0.08)) for C in Crange]
    testfull = figure(figsize=(8,2.5))
    subplot(1,3,1)
    plot(Crange, data1)
    ylabel("dC/dt")
    hlines(0.0, 0.0, 1.0, colors= "black")
    subplot(1,3,2)
    plot(Crange, data2)
    hlines(0.0, 0.0, 1.0, colors= "black")
    xlabel("C")
    subplot(1,3,3)
    plot(Crange, data3)
    hlines(0.0, 0.0, 1.0, colors= "black")
    tight_layout()
    # return testfull
    savefig(joinpath(abpath(), "figs/withoutlysogeny.png"))
end

function bifurcwol(s, p)
    @unpack r = p
    return ( -r - s ) / (r*s)
end

let
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data = [bifurcwol(s, BacPhagePar()) for s in srange]
    bifurcplot = figure(figsize=(5,4))
    plot(srange, data)
    ylabel("Ĉ")
    xlabel("s")
    ylim(-20, 20) #change coords
    hlines(0.0, st, en, colors= "black")
    hlines(1.0, st, en, colors= "black")
    # return bifurcplot
    savefig(joinpath(abpath(), "figs/bifurcwol.png"))
end

# With lysogeny - assuming r is tiny (0)
h(C) = ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + b * ( 1 - C )

SymPy.solve(h(C), C)

SymPy.simplify(diff(h(C),C))

SymPy.simplify(( 1 * s^2 * (1 - 1)  - b * (1 * s + 1)^2  + s * (1 - 2 * 1) * (1 * s + 1)) / (1 * s + 1)^2 )


SymPy.simplify(( (-b/(s*(b + 1))) * s^2 * ((-b/(s*(b + 1))) - 1)  - b * ((-b/(s*(b + 1))) * s + 1)^2  + s * (1 - 2 * (-b/(s*(b + 1))))*((-b/(s*(b + 1))) * s + 1)) / ((-b/(s*(b + 1)))*s + 1)^2)

# Stability when s is positive and when negative - identify bifurcation as change s from positive to negative - intuitively there should be a bifurcation because selection will push out gene (if strong enough)
#C*=1 stable when s > -b/(b+1)
#C*= -b/(s*(b + 1)) stable when s < -b/(b+1)


#### CHECK what type of function/graph is similar
# take limit of s to zero to see how figure changes shape
function modelwlrtiny(C, p)
    @unpack s, b = p
    return ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + b * ( 1 - C )
end

let
    Crange = 0.0:0.01:1.0
    data1 = [modelwlrtiny(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlrtiny(C, BacPhagePar(s = -0.02)) for C in Crange]
    data3 = [modelwlrtiny(C, BacPhagePar(s = -0.0099)) for C in Crange]
    data4 = [modelwlrtiny(C, BacPhagePar(s = 0.001)) for C in Crange]
    testfull = figure(figsize=(8,8))
    subplot(2,2,1)
    plot(Crange, data1)
    hlines(0.0, 0.0, 1.0, colors="black")
    ylabel("dC/dt")
    xlabel("C")
    subplot(2,2,2)
    plot(Crange, data2)
    xlabel("C")
    hlines(0.0, 0.0, 1.0, colors="black")
    subplot(2,2,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    hlines(0.0, 0.0, 1.0, colors="black")
    subplot(2,2,4)
    plot(Crange, data4)
    xlabel("C")
    hlines(0.0, 0.0, 1.0, colors="black")
    tight_layout()
    # return testfull
    savefig(joinpath(abpath(), "figs/withlysogenyrzero.png"))
end

let
    Crange = -5.0:0.01:5.0
    data1 = [modelwlrtiny(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlrtiny(C, BacPhagePar(s = -0.01)) for C in Crange]
    data3 = [modelwlrtiny(C, BacPhagePar(s = 1)) for C in Crange]
    testfull = figure(figsize=(8,3))
    subplot(1,3,1)
    plot(Crange, data1)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,2)
    plot(Crange, data2)
    xlabel("C")
    subplot(1,3,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    tight_layout()
    # return testfull
    savefig(joinpath(abpath(), "figs/withlysogenyrzero_generalgraph.png"))
end

let
    Crange = 0.0:0.001:1.0
    data1 = [modelwlrtiny(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlrtiny(C, BacPhagePar(s = -1)) for C in Crange]
    data3 = [modelwlrtiny(C, BacPhagePar(s = -0.9)) for C in Crange]
    testfull = figure(figsize=(8,3))
    subplot(1,3,1)
    plot(Crange, data1)
    ylabel("dC/dt")
    xlabel("C")
    # ylim(-1, 0)
    hlines(0.0, 0.0, 1.0, colors="black")
    subplot(1,3,2)
    plot(Crange, data2)
    xlabel("C")
    hlines(0.0, 0.0, 1.0, colors="black")
    subplot(1,3,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    hlines(0.0, 0.0, 1.0, colors="black")
    tight_layout()
    return testfull
    # savefig(joinpath(abpath(), "figs/withlysogenyrzero_generalgraph.png"))
end

function bifurcwlrtiny(s, p)
    @unpack b = p
    return -b / (s * ( 1 + b ) )
end

let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data = [bifurcwlrtiny(s, BacPhagePar()) for s in srange]
    bifurcplot = figure(figsize=(5,5))
    plot(srange, data)
    ylabel("Ĉ")
    xlabel("s")
    ylim(-1.1, 1.1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    # return bifurcplot
    savefig(joinpath(abpath(), "figs/bifurcwlrtiny.png"))
end

#If C* can go to infinity but we are bounded (proportion), how do we deal with bounds?

function simpeq(x)
    return 1 / (1 +x)
end

function medeq(x)
    return x / (1 + x)
end

function compeq(x)
    return  - x^2 / ( 1 + x )
end

let 
    xrange = -5.0:0.01:5.0
    data = [compeq(x) for x in xrange]
    test = figure()
    plot(xrange, data)
    return test
end

# With lysogeny
g(C) = r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + b * ( 1 - C )

SymPy.solve(g(C), C)

SymPy.simplify(diff(g(C),C))

SymPy.simplify(- ((1^2 * s^2) / (1 * s + 1)^2 ) - 2 * 1 * r - ((2 * 1 * s ) / (1 * s + 1)^2) - b + r + (s / (1 * s + 1)^2 ))

SymPy.solve(SymPy.simplify(- ((1^2 * s^2) / (1 * s + 1)^2 ) - 2 * 1 * r - ((2 * 1 * s ) / (1 * s + 1)^2) - b + r + (s / (1 * s + 1)^2 )), s)
# equilibrium at 1 is stable when s >-(b+r)/(b+r+1)




- ((C^2 * s^2) / (C * s + 1)^2 ) - 2 * C * r - ((2 * C * s ) / (C * s + 1)^2) - b + r + (s / (C * s + 1)^2 )
(-b * s - r - s)^2

SymPy.solve(s+2,s)
SymPy.solve(-(b*s + r + s - sqrt(b^2*s^2 - 2*b*r*s + 2*b*s^2 + r^2 + 2*r*s + s^2))/(2*r*s),s)
SymPy.solve(-(b*s + r + s + sqrt(b^2*s^2 - 2*b*r*s + 2*b*s^2 + r^2 + 2*r*s + s^2))/(2*r*s),s)

function modelwlr(C, p)
    @unpack s, r, b = p
    return r * C * ( 1 - C ) + ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + b * ( 1 - C )
end

let
    Crange = 0.0:0.01:1.0
    data1 = [modelwlr(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlr(C, BacPhagePar(s = -0.3)) for C in Crange]
    data3 = [modelwlr(C, BacPhagePar(s = -0.01)) for C in Crange]
    testfull = figure(figsize=(8,3))
    subplot(1,3,1)
    plot(Crange, data1)
    hlines(0.0, 0.0, 1.0, colors="black")
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,2)
    plot(Crange, data2)
    ylabel("dC/dt")
    xlabel("C")
    hlines(0.0, 0.0, 1.0, colors="black")
    subplot(1,3,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    hlines(0.0, 0.0, 1.0, colors="black")
    tight_layout()
    # return testfull
    savefig(joinpath(abpath(), "figs/withlysogeny.png"))
end

let
    Crange = -5.0:0.01:5.0
    data1 = [modelwlr(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlr(C, BacPhagePar(s = -0.01)) for C in Crange]
    data3 = [modelwlr(C, BacPhagePar(s = 0.5)) for C in Crange]
    testfull = figure(figsize=(8,3))
    subplot(1,3,1)
    plot(Crange, data1)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,2)
    plot(Crange, data2)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    tight_layout()
    # return testfull
    savefig(joinpath(abpath(), "figs/withlysogeny_generalgraph.png"))
end

let
    Crange = -5.0:0.01:5.0
    data1 = [modelwlr(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlr(C, BacPhagePar(s = -1)) for C in Crange]
    data3 = [modelwlr(C, BacPhagePar(s = -0.9)) for C in Crange]
    testfull = figure(figsize=(8,3))
    subplot(1,3,1)
    plot(Crange, data1)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,2)
    plot(Crange, data2)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    tight_layout()
    return testfull
    # savefig(joinpath(abpath(), "figs/withlysogeny_generalgraph.png"))
end

function equilwlr(s, p)
    @unpack b, r = p
    eq1 = (-(b*s + r + s) + sqrt(b^2*s^2 - 2*b*r*s + 2*b*s^2 + r^2 + 2*r*s + s^2))/(2*r*s)
    eq2 = (-(b*s + r + s) - sqrt(b^2*s^2 - 2*b*r*s + 2*b*s^2 + r^2 + 2*r*s + s^2))/(2*r*s)
    return [eq1, eq2]
end

function bioequilwlr(s, p)
    @unpack b, r = p
    eq2 = (-(b*s + r + s) - sqrt(b^2*s^2 - 2*b*r*s + 2*b*s^2 + r^2 + 2*r*s + s^2))/(2*r*s)
    return eq2
end

data = [equilwlr(s, BacPhagePar()) for s in -1:0.0001:1]
eq1 = zeros(length(-1:0.0001:1))
eq2 = zeros(length(-1:0.0001:1))
for i in 1:length(-1:0.0001:1)
    eq1[i] = data[i][1]
    eq2[i] = data[i][2]
end

let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data = [equilwlr(s, BacPhagePar()) for s in srange]
    eq1 = zeros(length(srange))
    eq2 = zeros(length(srange))
    for i in 1:length(srange)
        eq1[i] = data[i][1]
        eq2[i] = data[i][2]
    end
    bifurcplot = figure(figsize=(5,4))
    plot(srange, eq1, color = "green")
    plot(srange, eq2, color = "blue")
    ylabel("Ĉ")
    xlabel("s")
    ylim(-30, 30)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    # return bifurcplot
    savefig(joinpath(abpath(), "figs/bifurcwlr.png"))
end

let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data1 = [bioequilwlr(s, BacPhagePar()) for s in srange]
    data2 = [bioequilwlr(s, BacPhagePar(b=0.02)) for s in srange]
    data3 = [bioequilwlr(s, BacPhagePar(b=0.09)) for s in srange]
    bifurcplot = figure(figsize=(5,4))
    plot(srange, data1, color = "blue")
    plot(srange, data2, color = "green")
    plot(srange, data3, color = "red")
    ylabel("Ĉ")
    xlabel("s")
    xlim(-1.00, 0.2)
    ylim(0, 1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    return bifurcplot
    # savefig(joinpath(abpath(), "figs/bifurcwlr_changingb.png"))
end


#trying bifurcation of s with linear asexual reproduction function
function linequil(s, p)
    @unpack r, b = p
    return (-b + s) / (r)
end

let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data1 = [linequil(s, BacPhagePar()) for s in srange]
    data2 = [linequil(s, BacPhagePar(b=0.02)) for s in srange]
    data3 = [linequil(s, BacPhagePar(b=0.09)) for s in srange]
    bifurcplot = figure(figsize=(5,4))
    plot(srange, data1, color = "blue")
    plot(srange, data2, color = "green")
    plot(srange, data3, color = "red")
    ylabel("Ĉ")
    xlabel("s")
    xlim(-1.00, 1.0)
    ylim(0, 1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    return bifurcplot
    # savefig(joinpath(abpath(), "figs/bifurcwlr_changingb.png"))
end


#Horizontal Gene Transfer

function hgt(C, p)
    @unpack r, b = p
    return 100 * (r * C * ( 1 - C) + b * (1 - C))
end

let
    Crange = 0.0:0.01:1.0
    data1 = [hgt(C, BacPhagePar()) for C in Crange]
    test = figure()
    plot(Crange, data1)
    ylabel("HGT")
    xlabel("C")
    # return test
    savefig(joinpath(abpath(), "figs/HGT.png"))
end

# Selection function exploration
function conjug(C, r)
    # @unpack r = p
    return r * C #* (1-C)
end

let 
    Crange = 0.0:0.01:1.0
    data = [conjug(C, -0.1) for C in Crange]
    test = figure()
    plot(Crange, data)
    return test
end

function sel(C, p)
    @unpack s = p
    return (( s * C ) / ( 1 + s * C)) * (1 - C)

end

function sel2(C, half)
    #@unpack s = p
    return (1  / ( half + C))# *(1-C)
end


t(C) = (C  / ( (1/s) + C)) *(1-C)
SymPy.solve(diff(t(C),C),C)

let 
    Crange = 0.0:0.01:1.0
    data = [sel2(C, 0.5) for C in Crange]
    test = figure()
    plot(Crange, data)
    return test
end

function sel3(C, s, half)
    #@unpack s = p
    return ((s * C) / ( half + C)) *(1-C)
end

let 
    Crange = 0.0:0.01:1.0
    data = [sel3(C, 1, 0.4) for C in Crange]
    data2 = [sel3(C, 1, 0.1) for C in Crange]
    data3 = [sel3(C, -0.2, 0.4) for C in Crange]
    data4 = [sel3(C, -1.0, 0.4) for C in Crange]
    test = figure()
    plot(Crange, data)
    plot(Crange, data2)
    plot(Crange, data3)
    plot(Crange, data4)
    return test
end

#Simple linear selection function

## With lysogeny without conjugation
function bifurcwlwor_lin(s, p)
    @unpack b, r = p
    return -b / s
end

let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data = [bifurcwlwor_lin(s, BacPhagePar()) for s in srange]
    bifurcplot = figure(figsize=(5,5))
    plot(srange, data)
    ylabel("Ĉ")
    xlabel("s")
    ylim(-1.1, 1.1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    return bifurcplot
    # savefig(joinpath(abpath(), "figs/bifurcwlrtiny.png"))
end

##With lysogeny and conjugation
function modelwlr_lin(C, p)
    @unpack s, r, b = p
    return r * C * ( 1 - C ) +  s * C * ( 1 - C ) + b * ( 1 - C )
end

let
    Crange = 0.0:0.01:1.0
    data1 = [modelwlr(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlr(C, BacPhagePar(s = -0.3)) for C in Crange]
    data3 = [modelwlr(C, BacPhagePar(s = -0.01)) for C in Crange]
    testfull = figure(figsize=(8,3))
    subplot(1,3,1)
    plot(Crange, data1)
    hlines(0.0, 0.0, 1.0, colors="black")
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,2)
    plot(Crange, data2)
    ylabel("dC/dt")
    xlabel("C")
    hlines(0.0, 0.0, 1.0, colors="black")
    subplot(1,3,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    hlines(0.0, 0.0, 1.0, colors="black")
    tight_layout()
    # return testfull
    savefig(joinpath(abpath(), "figs/withlysogeny.png"))
end

let
    Crange = -5.0:0.01:5.0
    data1 = [modelwlr(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlr(C, BacPhagePar(s = -0.01)) for C in Crange]
    data3 = [modelwlr(C, BacPhagePar(s = 0.5)) for C in Crange]
    testfull = figure(figsize=(8,3))
    subplot(1,3,1)
    plot(Crange, data1)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,2)
    plot(Crange, data2)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    tight_layout()
    # return testfull
    savefig(joinpath(abpath(), "figs/withlysogeny_generalgraph.png"))
end

let
    Crange = -5.0:0.01:5.0
    data1 = [modelwlr(C, BacPhagePar(s = -1.1)) for C in Crange]
    data2 = [modelwlr(C, BacPhagePar(s = -1)) for C in Crange]
    data3 = [modelwlr(C, BacPhagePar(s = -0.9)) for C in Crange]
    testfull = figure(figsize=(8,3))
    subplot(1,3,1)
    plot(Crange, data1)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,2)
    plot(Crange, data2)
    ylabel("dC/dt")
    xlabel("C")
    subplot(1,3,3)
    plot(Crange, data3)
    xlabel("C")
    ylabel("dC/dt")
    tight_layout()
    return testfull
    # savefig(joinpath(abpath(), "figs/withlysogeny_generalgraph.png"))
end

function bifurcwlr_lin(s, p)
    @unpack b, r = p
    return -b / (r + s)
end

let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data = [bifurcwlr_lin(s, BacPhagePar()) for s in srange]
    bifurcplot = figure(figsize=(5,5))
    plot(srange, data)
    ylabel("Ĉ")
    xlabel("s")
    ylim(-1.1, 1.1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    return bifurcplot
    # savefig(joinpath(abpath(), "figs/bifurcwlrtiny.png"))
end

# With lytic taking carrier bacteria out
h(C) = r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C ) - γ * C + b * ( 1 - C )

SymPy.simplify(SymPy.solve(h(C), C))

@with_kw mutable struct BacPhagePar_remove
    r = 0.1
    s = 0.1
    b = 0.01
    γ = 0.01
end

function equilwlr_remove(s, p)
    @unpack b, r, γ = p
    eq1 = -(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))/(3*(-27*b/(2*r*s) + sqrt(-4*(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))^3 + (-27*b/(r*s) - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(r^2*s^2) + 2*(b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^2)/2 - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(2*r^2*s^2) + (b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^(1/3)) - (-27*b/(2*r*s) + sqrt(-4*(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))^3 + (-27*b/(r*s) - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(r^2*s^2) + 2*(b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^2)/2 - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(2*r^2*s^2) + (b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^(1/3)/3 - (b*s - r*s + r + s*γ + s)/(3*r*s)
    # eq2 = -(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))/(3*(-1/2 - sqrt(3)*I/2)*(-27*b/(2*r*s) + sqrt(-4*(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))^3 + (-27*b/(r*s) - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(r^2*s^2) + 2*(b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^2)/2 - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(2*r^2*s^2) + (b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^(1/3)) - (-1/2 - sqrt(3)*I/2)*(-27*b/(2*r*s) + sqrt(-4*(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))^3 + (-27*b/(r*s) - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(r^2*s^2) + 2*(b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^2)/2 - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(2*r^2*s^2) + (b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^(1/3)/3 - (b*s - r*s + r + s*γ + s)/(3*r*s)
    # eq3 = -(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))/(3*(-1/2 + sqrt(3)*I/2)*(-27*b/(2*r*s) + sqrt(-4*(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))^3 + (-27*b/(r*s) - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(r^2*s^2) + 2*(b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^2)/2 - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(2*r^2*s^2) + (b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^(1/3)) - (-1/2 + sqrt(3)*I/2)*(-27*b/(2*r*s) + sqrt(-4*(-3*(-b*s + b - r - s + γ)/(r*s) + (b*s - r*s + r + s*γ + s)^2/(r^2*s^2))^3 + (-27*b/(r*s) - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(r^2*s^2) + 2*(b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^2)/2 - 9*(-b*s + b - r - s + γ)*(b*s - r*s + r + s*γ + s)/(2*r^2*s^2) + (b*s - r*s + r + s*γ + s)^3/(r^3*s^3))^(1/3)/3 - (b*s - r*s + r + s*γ + s)/(3*r*s)
    return [eq1]#, eq2, eq3]
end


let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data = [equilwlr_remove(s, BacPhagePar_remove()) for s in srange]
    eq1 = eq2 = eq3 = zeros(length(srange))
    for i in 1:length(srange)
        eq1[i] = data[i][1]
        # eq2[i] = data[i][2]
        # eq3[i] = data[i][3]
    end
    bifurcplot = figure(figsize=(5,4))
    plot(srange, eq1, color = "green")
    # plot(srange, eq2, color = "blue")
    # plot(srange, eq3, color = "red")
    ylabel("Ĉ")
    xlabel("s")
    ylim(-30, 30)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, linestyles="dashed", colors= "black")
    return bifurcplot
    # savefig(joinpath(abpath(), "figs/bifurcwlr.png"))
end


#proove conjugation more than bacteriophage
@vars x

g(x) = -2*(x^3) + 4 * (x^2) -2*x

SymPy.solve(g(x), x)

eq2= 0.1:0.1:1.0

test2 = rand(Float64, (2, 2))
test2[1,2]
test2[2,1]

function conj_proof(eq1,eq2)
    return (2 * eq2 - eq2^2 - 2 * eq1 + eq1^2) / (3 * eq2^2 - 2 * eq2^3 - 3 * eq1^2 + 2 * eq1^3)
end

conj_proof(0.5,0.6)
conj_proof(0.6,0.5)

#seems to be symmetrical! what does this mean?

#rows are eq2
#columns are eq1
let
    eq1 = eq2 = 0.1:0.1:1.0
    test = zeros(10,10)
    for i in 1:10
        for j in 1:10
            test[i,j] = (2 * eq2[i] - eq2[i]^2 - 2 * eq1[j] + eq1[j]^2) / (3 * eq2[i]^2 - 2 * eq2[i]^3 - 3 * eq1[j]^2 + 2 * eq1[j]^3)
        end
    end
    return test
end


function dcdt(C,par)
    @unpack r, b, s = par
    return r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + b * ( 1 - C )
end

dcdt(0.4499999, BacPhagePar(s = 0.01))

function dHGTdt(C,dcdtval, par)
    @unpack r, b = par
    return 100 * (r * dcdtval - 2 * r * C * dcdtval - b * dcdtval)
end

dHGTdt(0.449999, dcdt(0.4499999, BacPhagePar(s = 0.01)), BacPhagePar(s = 0.01))

# Time series analysis


# "Potential" analysis

#

## Ideas
# White noise to red noise (selection)
# sine wave for selection
#investigate when r and s and l are really small
# investigate when selection places c equilibrium close to maximum of HGT and away from maximum
#slow fast analysis (with selection being much slower than growth rates of bacteria)