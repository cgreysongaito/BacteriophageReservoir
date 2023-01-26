include("packages.jl")
include("bacteriophage_commoncode.jl")

## Symbolic analysis

## With lysogeny without conjugation
@vars s b
@vars C
f(C) = s * C * ( 1 - C ) + b * ( 1 - C )

#Equilibria
SymPy.solve(f(C), C)

#Eigenvalue
SymPy.simplify(diff(f(C),C))

function bifurcwlwor_lin(s, p)
    @unpack b, r = p
    return -b / s
end

let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data = [bifurcwlwor_lin(s, BacPhagePar()) for s in srange]
    bifurcwlworplot = figure(figsize=(5,5))
    plot(srange, data)
    ylabel("Ĉ")
    xlabel("s")
    ylim(-1.1, 1.1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    return bifurcwlworplot
    # savefig(joinpath(abpath(), "figs/bifurcwlwor.png"))
end

##With lysogeny and conjugation
@vars r s b
@vars C
g(C) = r * C * (1 - C) +  s * C * ( 1 - C ) + b * ( 1 - C )

#Equilibria
SymPy.solve(g(C), C)

#Eigenvalue
SymPy.simplify(diff(g(C),C))
SymPy.simplify(-2*r - 2 *s - b +r +s) #C=1
SymPy.solve(SymPy.simplify(-2*r - 2 *s - b +r +s), s) #bifurcation point

function bifurcwlr(s, p)
    @unpack r, b = p
    return -b / (r+s)
end

let 
    st = -1.0
    en = 1.0 
    srange = st:0.0001:en
    data1 = [bifurcwlr(s, BacPhagePar()) for s in srange]
    data2 = [bifurcwlr(s, BacPhagePar(b=0.02)) for s in srange]
    data3 = [bifurcwlr(s, BacPhagePar(b=0.09)) for s in srange]
    bifurcwlrplot = figure(figsize=(5,4))
    plot(srange, data1, color = "blue")
    plot(srange, data2, color = "green")
    plot(srange, data3, color = "red")
    ylabel("Ĉ")
    xlabel("s")
    xlim(-1.00, 1.0)
    ylim(0, 1)
    hlines(0.0, st, en, linestyles="dashed", colors= "black")
    hlines(1.0, st, en, colors= "black")
    return bifurcwlrplot
    # savefig(joinpath(abpath(), "figs/bifurcwlr_changingb.png"))
end