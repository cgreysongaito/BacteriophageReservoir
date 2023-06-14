include("packages.jl")
include("bacteriophage_commoncode.jl")

## Symbolic analysis

# With conjugation without lysogeny
@vars r s
@vars C
f(C) = r * C * ( 1 - C ) + s * C * ( 1 - C )

#Equilibria
SymPy.solve(f(C), C)

#Eigenvalue
SymPy.simplify(diff(f(C),C))


## With lysogeny without conjugation
@vars s b
@vars C
g(C) = s * C * ( 1 - C ) + b * ( 1 - C )

#Equilibria
SymPy.solve(g(C), C)

#Eigenvalue
SymPy.simplify(diff(g(C),C))

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
h(C) = r * C * (1 - C) +  s * C * ( 1 - C ) + b * ( 1 - C )

#Equilibria
SymPy.solve(h(C), C)

#Eigenvalue
SymPy.simplify(diff(h(C),C))
SymPy.simplify(-2*r - 2 *s - b +r +s) #C=1
SymPy.solve(SymPy.simplify(-2*r - 2 *s - b +r +s), s) #bifurcation point