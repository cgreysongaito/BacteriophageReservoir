include("packages.jl")

##


function abpath()
    replace(@__DIR__, "scripts" => "")
end #function to create a absolute path as a character string

@with_kw mutable struct BacPhagePar
    r = 0.01
    s = 0.1
    l = 0.01
end

function bacphage!(du, u, p, t,)
    @unpack r, l = p
    C = u
    du[1] = r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + l * ( 1 - C )
    return
end

function bacphage(u, par)
    du = similar(u)
    bacphage!(du, u, p, 0.0)
    return du
end

## Symbolic analysis
@vars C

@vars r s l

# Without lysogeny
f(C) = r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C )

SymPy.solve(f(C), C)

# With lysogeny
g(C) = r * C * (1 - C) + ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + l * ( 1 - C )

SymPy.solve(g(C), C)

# With lysogeny - assuming r is tiny (0)
h(C) = ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + l * ( 1 - C )

SymPy.solve(h(C), C)

# Stability when s is positive and when negative - identify bifurcation as change s from positive to negative - intuitively there should be a bifurcation because selection will push out gene (if strong enough)


## Geometric analysis
# With lysogeny - assuming r is tiny (0)
#### CHECK what type of function/graph is similar
# take limit of s to zero to see how figure changes shape
function selec(C, p)
    @unpack s, l = p
    return ( s / ( 1 + s * C ) ) * C * ( 1 - C )
end 

function lyso(C, p)
    @unpack l = p
    return l * ( 1 - C )
end

function fullmodel(C, p)
    @unpack s, l = p
    return ( s / ( 1 + s * C ) ) * C * ( 1 - C ) + l * ( 1 - C )
end

let
    Crange = 0.0:0.01:1.0
    datasel = [selec(C, BacPhagePar(s = 0.018)) for C in Crange]
    datal = [lyso(C, BacPhagePar(s = 0.00018)) for C in Crange]
    test = figure()
    plot(Crange, datasel)
    plot(Crange, datal)
    return test
end

let
    Crange = -100.0:0.01:5.0
    datafull = [fullmodel(C, BacPhagePar(s = 0.001)) for C in Crange]
    testfull = figure()
    plot(Crange, datafull)
    hlines(0.0, -5.0, 5.0)
    return testfull
end

function equil(s, p)
    @unpack l = p
    return -l / (s * ( 1 + l ) )
end

let 
    srange = -5.0:0.01:5.0
    datas = [equil(s, BacPhagePar()) for s in srange]
    testeq = figure()
    plot(srange, datas)
    ylim(-1.0,1.0)
    return testeq
end


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
# Time series analysis


# "Potential" analysis



## Ideas
# White noise to red noise (selection)
# sine wave for selection
#investigate when r and s and l are really small