include("packages.jl")

@with_kw mutable struct BacPhageDelayPar
    r = 0.1
    s = 0.1
    b = 0.01
    τ = 10
end

function bacphagedelay!(du, u, p, t,)
    @unpack r, b, s = p
    du[1] = r * u[1] * (1 - u[1]) + ( s / ( 1 + s * u[1] ) ) * u[1] * ( 1 - u[1] ) + b * ( 1 - u[1] )
    return
end