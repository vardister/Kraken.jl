export pekeris_env
export one_layer_env
export one_layer_slope_env
export two_layer_slope_env
export munk_env

### Standard Pekeris
"""
    pekeris_env(; c0=1500.0, cb=1600.0, ρ0=1000.0, ρb=1500.0, depth=100.0)

Canonical Pekeris waveguide: one isovelocity water layer over a homogeneous fluid half-space.

Returns `(ssp, layers, sspHS)` in the KRAKEN `.env` record layout, ready for
[`UnderwaterEnv`](@ref). `ssp` columns are `[z, cp, cs, ρ, αp, αs]`.

This is the environment with a closed-form solution — see [`PekerisUnderwaterEnv`](@ref) and
`find_kr(::PekerisUnderwaterEnv, freq)` — which makes it the reference for checking the
finite-difference solver.

# Arguments
- `c0`, `ρ0`: water column sound speed (m/s) and density (kg/m³).
- `cb`, `ρb`: bottom half-space sound speed and density. `cb` must exceed `c0` for trapped modes.
- `depth`: water depth (m), i.e. the interface with the half-space.
"""
function pekeris_env(; c0::Real=1500.0, cb::Real=1600.0, ρ0::Real=1000.0, ρb::Real=1500.0, depth::Real=100.0)
    # Input validation
    # for (param_name, param_value) in zip([:c0, :cb, :ρ0, :ρb, :depth], [c0, cb, ρ0, ρb, depth])
    #     if !isfinite(param_value)
    #         throw(DomainError(param_name, "Parameter must be a finite real number"))
    #     end
    #     if param_value <= 0
    #         throw(DomainError(param_name, "Parameter must be positive"))
    #     end
    # end

    # Water column
    α0 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    freq = 100.0
    z0 = depth

    # Columns are [z, cp, cs, ρ, αp, αs] — the KRAKEN .env SSP record layout. Built through a
    # promoted element type so ForwardDiff.Dual parameters survive into the environment.
    T = promote_type(typeof(c0), typeof(cb), typeof(ρ0), typeof(ρb), typeof(depth))
    ssp = zeros(T, 2, 6)
    ssp[1, :] = [0.0 c0 0.0 ρ0 α0 0.0]
    ssp[2, :] = [depth c0 0.0 ρ0 α0 0.0]

    layers = [0.0 0.0 depth]

    T2 = promote_type(typeof(depth), typeof(cb), typeof(ρb), typeof(αb))
    sspHS = zeros(T2, 2, 6)
    sspHS[1, :] = [0.0 343.0 0.0 0.00121 0.0 0.0]
    sspHS[2, :] = [depth cb 0.0 ρb αb 0.0]

    # env_dict = Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
    return ssp, layers, sspHS
end

### Standard 1-layer sediment model with constant sound speeds
"""
    one_layer_env(; c0=1500.0, c1=1550.0, cb=1600.0, ρ0=1000.0, ρ1=1500.0, ρb=2000.0, h0=100.0, h1=20.0)

Water column over one isovelocity sediment layer over a fluid half-space.

Returns `(ssp, layers, sspHS)`. `h0` is the water depth, `h1` the sediment thickness; `c1`/`ρ1` are
the sediment properties. The multi-medium case that exercises the interface conditions in
[`AcousticProblemCache`](@ref).
"""
function one_layer_env(; c0=1500.0, c1=1550.0, cb=1600.0, ρ0=1000.0, ρ1=1500.0, ρb=2000.0, h0=100.0, h1=20.0)
    # Water column
    α0 = 0.0
    # sediment layer
    α1 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    freq = 100.0
    z0 = h0
    z1 = h0 + h1

    # Columns are [z, cp, cs, ρ, αp, αs] — the KRAKEN .env SSP record layout.
    ssp = [
        0.0 c0 0.0 ρ0 α0 0.0
        z0 c0 0.0 ρ0 α0 0.0
        z0+eps(z0) c1 0.0 ρ1 α1 0.0
        z1 c1 0.0 ρ1 α1 0.0
    ]

    layers = [
        0.0 0.0 z0
        0.0 0.0 z1
    ]

    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        z1 cb 0.0 ρb αb 0.0
    ]

    return ssp, layers, sspHS
end

### Standard 1-layey model with slope in sound speed
"""
    one_layer_slope_env(; c0=1500.0, c1_1=1550.0, c1_2=1580.0, cb=1600.0, ρ0=1000.0, ρ1=1500.0, ρb=2000.0, h0=100.0, h1=20.0)

Like [`one_layer_env`](@ref), but the sediment sound speed ramps linearly from `c1_1` at its top to
`c1_2` at its base.

Returns `(ssp, layers, sspHS)`.
"""
function one_layer_slope_env(;
    c0=1500.0, c1_1=1550.0, c1_2=1580.0, cb=1600.0, ρ0=1000.0, ρ1=1500.0, ρb=2000.0, h0=100.0, h1=20.0
)
    # Water column
    α0 = 0.0
    # sediment layer
    α1 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    freq = 100.0
    z0 = h0
    z1 = h0 + h1

    # Columns are [z, cp, cs, ρ, αp, αs] — the KRAKEN .env SSP record layout.
    ssp = [
        0.0 c0 0.0 ρ0 α0 0.0
        z0 c0 0.0 ρ0 α0 0.0
        z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
        z1 c1_2 0.0 ρ1 α1 0.0
    ]

    layers = [
        0.0 0.0 z0
        0.0 0.0 z1
    ]

    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        z1 cb 0.0 ρb αb 0.0
    ]

    return ssp, layers, sspHS
end

### Standard 2-layer model with slope in sound speed
"""
    two_layer_slope_env(; c0=1500.0, c1_1=1550.0, c1_2=1580.0, c2_1=1600.0, c2_2=1650.0, cb=1800.0, ρ0=1000.0, ρ1=1500.0, ρ2=1600.0, ρb=2000.0, h0=100.0, h1=20.0, h2=20.0)

Water column over two sediment layers, each with a linear sound-speed gradient, over a fluid
half-space.

Returns `(ssp, layers, sspHS)`. Layer `i` ramps from `ci_1` to `ci_2` over thickness `hi`.
"""
function two_layer_slope_env(;
    c0=1500.0,
    c1_1=1550.0,
    c1_2=1580.0,
    c2_1=1600.0,
    c2_2=1650.0,
    cb=1800.0,
    ρ0=1000.0,
    ρ1=1500.0,
    ρ2=1600.0,
    ρb=2000.0,
    h0=100.0,
    h1=20.0,
    h2=20.0,
)
    # Water column
    α0 = 0.0
    # sediment layer
    α1 = 0.0
    # sediment layer
    α2 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    freq = 100.0
    z0 = h0
    z1 = h0 + h1
    z2 = h0 + h1 + h2

    # Columns are [z, cp, cs, ρ, αp, αs] — the KRAKEN .env SSP record layout.
    ssp = [
        0.0 c0 0.0 ρ0 α0 0.0
        z0 c0 0.0 ρ0 α0 0.0
        z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
        z1 c1_2 0.0 ρ1 α1 0.0
        z1+eps(z1) c2_1 0.0 ρ2 α2 0.0
        z2 c2_2 0.0 ρ2 α2 0.0
    ]

    layers = [
        0.0 0.0 z0
        0.0 0.0 z1
        0.0 0.0 z2
    ]

    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        z2 cb 0.0 ρb αb 0.0
    ]

    return ssp, layers, sspHS
end

"""
    munk_env()

Munk canonical deep-water sound-speed profile, 5000 m of water over a fluid half-space.

Returns `(ssp, layers, sspHS)`, sampled every 100 m. The profile is

```math
c(z) = 1500 \\, [1 + ε (ẑ - 1 + e^{-ẑ})], \\quad ẑ = 2(z - 1300)/1300, \\quad ε = 0.00737
```

so its minimum is exactly 1500.0 m/s at the 1300 m sound channel axis, where `ẑ = 0`. That value and
location are the defining property of the profile — the environment tests assert both.
"""
function munk_env()
    function c(z)
        ϵ = 0.00737
        zhat = 2 * (z - 1300.0) / 1300
        return 1500.0 * (1.0 + ϵ * (zhat - 1.0 + exp(-zhat)))
    end

    # Water column
    ρ0 = 1000.0
    α0 = 0.0
    # bottom half-space
    αb = 0.0
    ρb = 1500.0
    cb = 1600.0

    # other
    freq = 100.0
    z0 = 0.0
    zvec = 0:100:5000
    cvec = c.(zvec)

    # `permutedims`, not `transpose`: transpose returns a lazy Transpose wrapper, and every other
    # standard environment returns a plain Matrix. UnderwaterEnvFORTRAN promotes its three inputs
    # to a common type, which fails outright on Transpose-vs-Matrix.
    ssp = permutedims(hcat([[x[1], x[2], 0.0, ρ0, α0, 0.0] for x in zip(zvec, cvec)]...))
    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        zvec[end] cb 0.0 ρb αb 0.0
    ]
    layers = [0.0 0.0 zvec[end]]

    return ssp, layers, sspHS
end
