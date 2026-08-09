### Load packages
using LinearAlgebra
using Statistics
using Roots
using LinearSolve
using UnPack
using DataInterpolations
import NaNMath as nm
# `ignore_derivatives` marks the mesh-refinement loop's convergence test as the discrete decision it
# is. The rules themselves live in `kraken_ad.jl`; this is the one place the seam needs the package.
using ChainRulesCore: ignore_derivatives

### Docs
using DocStringExtensions

# Exports
export SampledSSP, SampledDensity
export soundspeed, maxsoundspeed, density
export UnderwaterEnv, AcousticProblemProperties, UnderwaterEnvFORTRAN
export AcousticProblemCache, bisection, solve_for_kr, inverse_iteration, det_sturm, kraken_jl, find_kr, get_g
export finite_difference_coefficients, mode_eigenvector, normalize_mode

### The differentiable seam
#
# Reverse-mode AD (Milestone 4) does not trace the solver. Rules are attached at `solve_for_kr` and
# `inverse_iteration`, which makes everything they call *opaque* to the tape. That splits this file
# in two, and the split has to be maintained deliberately — moving a computation across it silently
# either breaks reverse mode or makes it needlessly slow.
#
# ON THE SEAM — traced by the AD backend, so these must stay allocating and free of any mutation of
# arrays that outlive the expression that made them:
#   * `soundspeed`, `density`, `get_thickness`             — parameters entering from the environment
#   * `a_element`, `e_element`                             — elementwise coefficient formulas
#   * `finite_difference_coefficients`                     — the whole assembly `(env, props) -> (a, e, λ)`
#   * `AcousticProblemCache(env, props)`                   — a thin wrapper over the above
#   * `get_g`                                              — bottom half-space term, also used inside the rules
#   * `get_Nz_vec`, `get_z_vec`, `AcousticProblemProperties` — the mesh; see the caveat below
#   * `richard_extrap`, `h_extrap_matrix`, and the mesh-refinement loop in `kraken_jl` — which is
#     why that loop grows `h_list`/`krs_all` with `vcat` instead of `push!`ing into them, and why
#     its convergence test sits under `ignore_derivatives`
#   * `normalize_mode`, `integral_trapz`                   — the mode's energy normalization
#   * `find_kr`, `inverse_iteration`                       — the loops over the two rule-bearing
#                                                            functions below; both had to stop
#                                                            filling preallocated arrays
#
# BEHIND THE SEAM — never traced; free to mutate the cache, branch on values, and iterate:
#   * `det_sturm`, `scale_const`                           — differentiated analytically by the rule
#   * `bisection`                                          — integer mode counting, `@non_differentiable`
#   * `solve_for_kr`                                       — rrule: implicit function theorem on `det_sturm`
#   * `create_finite_diff_matrix!`, `return_finite_diff_matrix!` — in-place cache updates
#   * `mode_eigenvector`                                   — rrule: eigenvector adjoint
#
# The mesh is only *half* discrete, and getting that wrong loses a gradient silently. `Nz_vec` is an
# integer point count and is genuinely piecewise constant in the parameters. But `Δz_vec = h ./ Nz_vec`
# is linear in the layer thickness, and `zn_vec` is built from it — so both carry real derivative
# information with respect to layer depths, and the discretization-error part of `∂kr/∂depth` flows
# through them. `AcousticProblemProperties` is therefore *on* the seam and its construction must stay
# non-mutating; only `Nz_vec` is treated as a constant.

### Main Types
### Sound Speed Profile
abstract type SoundSpeedProfile end
abstract type SampledSSP <: SoundSpeedProfile end

"""
Sound speed profile based on measurements at discrete depths `z` in meters and sound speed `c` in m/s.
"""
struct SampledSSP1D{T1,T2,T3} <: SampledSSP
    z::Vector{T1}
    c::Vector{T2}
    f::T3
    function SampledSSP1D(depth, c, f)
        interp = f(c, depth; extrapolation=ExtrapolationType.Constant)
        return new{eltype(depth),eltype(c),typeof(interp)}(-depth, c, interp)
    end
end

"""
    SampledSSP(depth, c)
    SampledSSP(depth, c, type::Symbol)

Constructor for `SampledSSP1D`.

    Create a sound speed profile based on measurements at discrete depths `z` in meters and sound speed `c` in m/s.
    Two options for interpolation are available: `:linear` and `:smooth`.
"""
SampledSSP(depth, c) = SampledSSP1D(depth, c, DataInterpolations.LinearInterpolation)
SampledSSP(depth, c, type::Symbol) = SampledSSP1D(depth, c, type)

function Base.show(io::IO, ssp::SampledSSP1D{T1,T2,T3}) where {T1,T2,T3}
    return print(io, "SampledSSP1D{", T1, ",", T2, ",", ssp.f, "}(", length(ssp.z), " points)")
end

### Density Profile
abstract type DensityProfile end
abstract type SampledDensity <: DensityProfile end

"""
$(TYPEDEF)
Density profile based on measurements at discrete depths `z` in meters and density `ρ` in kg/m³.
"""
struct SampledDensity1D{T1,T2,T3} <: SampledDensity
    z::Vector{T1}
    ρ::Vector{T2}
    f::T3
    # Constructor for Type inputs
    function SampledDensity1D(depth, ρ, f)
        interp = f(ρ, depth; extrapolation=ExtrapolationType.Constant)
        return new{eltype(depth),eltype(ρ),typeof(interp)}(-depth, ρ, interp)
    end
end

"""
    SampledDensity(depth, ρ)

Constructor for `SampledDensity1D`.

Create a density profile based on measurements at discrete depths `z` in meters and density `ρ` in kg/m³.
Two options for interpolation are available: `:linear` and `:smooth`.
"""
SampledDensity(depth, ρ) = SampledDensity1D(depth, ρ, DataInterpolations.LinearInterpolation)
SampledDensity(depth, ρ, type::Symbol) = SampledDensity1D(depth, ρ, type)

function Base.show(io::IO, ρint::SampledDensity1D{T1,T2,T3}) where {T1,T2,T3}
    # Mirrors the SampledSSP1D method above: the third slot is the interpolant, `f`. There is no
    # `.type` field on this struct — printing one threw on every `show` of a density profile.
    return print(io, "SampledDensity1D{", T1, ",", T2, ",", ρint.f, "}(", length(ρint.z), " points)")
end

### Underwater Environment

"""
Underater environment used for running the FORTRAN version of KRAKEN.
"""
struct UnderwaterEnvFORTRAN{T<:Real}
    ssp::Matrix{T}
    layers::Matrix{T}
    sspHS::Matrix{T}
end

"""
    UnderwaterEnvFORTRAN(ssp, layers, sspHS)

Constructor for `UnderwaterEnvFORTRAN`.
The input matrices are the same as the matrices used in the FORTRAN version of KRAKEN.
These are the sound speed profile `ssp`, the layer information `layers`, and the sound speed profile
at the surface and bottom half-space `sspHS`.
"""
function UnderwaterEnvFORTRAN(ssp, layers, sspHS)
    ssp, layers, sspHS = promote(ssp, layers, sspHS)
    return UnderwaterEnvFORTRAN{eltype(ssp)}(ssp, layers, sspHS)
end

function Base.show(io::IO, ::UnderwaterEnvFORTRAN{T}) where {T}
    return print(io, "UnderwaterEnvFORTRAN{$T}")
end

"""
Underwater environment containing the sound speed profile and density profile.
"""
struct UnderwaterEnv{T1<:SoundSpeedProfile,T2<:DensityProfile,T3<:Real}
    c::T1
    ρ::T2
    cb::T3
    ρb::T3
    h_vec::Vector{T3}
    layer_depth::Vector{T3}
    depth::T3
end

"""
    UnderwaterEnv(ssp, layers, sspHS)

Constructor for `UnderwaterEnv`.

Create an underwater environment based on the sound speed profile `ssp`, the layer information `layers`, and the sound speed profile
at the surface and bottom half-space `sspHS`.

`depth` — the interface with the bottom half-space — is taken from `layers[end, 3]` in both this
constructor and the `UnderwaterEnvFORTRAN` one. `layers` is what defines the media boundaries (it is
already the sole input to `get_thickness`, which sizes the finite-difference mesh); the `ssp` table
is only samples *within* those media. They coincide in every well-formed environment, but when they
disagree it is `ssp` that is short or long, not `layers` that is wrong.
"""
function UnderwaterEnv(ssp, layers, sspHS)
    c = SampledSSP(ssp[:, 1], ssp[:, 2])
    ρ = SampledDensity(ssp[:, 1], ssp[:, 4])
    ρb = sspHS[2, 4]
    cb = sspHS[2, 2]
    layer_thickness = get_thickness(layers)
    layer_depth = layers[:, 3]
    depth = layers[end, 3]
    return UnderwaterEnv{typeof(c),typeof(ρ),typeof(cb)}(c, ρ, cb, ρb, layer_thickness, layer_depth, depth)
end

"""
    UnderwaterEnv(krak_ssp::UnderwaterEnvFORTRAN{T}) where {T}

Constructor for `UnderwaterEnv` using the `UnderwaterEnvFORTRAN` struct.
"""
function UnderwaterEnv(krak_ssp::UnderwaterEnvFORTRAN{T}) where {T}
    c = SampledSSP(krak_ssp.ssp[:, 1], krak_ssp.ssp[:, 2])
    ρ = SampledDensity(krak_ssp.ssp[:, 1], krak_ssp.ssp[:, 4])
    ρb = krak_ssp.sspHS[2, 4]
    cb = krak_ssp.sspHS[2, 2]
    layer_thickness = get_thickness(krak_ssp.layers)
    layer_depth = krak_ssp.layers[:, 3]
    depth = krak_ssp.layers[end, 3]  # see the note on the (ssp, layers, sspHS) constructor
    return UnderwaterEnv{typeof(c),typeof(ρ),T}(c, ρ, cb, ρb, layer_thickness, layer_depth, depth)
end

function Base.show(io::IO, ::UnderwaterEnv{T1,T2,T3}) where {T1,T2,T3}
    return print(io, "UnderwaterEnv{$T1, $T2, $T3}")
end

### Sound Speed and Density Functions to extract values from profiles at a give depth from profiles
"""
    soundspeed(ssp::SoundSpeedProfile, x, y, z)

Get sound speed at location (`x`, `y`, `z`). If a sound speed profile is range
independent, `x` and `y` may be ignored. `z` is generally negative, since the
sea surface is the datum and z-axis points upwards.
"""
function soundspeed end

soundspeed(ssp::SampledSSP, z) = ssp.f(z)

"""
    maxsoundspeed(ssp::SoundSpeedProfile)

Get the maximum sound speed from the sound speed profile.
"""
function maxsoundspeed end

maxsoundspeed(ssp::SampledSSP) = maximum(ssp.c)

"""
    density(ρ::DensityProfile, x, y, z)

Get density at location (`x`, `y`, `z`). If a density profile is range
independent, `x` and `y` may be ignored. `z` is generally negative, since the
sea surface is the datum and z-axis points upwards.
"""
function density end

density(ρ::SampledDensity1D, z) = ρ.f(z)

### Finite Difference Scheme
### Functions that convert SSP information (similar to KRAKEN) to environment and problem structs
function get_thickness(layers::Matrix{<:Real})
    # `vcat(0.0, col)` rather than `vcat(0.0, col...)`: splatting the column routes the whole thing
    # through `Core._apply_iterate`, whose reverse-mode adjoint cannot be composed with the
    # `getindex` pullback feeding it. Same result, and layer thicknesses stay differentiable — which
    # they must be, since every depth derivative enters the mesh through them.
    a = vcat(0.0, layers[:, 3])
    return a[2:end] - a[1:(end - 1)]
end

"""
    get_Nz_vec(env::UnderwaterEnv, freq; n_per_wavelength=20, factor=1)

Get the number of mesh points and the mesh spacing for each layer of the `env` for building the finite-difference scheme.
This process is dependent on the frequency `f`.

# Arguments
- `env::UnderwaterEnv`: Underwater environment struct.
- `f::Real`: Frequency of the acoustic problem.
- `n_per_wavelength::Int=20`: Number of mesh points per wavelength.
- `factor::Int=1`: factor of the mesh spacing.
"""
function get_Nz_vec(env::UnderwaterEnv, freq; n_per_wavelength=20, factor=1)
    ω = 2π * freq
    @assert ω >= 0 "Frequency must be non-negative"
    @assert maxsoundspeed(env.c) < env.cb
    kr_max = ω / env.cb  # here we assume the bottom half-space sound speed is highest
    Lmin = 2π / kr_max   # The lowest wavelength available in the problem
    # 20 points per wavelength. `factor` is for Richardson extrapolation.
    Δz = Lmin / n_per_wavelength

    # `Nz_vec` is an integer count and so is piecewise constant in the parameters, but `Δz_vec` is
    # `h / Nz` — *linear* in the layer thickness, hence differentiable. Written as comprehensions
    # rather than a fill-in loop so reverse-mode AD can trace it; see the seam note at the top.
    Nz_vec = [n_mesh_points(h, Δz, factor) for h in env.h_vec]
    Δz_vec = env.h_vec ./ Nz_vec
    return Nz_vec, Δz_vec
end

"""
    n_mesh_points(h, Δz, factor)

Number of mesh points for a layer of thickness `h` at target spacing `Δz`, never fewer than 10.

Split out as its own function only so the AD rules can mark it non-differentiable in one place:
`ceil` has no derivative and reverse mode refuses to trace it, even though the result is an `Int` and
is genuinely piecewise constant in the parameters.
"""
n_mesh_points(h, Δz, factor) = max(10, ceil(Int, h / Δz) * factor)

"""
    get_z_vec(env::UnderwaterEnv, Nz_vec, Δz_vec)

Get the depth vector for each layer of the underwater environment according to the number of mesh points `Nz_vec`
 and mesh spacing `Δz_vec`.
"""
function get_z_vec(env::UnderwaterEnv, Nz_vec, Δz_vec)
    # Layer `i` runs from the previous layer's bottom to its own, and its mesh starts one step in
    # (there is no `z = 0` sample — see the Architecture Decisions). Non-mutating so reverse-mode AD
    # can trace the layer-depth dependence.
    z_starts = vcat(zero(eltype(env.layer_depth)), env.layer_depth[1:(end - 1)])
    return [layer_mesh(z_starts[i] + Δz_vec[i], env.layer_depth[i], Nz_vec[i]) for i in eachindex(Nz_vec)]
end

"""
    layer_mesh(z_top, z_bot, Nz)

`Nz` equally spaced depths from `z_top` to `z_bot` inclusive — `collect(range(z_top, z_bot, Nz))`,
no more.

It exists as a named function purely so `src/kraken_ad.jl` can hang an `rrule` on it. `range` builds
a `StepRangeLen` whose fields are `Base.TwicePrecision`, which reverse-mode AD cannot construct, and
the obvious workaround — open-coding the interpolation as a broadcast — is **not** acceptable here:
it moves interior mesh points by up to ~13 ulp, and that is enough to push `two_layer_slope` past the
1e-4 wavenumber tolerance it is cross-validated against. `range` stays; the derivative (which is just
linear interpolation, exactly) is supplied by the rule instead.
"""
layer_mesh(z_top, z_bot, Nz) = collect(range(z_top, z_bot, Nz))

"""
    AcousticProblemProperties(env::UnderwaterEnv, freq; factor=1, n_per_wavelength=20)

Properties of the acoustic problem based on the underwater environment `env` and frequency `freq`.
`factor` is a factor for the mesh spacing and `n_per_wavelength` is the number of mesh points per wavelength.
"""
struct AcousticProblemProperties{T<:Real,T2<:Real}
    freq::T
    Nz_vec::Vector{Int}
    Δz_vec::Vector{T2}
    zn_vec::Vector{Vector{T2}}
end

"""
    AcousticProblemProperties(env::UnderwaterEnv, freq; factor::Int=1, n_per_wavelength=20)

Get the properties of the acoustic problem based on the underwater environment `env` and frequency `freq`.
"""
function AcousticProblemProperties(env::UnderwaterEnv, freq; factor::Int=1, n_per_wavelength=20)
    if freq isa Int
        freq = float(freq)
    end
    Nz_vec, Δz_vec = get_Nz_vec(env, freq; factor=factor, n_per_wavelength=n_per_wavelength)
    zn_vec = get_z_vec(env, Nz_vec, Δz_vec)

    return AcousticProblemProperties{eltype(freq),eltype(Δz_vec)}(freq, Nz_vec, Δz_vec, zn_vec)
end

function Base.show(io::IO, props::AcousticProblemProperties{T,T2}) where {T,T2}
    return print(io, "AcousticProblemProperties{", T, ",", T2, "}(", length(props.Nz_vec), " layers)")
end

### Prepare vectors for the finite difference scheme
a_element(c, ρ, f, h) = @. (-2 + h^2 * (2pi * f / c)^2) / (h * ρ)
e_element(ρ, h) = @. 1 / (h * ρ)

"""
    get_g(kr, env::UnderwaterEnv, props::AcousticProblemProperties)

Get the value of `g` for the bottom half-space finite-difference element.

Only defined for `kr >= 2π * freq / cb`, i.e. for modes that are evanescent in the bottom
half-space — below that cutoff the vertical wavenumber in the bottom is real (a radiating,
leaky mode) and this real-valued formulation has no solution. `bisection` therefore only ever
searches `kr ∈ [ω/cb, max(ω/c)]`; calling `get_g` (or `det_sturm`) below the cutoff throws a
`DomainError` from `sqrt`.
"""
function get_g(kr, env::UnderwaterEnv, props::AcousticProblemProperties)
    g = sqrt(kr^2 - (2pi * props.freq / env.cb)^2) / env.ρb
    return g
end

"""
Cache for the acoustic problem vectors.
"""
mutable struct AcousticProblemCache{T}
    a_vec::T
    e_vec::T
    λ_scaling::T
    A::Tridiagonal
end

"""
    finite_difference_coefficients(env::UnderwaterEnv, props::AcousticProblemProperties)

Assemble the finite-difference coefficients for `env` on the mesh described by `props`. Returns
`(a_vec, e_vec, λ_scaling)`, the three vectors [`AcousticProblemCache`](@ref) stores.

This is the whole dependence of the discretized problem on the environment parameters, isolated as a
pure function: it allocates its results and mutates nothing — not its arguments, not a cache, not
even a locally allocated buffer. That is what puts it *on the differentiable seam* (see the note at
the top of this file). Reverse-mode AD traces this function, so keep it free of `setindex!`,
`push!`, and in-place broadcasts; the mutating hot path lives behind the rules in
`create_finite_diff_matrix!` and never needs to be traced.

- `a_vec` — main diagonal, from `a_element`; at each interface between two media the value
  is the average of the coefficients on either side.
- `e_vec` — off-diagonals, from `e_element`.
- `λ_scaling` — the factor multiplying `kr²` in the Sturm sequence: a two-point moving average of
  `e_vec * Δz²`, with the final entry halved for the bottom half-space boundary condition.

The mesh (`Nz_vec`, `Δz_vec`, `zn_vec`) is piecewise constant in the parameters and is treated as a
constant here — only `env` carries derivative information.
"""
function finite_difference_coefficients(env::UnderwaterEnv, props::AcousticProblemProperties)
    #TODO: create vector that generalizes well to different types
    T = promote_type(eltype(env.c.c), eltype(env.ρ.ρ), typeof(env.cb), typeof(props.freq))

    # Flatten the per-layer mesh once, so the coefficients are a single broadcast over the whole
    # column rather than a loop that fills slices of a preallocated array.
    zn = reduce(vcat, props.zn_vec)
    Δzn = reduce(vcat, [fill(props.Δz_vec[i], props.Nz_vec[i]) for i in eachindex(props.Nz_vec)])
    cn = soundspeed(env.c, zn)
    ρn = density(env.ρ, zn)

    a_raw = T.(a_element(cn, ρn, props.freq, Δzn))
    e_vec = T.(e_element(ρn, Δzn))

    # Interface conditions between layers: the last point of each layer but the last carries the
    # average of the two media's coefficients. `interfaces` are the cumulative mesh counts.
    interfaces = cumsum(props.Nz_vec)[1:(end - 1)]
    a_vec = if isempty(interfaces)
        a_raw
    else
        map(k -> k in interfaces ? (a_raw[k] + a_raw[k + 1]) / 2 : a_raw[k], eachindex(a_raw))
    end

    # λ_scaling: pairwise mean of e * Δz² over the column, with the bottom entry halved.
    s = e_vec .* Δzn .^ 2
    λ_scaling = vcat((s[1:(end - 1)] .+ s[2:end]) ./ 2, e_vec[end] * props.Δz_vec[end]^2 / 2)

    return a_vec, e_vec, λ_scaling
end

"""
    AcousticProblemCache(env::UnderwaterEnv, props::AcousticProblemProperties)

Prepare the vectors `a_vec`, `e_vec`, and `λ_scaling` for the acoustic problem. Return an
`AcousticProblemCache` struct.

The coefficients themselves come from [`finite_difference_coefficients`](@ref); this constructor
only wraps them and builds the `Tridiagonal` view over `a_vec`/`e_vec` that root-finding and inverse
iteration mutate in place.
"""
function AcousticProblemCache(env::UnderwaterEnv, props::AcousticProblemProperties)
    a_vec, e_vec, λ_scaling = finite_difference_coefficients(env, props)
    # `a_vec` is shared with the matrix, not copied into it — that is what lets
    # `create_finite_diff_matrix!` update `A` by writing to `cache.a_vec`.
    A = Tridiagonal(e_vec[2:end], a_vec, e_vec[2:end])
    return AcousticProblemCache(a_vec, e_vec, λ_scaling, A)
end

function Base.show(io::IO, ::AcousticProblemCache{T}) where {T}
    return print(io, "AcousticProblemCache{$T}")
end

### Bisection and Sturm's Sequence

# Function to scale the Sturm sequence
function scale_const(p1, p2, Φ=1e20, Γ=1e-20)
    w = max(abs(p1), abs(p2))
    if w > Φ
        return Γ
    elseif 0 < w < Γ
        return Φ
    else
        return 1.0
    end
end

"""
    det_sturm(kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; scale=true)

Evaluate the Sturm sequence of the finite-difference system at trial wavenumber `kr`.

Returns `(det, mode_num)`: the determinant of the tridiagonal system, and the number of modes with
wavenumber above `kr` (the count decreases as `kr` sweeps up through the trapped band). `scale=true`
rescales the sequence whenever it would overflow or underflow; the factor is piecewise constant in
both `kr` and the environment parameters, so it cancels in any derivative of the root and does not
affect the mode count.

Only defined for `kr >= 2π * freq / cb` — see [`get_g`](@ref).
"""
function det_sturm(kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; scale=true)
    local p2, p1, p0, λ
    mode_count = 0
    g = get_g(kr, env, props)

    # Calculate the Sturm Sequence.
    k = 1
    p0 = 0.0
    p1 = 1.0
    for i in eachindex(props.Nz_vec)
        Nz = props.Nz_vec[i]
        for j in 1:Nz
            a = cache.a_vec[k]
            e = cache.e_vec[k]
            λ = kr^2 * cache.λ_scaling[k]
            k += 1
            # If we reached the last element of the last layer
            if (i == length(props.Nz_vec)) && (j == Nz)
                p2 = (λ - (0.5 * a - g)) * p1 - e^2 * p0
                if scale
                    s = scale_const(p1, p2)
                    p1 *= s
                    p2 *= s
                end
                if p1 * p2 < 0
                    mode_count += 1
                end
            else
                # Else, we're in the middle of the layers
                p2 = (λ - a) * p1 - e^2 * p0
                # rescale the sequence
                if scale
                    s = scale_const(p1, p2)
                    p1 *= s
                    p2 *= s
                end
                # count the modes
                if p1 * p2 < 0
                    mode_count += 1
                end
                p0 = p1
                p1 = p2
            end
        end
    end
    return (det=p2, mode_num=mode_count)
end

"""
    bisection(env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache)

Bisection method to find the intervals where the roots (wavenumbers) lie.

Searches `kr ∈ [ω/(0.9999 cb), max(ω/c)]`, the band in which modes are trapped. Returns an
`n_modes × 2` matrix of `[left right]` brackets, or `nothing` when the environment supports no
trapped modes at this frequency (e.g. water too shallow relative to the wavelength).
"""
function bisection(env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache)
    ω = 2pi * props.freq
    kr_max = maximum(ω ./ env.c.c)
    kr_min = ω / (0.9999 * env.cb)  # multiplying by 0.9999 so I don't touch the wavenumber boundary when root finding
    kr_min, kr_max = promote(kr_min, kr_max)
    n_max = last(det_sturm(kr_min, env, props, cache))
    if n_max == 0
        return nothing
    end
    n_min = last(det_sturm(kr_max, env, props, cache))

    # Initialize arrays
    kLeft = fill(kr_min, n_max + 1)
    kRight = fill(kr_max, n_max + 1)

    # Main loop
    k1 = kr_min
    k2 = kr_max
    if n_max > 1
        for mm in 1:(n_max - 1)
            # ii = 0
            if kLeft[mm] == kr_min
                k2 = kRight[mm]
                k1 = max(maximum(kLeft[(mm + 1):end]), kr_min)

                for _ in 1:50
                    # ii += 1
                    kmid = sqrt(mean([k1^2, k2^2]))
                    nmid = last(det_sturm(kmid, env, props, cache))
                    Δn = nmid - n_min

                    if Δn < mm
                        k2 = kmid
                        kRight[mm] = kmid
                    else
                        # Reaching here means Δn >= mm >= 1, so `kLeft[Δn]` is in bounds, and
                        # Δn <= n_max - n_min <= n_max, so `kRight[Δn + 1]` is too (both vectors
                        # have n_max + 1 entries). Keep that invariant in mind before changing the
                        # branch condition — indexing here is not otherwise guarded.
                        k1 = kmid
                        if kRight[Δn + 1] >= kmid
                            kRight[Δn + 1] = kmid
                        end
                        if kLeft[Δn] <= kmid
                            kLeft[Δn] = kmid
                        end
                    end

                    if kLeft[mm] != kr_min # if the the min wavenumber changed, we're done
                        # verbose && println("Mode $mm: Took $ii iterations")
                        break
                    end
                end
            end
        end
    end
    intervals = [kLeft[1:(end - 1)] kRight[1:(end - 1)]]
    if !isempty(intervals)
        # intervals[end, 1] += eps(kr_min) # to avoid solvers to get complex roots
        # intervals[1, 2] -= eps(kr_min) # to avoid solvers to get complex roots
    else
        println("Wavenumber intervals are empty!")
    end
    return intervals
end

"""
    find_kr(env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; method=ITP(), kwargs...)

Find the roots of the acoustic problem. `method` is any bracketing solver accepted by
`IntervalNonlinearProblem`. Returns an empty vector when `bisection` finds no trapped modes.
"""
function find_kr(
    env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; method=ITP(), kwargs...
)
    kr_spans = bisection(env, props, cache)
    if isnothing(kr_spans)
        return Vector{eltype(props.freq)}()
    end
    # `map` rather than filling a preallocated vector: this function is on the seam, and the
    # `setindex!` the loop used to do is exactly what reverse mode cannot follow. `bisection` is
    # `@non_differentiable`, so indexing its result carries no derivative either way.
    return map(i -> solve_for_kr(kr_spans[i, :], env, props, cache; method=method, kwargs...), axes(kr_spans, 1))
end

"""
Solve for the roots of the acoustic problem.
"""
function solve_for_kr(span, env, props, cache; method=ITP(), kwargs...)
    function f(u, p)
        return first(det_sturm(u, env, props, cache))
    end
    prob = IntervalNonlinearProblem{false}(f, span)
    sol = solve(prob, method; kwargs...)
    return sol.u
end

### Inverse Iteration

"""
    integral_trapz(y, x)

Trapezoidal integral of samples `y` taken at (not necessarily uniform) abscissae `x`.

Written out rather than delegated to `Integrals.SampledIntegralProblem` because it sits *on the
differentiable seam* — [`normalize_mode`](@ref) calls it, and reverse mode has to trace through it
to reach the density profile and the mesh coordinates. The formula is the same one the package's
`TrapezoidalRule` applies, and the two agree to the last bit on the mode normalizations here.
"""
function integral_trapz(y, x)
    return sum((x[2:end] .- x[1:(end - 1)]) .* (y[1:(end - 1)] .+ y[2:end]) ./ 2)
end

function create_finite_diff_matrix!(kr, env, props, cache)
    g = get_g(kr, env, props)

    # Update the diagonal elements.
    # Spelled `.= x .- y` rather than `.-=`: `@views` on an updating broadcast is a syntax error on
    # Julia 1.10 ("invalid let syntax"), which is the lower bound declared in Project.toml. Same
    # allocation behaviour — @views makes both sides views either way.
    cache.a_vec[end] = 0.5 * cache.a_vec[end] - kr^2 .* cache.λ_scaling[end] - g
    @views cache.a_vec[1:(end - 1)] .= cache.a_vec[1:(end - 1)] .- kr^2 .* cache.λ_scaling[1:(end - 1)]

    # The Tridiagonal matrix will automatically reflect these changes
    # since it's using views of the vectors
    return nothing
end

function return_finite_diff_matrix!(kr, env, props, cache)
    g = get_g(kr, env, props)
    cache.a_vec[end] = 2 * (cache.a_vec[end] + kr^2 .* cache.λ_scaling[end] + g)
    @views cache.a_vec[1:(end - 1)] .= cache.a_vec[1:(end - 1)] .+ kr^2 .* cache.λ_scaling[1:(end - 1)]  # see above
    # The Tridiagonal matrix will automatically reflect these changes
    # since it's using views of the vectors
    return nothing
end

"""
    mode_eigenvector(kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; reltol=0.01)

Inverse-iterate the finite-difference system at wavenumber `kr` for the mode shape it supports.

Returns `(kr_new, v)`: the refined wavenumber estimate produced by the iteration, and the
eigenvector, scaled to `‖v‖₂ = 1` and sign-fixed so that `v[1] ≥ 0`. Physical normalization is
[`normalize_mode`](@ref)'s job.

This is the second of the two functions *behind the differentiable seam* (see the note at the top of
this file): it mutates the cache, branches on `argmax`, and iterates, and reverse mode never traces
any of it — `src/kraken_ad.jl` states the eigenvector adjoint directly instead. The cache is left
exactly as it was found.
"""
function mode_eigenvector(
    kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; reltol=0.01
)
    local kr_new, w0, w1
    N = sum(props.Nz_vec)
    # Initialization
    w0 = normalize(ones(eltype(kr), N))
    w1 = similar(w0)
    # Create the finite-difference matrix
    kr_try = kr - 1e3 * eps(kr)
    create_finite_diff_matrix!(kr_try, env, props, cache) # Generate the tridigonal finite-diff matrix with the new kr
    # Inversete iteration
    for ii in 1:50 # We typically don't need more than 50 iterations
        w1 .= cache.A \ w0
        m = argmax(abs.(w1))
        kr_new = w0[m] / w1[m] + kr_try
        normalize!(w1)
        if relative_error(w0, w1) < reltol # Default is 1% relative tolerance
            w0 .= w1
            break # If the relative error is small enough, we're done with inverse iteration
        end
        w0 .= w1
    end
    # Inverse iteration complete
    w0 = ifelse(w0[1] < 0, w0 .* -1, w0) # Ensure the first element is positive for consistency between modes
    return_finite_diff_matrix!(kr_try, env, props, cache) # Reset the cache
    return kr_new, w0
end

"""
    normalize_mode(v, kr, env::UnderwaterEnv, props::AcousticProblemProperties)

Scale a mode shape so that its total energy — water column plus bottom half-space — is one:

    ∫ v(z)² / ρ(z) dz  +  v(D)² / (2 ρb √(kr² - (ω/cb)²))  =  1

`v` is [`mode_eigenvector`](@ref)'s unit-2-norm eigenvector; the scaling is what turns it into the
mode function the field sum expects.

Kept *on the differentiable seam* deliberately. The normalization is the only part of a mode shape
that depends on the density profile and the mesh coordinates other than through the finite-difference
coefficients, and leaving it traced means the backend reaches `env.ρ` through the same
[`density`](@ref) call it already differentiates in [`finite_difference_coefficients`](@ref), rather
than the eigenvector rule having to carry an interpolant adjoint of its own.
"""
function normalize_mode(v, kr, env::UnderwaterEnv, props::AcousticProblemProperties)
    zn = reduce(vcat, props.zn_vec)
    ρn = density(env.ρ, zn)
    amp1 = integral_trapz(abs2.(v) ./ ρn, zn) # Amplitude of the waveguide
    amp2 = v[end]^2 / (2 * env.ρb * sqrt(kr^2 - (2pi * props.freq / env.cb)^2)) # Same for the bottom half-space
    return v ./ sqrt(amp1 + amp2)
end

"""
    inverse_iteration(kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; kwargs...)

Performs inverse iteration to find the corresponding modal depth function ψₘ for a given wavenumber kᵣ

Returns `(kr_new, ψ)`. Composed of [`mode_eigenvector`](@ref), which is opaque to reverse-mode AD and
carries its own rule, and [`normalize_mode`](@ref), which is traced.
"""
function inverse_iteration(
    kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; reltol=0.01
)
    kr_new, v = mode_eigenvector(kr, env, props, cache; reltol=reltol)
    return kr_new, normalize_mode(v, kr_new, env, props)
end

function inverse_iteration(
    kr_vec::Vector, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; kws...
)
    if isempty(kr_vec)
        return kr_vec, Matrix{eltype(kr_vec)}(undef, length(cache.a_vec), 0)
    end
    # `map` plus `reduce(hcat, ...)` rather than filling a preallocated matrix: this method is on the
    # seam, and writing columns into a buffer is exactly the mutation reverse mode cannot follow.
    solved = map(kr -> inverse_iteration(kr, env, props, cache; kws...), kr_vec)
    return first.(solved), reduce(hcat, last.(solved))
end

### Full KRAKEN solve with Richardson's Extrapolation
"""
    h_extrap_matrix(hs)

The Vandermonde matrix of the Richardson extrapolation: row `i` is `[1, h_i², h_i⁴, …]` for the mesh
spacings `hs`, square by construction (one column per level, one row per level).

Built by broadcasting rather than assembled row by row, because `kraken_jl` is on the differentiable
seam and writing rows into a preallocated matrix is exactly the mutation reverse mode cannot follow.
"""
h_extrap_matrix(hs) = hs .^ permutedims(0:2:(2 * length(hs) - 2))

function richard_extrap(h_meshes, krs_meshes)
    # Mesh spacings and wavenumbers do not have to carry the same parameters — differentiating with
    # respect to `cb` alone makes `h_meshes` a plain `Float64` matrix and `krs_meshes` a `Dual`
    # vector — so promote explicitly rather than handing a mixed-precision system to `LinearSolve`.
    T = promote_type(eltype(h_meshes), eltype(krs_meshes))
    sol = solve(LinearProblem(T.(h_meshes), T.(krs_meshes))).u
    return sqrt(sol[1])
end

"""
    kraken_jl(env, freq; n_meshes=5, rmax=10_000, method=ITP(), dont_break=false, abstol=1e-6, reltol=1e-6)

Solve the normal-mode problem for environment `env` at frequency `freq` (Hz). Top-level entry point.

Discretizes the depth-separated wave equation on a sequence of successively finer meshes, brackets
each mode with [`bisection`](@ref), refines it with [`solve_for_kr`](@ref), recovers the mode shapes
with [`inverse_iteration`](@ref), and Richardson-extrapolates the squared wavenumbers across mesh
levels until they converge.

Returns a `NormalModeSolution` with fields `kr` (horizontal wavenumbers, descending), `modes` (one
column per mode, sampled on the finest mesh), `env`, and `props`. Returns an empty solution when the
environment supports no trapped modes at this frequency.

# Keyword arguments
- `n_meshes`: maximum number of mesh refinement levels to extrapolate across.
- `rmax`: maximum range of interest (m). Sets the convergence criterion — refinement stops once the
  wavenumber change would shift phase by less than one radian over `rmax`.
- `method`: bracketing solver passed to `IntervalNonlinearProblem`.
- `dont_break`: run all `n_meshes` levels even after convergence. Useful for studying mesh error.
- `abstol`, `reltol`: tolerances forwarded to the root solver.

# Differentiability

`kraken_jl` is differentiable end to end in both directions — forward mode straight through, reverse
mode via the rules in `kraken_ad.jl`. Two things about the returned gradient are worth being explicit
about, because neither is an approximation that can be tightened away:

  * **It is conditional on the mesh schedule this call selected.** The number of refinement levels is
    chosen by the convergence test above, which is an integer decision in the parameters; it is
    wrapped in `ignore_derivatives` so the gradient is taken at *fixed* schedule. That is the correct
    and standard treatment, but it means the gradient is discontinuous at parameter values where the
    level count changes — approach such a point from either side and the two one-sided gradients
    differ by the extrapolation's own error term. Pass `dont_break=true` for a schedule that is fixed
    by construction if a smooth objective matters more than the cost.
  * **The mode shapes come from the coarsest mesh**, as they do in the primal — only the wavenumbers
    are extrapolated. So `∂modes/∂θ` is the derivative of the level-1 solve.

# Example
```julia
env = UnderwaterEnv(pekeris_env()...)
sol = kraken_jl(env, 100.0)
sol.kr        # 5 wavenumbers
sol.modes     # mode shapes, one column each
```
"""
function kraken_jl(env, freq; n_meshes=5, rmax=10_000, method=ITP(), dont_break=false, abstol=1e-6, reltol=1e-6)
    # convert frequency to float if needed
    if freq isa Int
        freq = float(freq)
    end
    # generate all problem properties for every mesh
    props = AcousticProblemProperties(env, freq)
    # The mesh spacings, one per refinement level, grown by `vcat` as the loop below adds levels.
    # This used to be an `n_meshes × n_meshes` matrix filled a row at a time; the rows are now
    # broadcast into shape by `h_extrap_matrix` only when the extrapolation needs them, because
    # `setindex!` on a matrix that outlives the expression is not traceable in reverse mode.
    h_list = [props.Δz_vec[1]]

    # First Mesh (i_power = 1)
    cache = AcousticProblemCache(env, props)
    krs = find_kr(env, props, cache; method=method, abstol=abstol, reltol=reltol)
    if isempty(krs)
        return NormalModeSolution(krs, Matrix{eltype(krs)}(undef, 0, 0), env, props)
    end

    # Inverse Iteration
    krs_coarse, modes = inverse_iteration(krs, env, props, cache)
    # If we want only one mesh calculation (n_meshes = 1), return the result
    if n_meshes == 1
        return NormalModeSolution(krs_coarse, modes, env, props)
    end

    # Richardson Extrapolation from here on out if n_mesh > 1
    # Initialize
    #TODO: reuse kr_coarse for initial value for root finding for higher meshes
    M = length(krs_coarse)
    # `vcat` rather than `push!`, for the same reason `h_list` above is not a preallocated matrix.
    # `n_meshes` is 5 by default, so rebuilding these two vectors per level is free.
    krs_all = [krs_coarse .^ 2]

    # Richardson's extrapolation
    krs_old = krs_coarse
    rich_krs = krs_coarse
    for i_power in 2:n_meshes
        factor = 2^(i_power - 1)
        props_new = AcousticProblemProperties(env, freq; factor=factor)
        cache = AcousticProblemCache(env, props_new)
        krs_new = find_kr(env, props_new, cache; method=method)
        if length(krs_new) < M
            M = length(krs_new)
        end
        h_list = vcat(h_list, [props_new.Δz_vec[1]])
        krs_all = vcat(krs_all, [krs_new .^ 2])
        h_meshes = h_extrap_matrix(h_list)
        rich_krs = map(mm -> richard_extrap(h_meshes, map(ii -> krs_all[ii][mm], 1:i_power)), 1:M)
        # Whether to refine again is a *discrete* decision — it selects a mesh schedule, it is not a
        # quantity anything downstream depends on smoothly. Under `ignore_derivatives` it builds no
        # tape, so the gradient is the one conditional on the schedule this run selected. See the
        # note in the docstring.
        stop = ignore_derivatives() do
            # Check if the difference is less than the tolerance
            errs = abs.(rich_krs[1:M] - krs_old[1:M])
            err = errs[round(Int, 2 * M / 3)] # apparently this is used in KRAKEN to check for convergence
            # If the difference is less than the tolerance, or we've reached the maximum number of
            # meshes, interpolate krs_meshes with h_meshes and return the result
            return !dont_break && err * rmax < 1
        end
        stop && break
        krs_old = krs_new
    end

    return NormalModeSolution(rich_krs[1:M], modes[:, 1:M], env, props)
end

struct NormalModeSolution{T1,T2}
    kr::T1
    modes::T2
    env::UnderwaterEnv
    props::AcousticProblemProperties
    function NormalModeSolution(kr, modes, env, props)
        return new{typeof(kr),typeof(modes)}(kr, modes, env, props)
    end
end

function Base.show(io::IO, ρint::NormalModeSolution{T1,T2}) where {T1,T2}
    return print(io, "NormalModeSolution{", eltype(T1), "}(", length(ρint.kr), " modes)")
end

### Helper functions
function relative_error(v1, v2)
    return mean(abs.((v1 .- v2) ./ v1))
end
