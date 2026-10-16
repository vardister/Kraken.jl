using NonlinearSolve
using UnPack
using StaticArrays
# using Infiltrator
using DocStringExtensions
using Roots
import NaNMath as nm

export PekerisUnderwaterEnv, find_kr, get_modal_function, pressure_f, find_spans, get_modal_function_values

"""
    PekerisUnderwaterEnv(c1, cb, ρ1, ρb, depth)

The Pekeris waveguide in closed form: an isovelocity water layer over a homogeneous fluid half-space.
This is the separate, analytic solve path: [`find_kr`](@ref) on one of these solves the dispersion
relation directly instead of going through the finite-difference machinery, which makes it the
reference the general solver is checked against. It is lossless.

# Arguments
- `c1`, `ρ1`: water sound speed (m/s) and density (kg/m³).
- `cb`, `ρb`: half-space sound speed and density. `cb > c1` is required for trapped modes and
  asserted.
- `depth`: water depth (m).

The same waveguide for the finite-difference solver is [`pekeris_env`](@ref).

# Example
```jldoctest
julia> env = PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0);

julia> env.depth
100.0
```
"""
struct PekerisUnderwaterEnv{T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}
    c1::T1
    cb::T2
    ρ1::T3
    ρb::T4
    depth::T5

    function PekerisUnderwaterEnv(c1::T1, cb::T2, ρ1::T3, ρb::T4, depth::T5) where {T1,T2,T3,T4,T5}
        @assert cb > c1 "Bottom sound speed must be greater than water sound speed"
        return new{T1,T2,T3,T4,T5}(c1, cb, ρ1, ρb, depth)
    end
end

"""
    find_kr(env::PekerisUnderwaterEnv, freq; n_points=2_000, method=Roots.A42(), kwargs...)

The horizontal wavenumbers of the Pekeris waveguide at frequency `freq` (Hz), in descending order,
from its dispersion relation `tan(kz·D) + (ρb/ρ1)·kz/γ = 0`.

The trapped band `ω/cb < kr < ω/c1` is sampled at `n_points` points, every sign change from positive
to negative is taken as a bracket, and each is refined with `Roots.find_zero`. The dispersion
function also changes sign at the poles of `tan`, from negative to positive, which is why only one
direction counts. Two modes closer together than the sampling step would be missed; `n_points` is
what to raise if that ever happens.

# Arguments
- `env`: the [`PekerisUnderwaterEnv`](@ref).
- `freq`: frequency (Hz).
- `n_points`: samples of the trapped band used to bracket the roots.
- `method`: the `Roots.jl` method passed to `find_zero`.
- `kwargs...`: passed to `find_zero`.

Returns an empty vector when no mode is trapped.

# Example
The finite-difference solve of the same waveguide agrees to the digits shown:
```jldoctest
julia> env = PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0);

julia> round.(find_kr(env, 100.0); digits=4)
5-element Vector{Float64}:
 0.4179
 0.415
 0.41
 0.4029
 0.3938

julia> round.(kraken_jl(UnderwaterEnv(pekeris_env()...), 100.0).kr; digits=4) == round.(find_kr(env, 100.0); digits=4)
true
```
"""
function find_kr(env::PekerisUnderwaterEnv, f; n_points=2_000, method=Roots.A42(), kwargs...)
    T = promote_type(typeof(env.c1), typeof(env.cb), typeof(env.ρ1), typeof(env.ρb), typeof(env.depth), typeof(f))
    f = convert(T, f)
    ω = 2π * f
    kr_min, kr_max = extrema([ω / env.c1, ω / env.cb])
    p = SA[env.c1, env.cb, env.ρ1, env.ρb, env.depth]
    function func(kr)
        @. nm.tan(nm.sqrt((ω / env.c1)^2 - kr^2) * env.depth) +
            (env.ρb / env.ρ1) * nm.sqrt((ω / env.c1)^2 - kr^2) / nm.sqrt(kr^2 - (ω / env.cb)^2)
    end

    kr_try = range(kr_min + eps(kr_min), kr_max - eps(kr_min); length=n_points)
    kr0 = collect(find_spans(func, kr_try))
    # kr0 = [0.39380573744004627, 0.4028554311353545, 0.40996684021864305, 0.4149566134428289, 0.41790332991380486]
    # @show kr0
    if isempty(kr0)
        return Vector{eltype(kr_min)}()
    end

    # prob = NonlinearProblem{false}(func, kr0, p)
    sol = zeros(eltype(f), length(kr0))
    for i in eachindex(kr0)
        sol[i] = find_zero(kr -> func(kr), kr0[i], method; kwargs...)
    end
    return reverse(sol)
end

"""
    find_spans(func, x)

The consecutive pairs `(x[i], x[i+1])` across which `func` goes from positive to negative, the
brackets `find_kr(::PekerisUnderwaterEnv, freq)` refines. Returns a lazy generator of tuples, or an
empty vector when there are none.

# Arguments
- `func`: a scalar function.
- `x`: increasing sample points.

# Example
```jldoctest
julia> collect(find_spans(cos, 0:0.5:5))   # cos falls through zero at π/2 and rises at 3π/2
1-element Vector{Tuple{Float64, Float64}}:
 (1.5, 2.0)
```
"""
function find_spans(func, x)
    idxes = [i for i in 1:(length(x) - 1) if (func(x[i]) > 0 && func(x[i + 1]) < 0)]
    if isempty(idxes)
        return Vector{eltype(x)}()
    end
    return ((x[i], x[i + 1]) for i in idxes)
end

"""
    solve_pekeris_equation(span, func_to_min, p; method = ITP())

    Solve the Pekeris equation.

    # Arguments
    - `span::Vector{Real}`: The span to search.
    - `func_to_min::Function`: The function to minimize.
    - `p::Vector{Real}`: The parameters.
    - `method::NonlinearSolver=ITP()`: The method to use for the search.

    # Returns
    - `Real`: The solution.
"""
function solve_pekeris_equation(span, func_to_min, p; method=ITP())
    kr0 = mean(span)
    prob = NonlinearProblem(func_to_min, kr0, p)
    sol = solve(prob, TrustRegion())
    return sol.u
end

"""
    get_modal_function_values(env, krs, freq, zr, zs)

Evaluate the normalized Pekeris mode functions at a receiver depth and a source depth.

Returns `(modes_zr, modes_zs)`, two vectors with one value per wavenumber in `krs`. Below `depth`
the modes decay exponentially into the half-space.

# Arguments
- `env`: the [`PekerisUnderwaterEnv`](@ref).
- `krs`: its wavenumbers at `freq`, from `find_kr(env, freq)`.
- `freq`: frequency (Hz).
- `zr`, `zs`: receiver and source depths (m, positive downwards).

# Example
```jldoctest
julia> env = PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0);

julia> krs = find_kr(env, 100.0);

julia> ψr, ψs = get_modal_function_values(env, krs, 100.0, 0.0, 50.0);

julia> all(iszero, ψr)   # every mode vanishes at the pressure-release surface
true
```
"""
function get_modal_function_values(env, krs, freq, zr, zs)
    @unpack c1, cb, ρ1, ρb, depth = env
    ω = 2π * freq
    kzw = sqrt.((ω / c1)^2 .- krs .^ 2) # vertical wavenumber in the water column
    kzb = sqrt.(krs .^ 2 .- (ω / cb)^2) # vertical wavenumber in the bottom half-space
    # The amplitude `amplitude` of the modes
    amplitude = sqrt.(
        (4 .* kzb .* kzw .* ρ1 .* ρb) ./
        (2 .* kzw .* ρ1 .* sin.(depth .* kzw) .^ 2 .- kzb .* ρb .* (-2 .* depth .* kzw .+ sin.(2 .* depth .* kzw))),
    )
    # The mode function
    function Ψ(z, mode)
        return amplitude[mode] * sin(kzw[mode] * z) * (z .<= depth) +
               amplitude[mode] * sin(kzw[mode] * depth) * exp(-kzb[mode] * (z - depth)) * (z .> depth)
    end
    modes_zr = [Ψ(zr, mode) for mode in eachindex(krs)]
    modes_zs = [Ψ(zs, mode) for mode in eachindex(krs)]
    return modes_zr, modes_zs
end

"""
    get_modal_function(env, krs, freq)

Build the Pekeris mode function `Ψ(z, mode)` for wavenumbers `krs` at frequency `freq`.

Returns a closure, so it can be evaluated at any depth: sinusoidal in the water column and
exponentially decaying below the interface, normalized the same way as
[`get_modal_function_values`](@ref) (which is the right call when you only need a fixed pair of
depths).

# Arguments
- `env`: the [`PekerisUnderwaterEnv`](@ref).
- `krs`: its wavenumbers at `freq`, from `find_kr(env, freq)`.
- `freq`: frequency (Hz).

The closure takes a depth `z` (m, positive downwards) and a mode index.

# Example
```jldoctest
julia> env = PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0);

julia> krs = find_kr(env, 100.0);

julia> Ψ = get_modal_function(env, krs, 100.0);

julia> Ψ(50.0, 2) == last(get_modal_function_values(env, krs, 100.0, 0.0, 50.0))[2]
true
```
"""
function get_modal_function(env, krs, freq)
    @unpack c1, cb, ρ1, ρb, depth = env
    ω = 2π * freq
    kzw = sqrt.((ω / c1)^2 .- krs .^ 2) # vertical wavenumber in the water column
    kzb = sqrt.(krs .^ 2 .- (ω / cb)^2) # vertical wavenumber in the bottom half-space
    # The amplitude `amplitude` of the modes
    amplitude = sqrt.(
        (4 .* kzb .* kzw .* ρ1 .* ρb) ./
        (2 .* kzw .* ρ1 .* sin.(depth .* kzw) .^ 2 .- kzb .* ρb .* (-2 .* depth .* kzw .+ sin.(2 .* depth .* kzw))),
    )
    # The mode function
    function Ψ(z, mode)
        return amplitude[mode] * sin(kzw[mode] * z) * (z .<= depth) +
               amplitude[mode] * sin(kzw[mode] * depth) * exp(-kzb[mode] * (z - depth)) * (z .> depth)
    end
    return Ψ
end

"""
    pressure_f(env::PekerisUnderwaterEnv, krs, freq, r, zs, zr; t0=0.1, max_modes=Inf)

The complex pressure of the Pekeris waveguide at one frequency, range and depth pair, as a normal
mode sum. Returns zero at `freq = 0` or when there are no modes.

The result carries an extra phase `exp(2πi·freq·(r/min(c1, cb) - t0))`, which shifts the arrival
to `t0` seconds into a time window when the spectrum is inverse-transformed. At a single frequency
it changes the phase and leaves `abs` alone.

# Arguments
- `env`: the [`PekerisUnderwaterEnv`](@ref).
- `krs`: its wavenumbers at `freq`, from `find_kr(env, freq)`.
- `freq`: frequency (Hz).
- `r`: range (m).
- `zs`, `zr`: source and receiver depths (m, positive downwards).
- `t0`: where in the time window the arrival lands (s).
- `max_modes`: sum only the first this many modes.

# Example
```jldoctest
julia> env = PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0);

julia> krs = find_kr(env, 100.0);

julia> p = pressure_f(env, krs, 100.0, 5_000.0, 25.0, 50.0);

julia> p isa Complex, abs(pressure_f(env, krs, 100.0, 5_000.0, 25.0, 50.0; t0=0.0)) ≈ abs(p)
(true, true)
```
"""
function pressure_f(env::PekerisUnderwaterEnv, krs, freq, r, zs, zr; t0=0.1, max_modes=Inf)
    if freq == 0
        return 0.0 + 0.0im
    end
    if length(krs) == 0
        return 0.0 + 0.0im
    end
    @unpack c1, cb, ρ1, ρb, depth = env
    # `get_modal_function(env, krs, freq)` returns a *function* of (z, mode); the values at the
    # receiver and source depths come from `get_modal_function_values`, which is the 5-argument
    # method this call always wanted.
    ψ_zr, ψ_zs = get_modal_function_values(env, krs, freq, zr, zs)
    t_offset = r / min(c1, cb) - t0  # align window correctly in time
    n_modes = Int(min(max_modes, length(krs)))
    pf =
        -1im * exp(-1im * π / 4) / (ρ1 * sqrt(8π * r)) *
        sum([ψ_zs[mode] * ψ_zr[mode] * exp(-1im * krs[mode] * r) / sqrt(krs[mode]) for mode in range(1, n_modes)])
    return pf * exp(2im * pi * freq * t_offset)
end
