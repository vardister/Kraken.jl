### Load packages
using LinearAlgebra
using Statistics
using Roots
using LinearSolve
using UnPack
using Integrals
using DataInterpolations: LinearInterpolation, CubicSpline
import NaNMath as nm

### Docs
using DocStringExtensions

## Debugging
# using Infiltrator

# Exports
export SampledSSP, SampledDensity
export UnderwaterEnv, AcousticProblemProperties, UnderwaterEnvFORTRAN
export AcousticProblemCache, bisection, solve_for_kr, inverse_iteration, det_sturm, kraken_jl, find_kr, get_g

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
SampledSSP(depth, c) = SampledSSP1D(depth, c, LinearInterpolation)
SampledSSP(depth, c, type::Symbol) = SampledSSP1D(depth, c, type)

function Base.show(io::IO, ρint::SampledSSP1D{T1,T2,T3}) where {T1,T2,T3}
    return print(io, "SampledSSP1D{", T1, ",", T2, ",", ρint.type, "}(", length(ρint.z), " points)")
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
SampledDensity(depth, ρ) = SampledDensity1D(depth, ρ, LinearInterpolation)
SampledDensity(depth, ρ, type::Symbol) = SampledDensity1D(depth, ρ, type)

function Base.show(io::IO, ρint::SampledDensity1D{T1,T2,T3}) where {T1,T2,T3}
    return print(io, "SampledDensity1D{", T1, ",", T2, ",", ρint.type, "}(", length(ρint.z), " points)")
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
    depth = krak_ssp.ssp[end, 1]
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

soundspeed(ssp::SampledSSP1D, z) = ssp.f(z)

"""
	maxsoundspeed(ssp::SoundSpeedProfile)

Get the maximum sound speed from the sound speed profile.
"""
function maxsoundspeed end

maxsoundspeed(ssp::SampledSSP1D) = maximum(ssp.c)

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
    a = vcat(0.0, layers[:, 3]...)
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
    Nz_vec = zeros(Int, length(env.h_vec))
    Δz_vec = zeros(eltype(env.h_vec), length(env.h_vec))

    for (i, h) in enumerate(env.h_vec)
        Lmin = 2π / kr_max # The lowest wavelength available in the problem
        # 20 points per wavelength. h_power is for richardson extrapolation
        Δz = (Lmin / n_per_wavelength)
        Nz = ceil(Int, h / Δz) * factor
        Nz = max(10, Nz) # Minimum of 10 points
        Δz_new = h / Nz
        Nz_vec[i] = Nz
        Δz_vec[i] = Δz_new
    end
    return Nz_vec, Δz_vec
end

"""
	get_z_vec(env::UnderwaterEnv, Nz_vec, Δz_vec)

Get the depth vector for each layer of the underwater environment according to the number of mesh points `Nz_vec`
 and mesh spacing `Δz_vec`.
"""
function get_z_vec(env::UnderwaterEnv, Nz_vec, Δz_vec)
    zn_all = Vector{typeof(env.layer_depth)}(undef, length(Nz_vec))
    z0 = 0.0
    for (i, Nz) in enumerate(Nz_vec)
        Δz = Δz_vec[i]
        z_layer = env.layer_depth[i]
        zn = range(z0 + Δz, z_layer, Nz)
        zn_all[i] = collect(zn)
        z0 = z_layer
    end
    return zn_all
end

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
"""
function get_g(kr, env::UnderwaterEnv, props::AcousticProblemProperties)
    g = sqrt(kr^2 - (2pi * props.freq / env.cb)^2) / env.ρb
    return g
end

### Bisection and Sturm's Sequence
function moving_average!(vs, n)
    vs[1:(end - n + 1)] .= [sum(@view vs[i:(i + n - 1)]) / n for i in 1:(length(vs) - (n - 1))]
    return nothing
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
	AcousticProblemCache(env::UnderwaterEnv, props::AcousticProblemProperties)

Prepare the vectors `a_vec`, `e_vec`, and `scaling_factor` for the acoustic problem. Return an `AcousticProblemCache` struct.
"""
function AcousticProblemCache(env::UnderwaterEnv, props::AcousticProblemProperties)
    Ntotal = sum(props.Nz_vec)
    Ni = prepend!(accumulate(+, props.Nz_vec), 0)
    #TODO: create vector that generalizes well to different types
    T = promote_type(eltype(env.c.c), eltype(env.ρ.ρ), typeof(env.cb), typeof(props.freq))
    a_vec = zeros(T, Ntotal)
    e_vec = similar(a_vec)
    scaling_factor = similar(a_vec)
    for i in eachindex(props.zn_vec)
        zn = props.zn_vec[i]
        Δz = props.Δz_vec[i]
        cn = soundspeed(env.c, zn)
        ρn = density(env.ρ, zn)

        a_vec[(Ni[i] + 1):Ni[i + 1]] .= a_element(cn, ρn, props.freq, Δz)
        e_vec[(Ni[i] + 1):Ni[i + 1]] .= e_element(ρn, Δz)

        scaling_factor[(Ni[i] + 1):Ni[i + 1]] .= e_vec[(Ni[i] + 1):Ni[i + 1]] .* (Δz^2)
    end
    # Interface conditions between layers
    if length(props.zn_vec) > 1
        loc = 0
        for i in 1:(length(props.zn_vec) - 1)
            loc += props.Nz_vec[i]
            a_vec[loc] = 0.5 * (a_vec[loc] + a_vec[loc + 1])
        end
    end

    moving_average!(scaling_factor, 2)
    scaling_factor[end] = e_vec[end] * props.Δz_vec[end]^2 / 2
    # Construct the Tridiagonal matrix
    A = Tridiagonal(e_vec[2:end], a_vec, e_vec[2:end])
    return AcousticProblemCache(a_vec, e_vec, scaling_factor, A)
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
	det_sturm(
		kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache;
		stop_at_k = nothing, return_det = false, scale = true)
"""
function det_sturm(
    kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; stop_at_k=nothing, scale=true
)
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
                if stop_at_k !== nothing && k == stop_at_k
                    p2, mode_count
                end
            end
        end
    end
    return (det=p2, mode_num=mode_count)
end

"""
	bisection(env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache)

Bisection method to find the intervals where the roots (wavenumbers) lie.
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
	find_kr(env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; method=Roots.A42, kwargs...)

Find the roots of the acoustic problem.
"""
function find_kr(
    env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; method=Roots.A42(), kwargs...
)
    kr_spans = bisection(env, props, cache)
    if isnothing(kr_spans)
        return Vector{eltype(props.freq)}()
    end
    krs = zeros(eltype(kr_spans), size(kr_spans, 1))
    i = 1
    for span in eachrow(kr_spans)
        krs[i] = solve_for_kr(span, env, props, cache; method=method, kwargs...)[1]
        i += 1
    end
    return krs
end

"""
Solve for the roots of the acoustic problem.
"""
function solve_for_kr(span, env, props, cache; method=Roots.A42(), kwargs...)
    function f(u, p=nothing)
        return first(det_sturm(u, env, props, cache))
    end
    sol = find_zero(f, span, method)
    return sol
end

### Inverse Iteration

function integral_trapz(y, x)
    problem = SampledIntegralProblem(y, x)
    method = TrapezoidalRule()
    return solve(problem, method).u
end


function create_finite_diff_matrix!(kr, env, props, cache)
    g = get_g(kr, env, props)

    # Update the diagonal elements
    cache.a_vec[end] = 0.5 * cache.a_vec[end] - kr^2 .* cache.λ_scaling[end] - g
    @views cache.a_vec[1:(end - 1)] .-= kr^2 .* cache.λ_scaling[1:(end - 1)]

    # The Tridiagonal matrix will automatically reflect these changes
    # since it's using views of the vectors
    return nothing
end

function return_finite_diff_matrix!(kr, env, props, cache)
    g = get_g(kr, env, props)
    cache.a_vec[end] = 2 * (cache.a_vec[end] + kr^2 .* cache.λ_scaling[end] + g)
    @views cache.a_vec[1:(end - 1)] .+= kr^2 .* cache.λ_scaling[1:(end - 1)]
    # The Tridiagonal matrix will automatically reflect these changes
    # since it's using views of the vectors
    return nothing
end

"""
    inverse_iteration(kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; kwargs...)

Performs inverse iteration to find the corresponding modal depth function ψₘ for a given wavenumber kᵣ
"""
function inverse_iteration(
    kr, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; reltol=0.01
)
    local kr_new, w0, w1, amp1, amp2
    zn = vcat(props.zn_vec...)
    ρn = density(env.ρ, zn)
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
    amp1 = integral_trapz(abs2.(w0) ./ ρn, zn) # Amplitude calculation for the waveguide
    amp2 = w0[end]^2 / (2 * env.ρb * sqrt(kr_new^2 - (2pi * props.freq / env.cb)^2)) # Same but for bottom half-space
    w0 ./= sqrt(amp1 + amp2) # Normalize the mode function
    return_finite_diff_matrix!(kr_try, env, props, cache) # Reset the cache
    return kr_new, w0
end

function inverse_iteration(
    kr_vec::Vector, env::UnderwaterEnv, props::AcousticProblemProperties, cache::AcousticProblemCache; kws...
)
    # Initialize containers
    modes = zeros(eltype(kr_vec), length(cache.a_vec), length(kr_vec))
    kr_vec_new = similar(kr_vec)
    # Loop over the wavenumbers
    for (ii, kr) in enumerate(kr_vec)
        kr_vec_new[ii], modes[:, ii] = inverse_iteration(kr, env, props, cache; kws...)
    end
    # Return the new kr_vec and modes
    return kr_vec_new, modes
end

### Full KRAKEN solve with Richardson's Extrapolation
h_extrap(h, Nh) = [h^pow for pow in 0:2:(2Nh - 2)]

function richard_extrap(h_meshes, krs_meshes)
    sol = solve(LinearProblem(h_meshes, krs_meshes)).u
    return sqrt(sol[1])
end

function kraken_jl(env, freq; n_meshes=5, rmax=10_000, method=A42(), dont_break=false, abstol=1e-6, reltol=1e-6)
    # First mesh first
    local rich_krs
    # convert frequency to float if needed
    if freq isa Int
        freq = float(freq)
    end
    # generate all problem properties for every mesh
    props = AcousticProblemProperties(env, freq)
    h_meshes = zeros(eltype(env.c.c), n_meshes, n_meshes)
    h_meshes[1, :] = h_extrap(props.Δz_vec[1], n_meshes)

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
    krs_all = Vector{Vector{eltype(krs_coarse)}}()
    push!(krs_all, krs_coarse .^ 2)

    # Richardson's extrapolation
    krs_old = krs_coarse
    for i_power in 2:n_meshes
        factor = 2^(i_power - 1)
        props_new = AcousticProblemProperties(env, freq; factor=factor)
        h_meshes[i_power, :] = h_extrap(props_new.Δz_vec[1], n_meshes)
        cache = AcousticProblemCache(env, props_new)
        krs_new = find_kr(env, props_new, cache; method=method)
        if length(krs_new) < M
            M = length(krs_new)
        end
        push!(krs_all, krs_new .^ 2)
        rich_krs = [richard_extrap(h_meshes[1:i_power, 1:i_power], [krs_all[ii][mm] for ii in 1:i_power]) for mm in 1:M]
        # Check if the difference is less than the tolerance
        errs = abs.(rich_krs[1:M] - krs_old[1:M])
        err = errs[round(Int, 2 * M / 3)] # apparently this is used in KRAKEN to check for convergence
        # If the difference is less than the tolerance, or we've reached the maximum number of meshes
        # interpolate krs_meshes with h_meshes and return the result
        if !dont_break && err * rmax < 1
            break
        end
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