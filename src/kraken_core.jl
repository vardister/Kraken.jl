# This file contains the core types and functions for the Kraken package.
# 
# The main types are:
# - `UnderwaterEnv`: Contains the sound speed profile and density profile of the underwater environment.
# - `AcousticProblemProperties`: Contains the properties of the acoustic problem.
# 
# The main functions are:
# - `prepare_vectors`: Prepares the vectors for the acoustic problem.
# - `bisection`: Finds the intervals for the roots of the acoustic problem.
# - `solve_for_kr`: Solves for the roots of the acoustic problem.
# - `inverse_iteration`: Solves for the modes of the acoustic problem.
# - `to_diff`: Converts the acoustic problem into a differentiable function.
# - `find_cg`: Finds the group velocity of the acoustic problem.
# 
### Load packages
using LinearAlgebra
using Statistics
using NonlinearSolve
using UnPack
using Integrals
using LinearSolve
using DataInterpolations: LinearInterpolation, CubicSpline


# For debuggina and developement
using ProtoStructs
using Infiltrator


### Main Types

### Sound Speed Profile
abstract type SoundSpeedProfile end
abstract type SampledSSP <: SoundSpeedProfile end

"""
Sound speed profile based on measurements at discrete depths.
"""
struct SampledSSP1D{T1, T2, T3} <: SampledSSP
    z::Vector{T1}
    c::Vector{T2}
    interp::Symbol
    f::T3
    function SampledSSP1D(depth, c, interp::Symbol)
        f = interp === :smooth ? CubicSplineInterpolation(c, depth; extrapolate=true) :
            interp === :linear ? LinearInterpolation(c, depth; extrapolate=true) :
            throw(ArgumentError("Unknown interpolation. Choose from :linear or :smooth"))
        return new{eltype(depth), eltype(depth), typeof(f)}(-depth, c, interp, f)
    end
end

"""
Constructor for `SampledSSP1D`.
"""
SampledSSP(depth, c) = SampledSSP1D(depth, c, :linear)
SampledSSP(depth, c, interp::Symbol) = SampledSSP1D(depth, c, interp)


function Base.show(io::IO, ρint::SampledSSP1D{T1,T2,T3}) where {T1,T2,T3}
    print(io, "SampledSSP1D{", T1, ",", T2, ",", ρint.interp, "}(", length(ρint.z), " points)")
end


### Density Profile
abstract type DensityProfile end
abstract type SampledDensity <: DensityProfile end

"""
Density profile based on measurements at discrete depths.

# Fields
- `z::Vector`: Depths at which the density is measured.
- `ρ::Vector`: Density values at the depths.
- `interp::Symbol`: Interpolation method used.
- `f::Function`: Interpolation function.
"""
struct SampledDensity1D{T1, T2, T3} <: SampledDensity
    z::Vector{T1}
    ρ::Vector{T2}
    interp::Symbol
    f::T3
    function SampledDensity1D(depth, ρ, interp)
        f = interp === :smooth ? CubicSpline(ρ, depth; extrapolate=true) :
            interp === :linear ? LinearInterpolation(ρ, depth; extrapolate=true) :
            throw(ArgumentError("Unknown interpolation. Choose from :linear or :smooth"))
        return new{eltype(depth), eltype(ρ), typeof(f)}(-depth, ρ, interp, f)
    end
end

SampledDensity(depth, ρ) = SampledDensity1D(depth, ρ, :linear)
SampledDensity(depth, ρ, interp::Symbol) = SampledDensity1D(depth, ρ, interp)

function Base.show(io::IO, ρint::SampledDensity1D{T1,T2,T3}) where {T1,T2,T3}
    print(io, "SampledDensity1D{", T1, ",", T2, ",", ρint.interp, "}(", length(ρint.z), " points)")
end


### Underwater Environment

"""
Underwater environment containing the sound speed profile and density profile.
"""
struct UnderwaterEnv{SSPType<:SoundSpeedProfile, DensityType <: DensityProfile}
    c::SSPType
    ρ::DensityType
    ρb::Float64
    cb::Float64
    h_vec::Vector{Float64}
    layer_depth::Vector{Float64}
    depth::Float64

    function UnderwaterEnv(ssp, layers, sspHS)
        c = SampledSSP(ssp[:, 1], ssp[:, 2])
        ρ = SampledDensity(ssp[:, 1], ssp[:, 3])
        ρb = sspHS[2, 3]
        cb = sspHS[2, 2]
        layer_thickness = get_thickness(layers)
        layer_depth = layers[:, 3]
        depth = sspHS[2, 1]
        return new{typeof(c), typeof(ρ)}(c, ρ, ρb, cb, layer_thickness, layer_depth, depth)
    end
end

### Functions that convert SSP information (similar to KRAKEN) to environment and problem structs
get_thickness(layers::Matrix{Float64}) = abs.(accumulate(-, layers[:, 3]))

function get_Nz_vec(env, ω; n_per_wavelength = 20, ipower = 1)
    @unpack h_vec, cb, c = env
    all_c = vcat(c.c, cb)
    kr_max = maximum(ω ./ all_c)
    Nz_vec = Int[]
    Δz_vec = Float64[]

    for h in h_vec
        Δz_power = ipower # Initial mesh multiplier
        n_per_wavelength = n_per_wavelength # The number of depth mesh points in a wavelength
        Lmin = 2π / kr_max # The lowest wavelength available in the problem
        # 20 points per wavelength. h_power is for richardson extrapolation
        Δz = (Lmin / n_per_wavelength)^Δz_power
        Nz = ceil(Int, h / Δz)
        Δz_new = h / Nz
        push!(Nz_vec, Nz)
        push!(Δz_vec, Δz_new)
    end
    return Nz_vec, Δz_vec
end

function get_z_vec(env::UnderwaterEnv, Nz_vec, Δz_vec)
    @unpack layer_depth = env
    zn_all = Vector{Vector{Float64}}()
    z0 = 0.0
    for (i, Nz) in enumerate(Nz_vec)
        Δz = Δz_vec[i]
        z_layer = layer_depth[i]
        zn = range(z0 + Δz, z_layer, Nz)
        push!(zn_all, zn)
        z0 = z_layer
    end
    return zn_all
end

struct AcousticProblemProperties{T<:Real}
    freq::T
    ω::T
    Nz_vec::Vector{Int}
    Δz_vec::Vector{Float64}
    zn_vec::Vector{Vector{Float64}}

    function AcousticProblemProperties(env::UnderwaterEnv, freq; ipower::Int = 1, n_per_wavelength = 20)
        ω = 2π * freq
        Nz_vec, Δz_vec = get_Nz_vec(env, ω; ipower = ipower, n_per_wavelength=n_per_wavelength)
        zn_vec = get_z_vec(env, Nz_vec, Δz_vec)
        return new{eltype(freq)}(freq, ω, Nz_vec, Δz_vec, zn_vec)
    end
end


### Prepare vectors

a_element_norm(c, ρ, ω, h) = (-2 + h^2 * (ω / c)^2) / (h * ρ)
e_element(ρ, h) = 1 / (h * ρ)

function get_g(kr, env::UnderwaterEnv, props::AcousticProblemProperties)
    @unpack cb, ρb = env
    @unpack ω = props
    g = sqrt(kr^2 - (ω / cb)^2) / ρb
    return g
end

moving_average(vec, len) = [mean(vec[i:(i+len-1)]) for i in 1:(length(vec)-len+1)]

function prepare_vectors(env, props)
    @unpack c, ρ, ρb, cb = env
    @unpack ω, zn_vec, Δz_vec, Nz_vec= props
    Ntotal = sum(Nz_vec)
    Ni = prepend!(accumulate(+, Nz_vec), 0)
    a_vec = zeros(eltype(ω), Ntotal)
    e_vec = similar(a_vec)    
    scaling_factor = similar(a_vec)

    for i in eachindex(props.zn_vec)
        zn = props.zn_vec[i]
        Δz = props.Δz_vec[i]
        cn = env.c.f(zn)
        ρn = env.ρ.f(zn)

        a_vec[Ni[i]+1:Ni[i+1]] = a_element_norm.(cn, ρn, ω, Δz)
        e_vec[Ni[i]+1:Ni[i+1]] = e_element.(env.ρ.f(zn), Δz)

        scaling_factor[Ni[i]+1:Ni[i+1]] .= e_vec[Ni[i]+1:Ni[i+1]] .* (Δz ^ 2)
    end
    # Interface conditions between layers
    if length(props.zn_vec) > 1
        for i in 1:(length(props.zn_vec)-1)
            a_vec[Nz_vec[i]] = 0.5 * (a_vec[Nz_vec[i]] + a_vec[Nz_vec[i]+1])
        end
    end
    # λ scaling
    scaling_factor = moving_average(scaling_factor, 2)
    λ_scaling = append!(scaling_factor, e_vec[end] * props.Δz_vec[end]^2 / 2)
    return a_vec, e_vec, λ_scaling
end

### Bisection and Sturm's Sequence
"""
Calculate the Sturm sequence for the acoustic problem.
"""

function det_sturm(kr, env, props, a_vec, e_vec, λ_scaling; stop_at_k = nothing)
    # If A is 1x1, no need to calculate determinant
    mode_count = 0
    if length(a_vec) == 1
        return mode_count, 0.0
    end

    @unpack Nz_vec, Δz_vec, zn_vec = props
    @unpack ρ = env
    g = get_g(kr, env, props)

    local p2, p1, p0
    # Calculate the Sturm Sequence.
    k = 1
    p0 = 0.0
    p1 = 1.0
    for i in eachindex(Nz_vec)
        Nz = Nz_vec[i]

        for j in 1:Nz
            a = a_vec[k]
            e = e_vec[k]
            λ = kr^2 * λ_scaling[k]
            k += 1

            # If we reached the last element of the last layer
            if i == length(Nz_vec) && j == Nz
                p2 = (λ - (0.5 * a - g)) * p1 - e^2 * p0
                s = scale_const(p1, p2)
                p1 *= s
                p2 *= s
                if p1 * p2 < 0
                    mode_count += 1
                end
                return mode_count, p2
            end

            # Else, we're in the middle of the layers
            p2 = (λ - a) * p1 - e^2 * p0
            # rescale the sequence
            s = scale_const(p1, p2)
            p1 *= s
            p2 *= s
            # count the modes
            if p1 * p2 < 0
                mode_count += 1
            end
            p0 = p1
            p1 = p2
            if stop_at_k !== nothing && k == stop_at_k
                return mode_count, p2
            end
        end
    end
    return mode_count, p2
end

"""
Solve the acoustic problem using bisection method.
"""
function bisection(env, props, a_vec, e_vec, λ_scaling; verbose = false)
    @unpack ω, zn_vec = props
    @unpack c, cb = env
    kr_max = maximum(ω ./ c.c)
    kr_min = ω / cb
    n_max = first(det_sturm(kr_min, env, props, a_vec, e_vec, λ_scaling))
    n_min = first(det_sturm(kr_max, env, props, a_vec, e_vec, λ_scaling))

    # Initialize arrays
    kLeft = fill(kr_min, n_max + 1)
    kRight = fill(kr_max, n_max + 1)

    # Main loop
    k1 = kr_min
    k2 = kr_max

    for mm in 1:(n_max-1)
        ii = 0
        if kLeft[mm] == kr_min
            k2 = kRight[mm]
            k1 = max(maximum(kLeft[(mm+1):end]), kr_min)

            for _ in 1:50
                ii += 1
                kmid = sqrt(mean([k1^2, k2^2]))
                nmid = first(det_sturm(kmid, env, props, a_vec, e_vec, λ_scaling))
                Δn = nmid - n_min

                if Δn < mm
                    k2 = kmid
                    kRight[mm] = kmid
                else
                    k1 = kmid
                    if kRight[Δn+1] >= kmid
                        kRight[Δn+1] = kmid
                    end
                    if kLeft[Δn] <= kmid
                        kLeft[Δn] = kmid
                    end
                end

                if kLeft[mm] != kr_min # if the the min wavenumber changed, we're done
                    verbose && println("Mode $mm: Took $ii iterations")
                    break
                end
            end
        end
    end
    intervals = [kLeft[1:(end-1)] kRight[1:(end-1)]]
    return intervals
end

### Solve for kr

"""
Find the roots of the acoustic problem.
"""
function find_kr(env, props, a_vec, e_vec, λ_scaling; method = ITP(), kwargs...)
    @unpack ω = props
    kr_spans = bisection(env, props, a_vec, e_vec, λ_scaling)
    if isempty(kr_spans)
        return Vector{typeof(ω)}()
    end
    krs = solve_for_kr(kr_spans, env, props, a_vec, e_vec, λ_scaling; method = method, kwargs...)
    return krs
end

"""
Solve for the roots of the acoustic problem.
"""
function solve_for_kr(kr_spans, env, props, a_vec, e_vec, λ_scaling; method = ITP(), kwargs...)
    @unpack ω = props
    krs = ones(eltype(ω), size(kr_spans, 1))
    for (i, uspan) in enumerate(eachrow(kr_spans))
        prob = IntervalNonlinearProblem(
            (kr_try, mm) -> last(det_sturm(kr_try, env, props, a_vec, e_vec, λ_scaling)),
            uspan,
        )
        krs[i] = solve(prob, method, kwargs...).u  # sol.u is the solution itself
    end
    return krs
end


### Inverse Iteration
function inverse_iteration(kr, env, props, a_vec, e_vec, scaling; tol = 1e-3, verbose = false)
    @unpack Nz_vec, zn_vec = props
    zn = Iterators.flatten(zn_vec) # it's all iterators joined together
    zn = collect(zn)
    ρn = env.ρ.f(zn)
    N = sum(Nz_vec)
    kr_try = kr - 1e3 * eps(kr)
    λ_try = kr_try^2 .* scaling
    w0 = normalize(ones(eltype(kr), N))
    w1 = similar(w0)

    local kr_new
    g = get_g(kr_try, env, props)
    d_vec = similar(a_vec)
    d_vec[1:end-1] .= a_vec[1:(end-1)] .- λ_try[1:(end-1)]
    d_vec[end] = 0.5 * a_vec[end] - λ_try[end] - g

    A = SymTridiagonal(d_vec, e_vec[2:end])
    for ii in 1:200
        prob = LinearProblem(A, w0)
        w1 .= solve(prob, nothing).u
        # Improve the estimate of the wavenumber
        _, m = findmax(abs.(w1))
        kr_new = w0[m] / w1[m] + kr_try
        w1 ./= norm(w1)
        if norm(abs.(w1) .- abs.(w0)) < tol
            verbose && println("Took $ii iterations to converge")
            break
        end
        w0 .= w1
    end
    if w0[1] < 0
        w0 .*= -1
    end
    # normalize the mode
    amp1 = integral_trapz(abs2.(w0) ./ ρn, zn)
    amp2 = w0[end]^2 / (2 * env.ρb * sqrt(kr_new^2 - (props.ω / env.cb)^2))
    w0 .= w0 / sqrt(amp1 + amp2)
    prepend!(w0, 0.0)

    # restore the matrix to its original form
    # a_vec[1:(end-1)] .+= λ_try[1:(end-1)]
    # a_vec[end] = 2 * (a_vec[end] + λ_try[end] + g)

    return kr_new, w0
end

function inverse_iteration(kr_vec::Vector, env, props, a_vec, e_vec, scaling; kws...)
    zn = Iterators.flatten(props.zn_vec) # it's all iterators joined together
    zn = collect(zn)
    modes = zeros(eltype(kr_vec), length(zn) + 1, length(kr_vec))
    kr_vec_new = similar(kr_vec)
    for (i, kr) in enumerate(kr_vec)
        kr_vec_new[i], modes[:, i] = inverse_iteration(kr, env, props, a_vec, e_vec, scaling; kws...)
    end
    # @infiltrate
    return kr_vec, modes
end

### Full KRAKEN solve with Richardson's Extrapolation
h_extrap(h, Nh) = [h^pow for pow in 0:2:2Nh-2]
function kraken_jl(env, freq; n_meshes = 5, rmax = 10_000, method = ITP(), dont_break = false)
    # First mesh first
    props = AcousticProblemProperties(env, freq; ipower = 1)
    a_vec, e_vec, λ_scaling = prepare_vectors(env, props)
    krs = find_kr(env, props, a_vec, e_vec, λ_scaling)
    if isempty(krs)
        return krs, Matrix{eltype(krs)}(undef, 0, 0)
    end
    krs_new, modes_coarse = inverse_iteration(krs, env, props, a_vec, e_vec, λ_scaling)
    # If there is only one mesh, return the result
    if n_meshes == 1
        return krs_new, modes_coarse
    end

    # Richardson Extrapolation from here on out if n_mesh > 1
    # Initialize
    #TODO: reuse kr_coarse for initial value for root finding for higher meshes
    M = length(krs_new)
    rich_krs = zeros(eltype(freq), M)
    krs_meshes = zeros(eltype(krs_new), n_meshes, M)
    krs_meshes[1, :] = krs_new .^ 2
    h_meshes = zeros(eltype(krs_new), n_meshes, n_meshes)
    # Richardson's extrapolation
    h_meshes[1, :] .= h_extrap(props.Δz_vec[1], n_meshes)
    krs_old = krs_new
    for i_power in 2:n_meshes
        factor = 2^(i_power - 1)
        props = AcousticProblemProperties(env, freq; ipower = factor)
        a_vec, e_vec, λ_scaling = prepare_vectors(env, props)
        krs_new = find_kr(env, props, a_vec, e_vec, λ_scaling)
        # If the number of modes has changed, update M
        if length(krs_new) < M
            M = length(krs_new)
            krs_meshes = krs_meshes[:, 1:M] .^ 2
        end
        krs_meshes[i_power, :] = krs_new[1:M] .^ 2
        h_meshes[i_power, :] .= h_extrap(props.Δz_vec[1], n_meshes)
        # interpolate krs_meshes with h_meshes
        for i in 1:M
            # @infiltrate
            y = solve(LinearProblem(h_meshes[1:i_power, 1:i_power], krs_meshes[1:i_power, i])).u
            rich_krs[i] = sqrt(y[1])
        end
        # Check if the difference is less than the tolerance
        errs = abs.(rich_krs - krs_old)
        err = errs[round(Int, 2 * M / 3)] # apparently this is used in KRAKEN to check for convergence
        # println("Current tol: $err at mesh $i_power")
        # If the difference is less than the tolerance, or we've reached the maximum number of meshes
        # interpolate krs_meshes with h_meshes and return the result
        if !dont_break && err * rmax < 1
            break
        end
        krs_old = krs_new
    end

    return rich_krs, modes_coarse
end