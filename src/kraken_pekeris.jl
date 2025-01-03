using NonlinearSolve
using UnPack
using StaticArrays
using Infiltrator
using DocStringExtensions
import NaNMath as nm

export PekerisUnderwaterEnv, find_kr, get_modal_function, pressure_f, find_spans

"""
	PekerisUnderwaterEnv

	A type representing the underwater environment in the Pekeris model.

	# Fields
	- `c1::T1`: The sound speed in the water column.
	- `cb::T1`: The sound speed in the bottom half-space.
	- `ρ1::T1`: The density of the water column.
	- `ρb::T1`: The density of the bottom half-space.
	- `depth::T1`: The depth of the water column.

	# Constructor
	```julia
	PekerisUnderwaterEnv(c1::T1, cb, ρ1, ρb, depth) where {T1}
	```

	# Example
	```julia
	env = PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 2000.0, 100.0)
	```
"""
struct PekerisUnderwaterEnv{T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real}
	c1::T1
	cb::T2
	ρ1::T3
	ρb::T4
	depth::T5

	function PekerisUnderwaterEnv(c1::T1, cb::T2, ρ1::T3, ρb::T4, depth::T5) where {T1, T2, T3, T4, T5}
		@assert cb > c1 "Bottom sound speed must be greater than water sound speed"
		new{T1, T2, T3, T4, T5}(c1, cb, ρ1, ρb, depth)
	end
end

"""
    find_kr(env, freq, p; n_points = 5_000, method = NewtonRaphson())

    Find the horizontal wavenumbers for the Pekeris model.

    # Arguments
    - `env::PekerisUnderwaterEnv`: The underwater environment.
    - `freq::Real`: The frequency of the signal.
    - `p::Vector{Real}`: The Pekeris parameters.
    - `n_points::Int=2_000`: The number of points to use in the search.
    - `method::NonlinearSolver=NewtonRaphson()`: The method to use for the search.

    # Returns
    - `Vector{Real}`: The horizontal wavenumbers.
"""
function find_kr(env::PekerisUnderwaterEnv, freq, p; n_points = 2_000, method = NewtonRaphson())
	ω = 2π * freq
	kr_min, kr_max = extrema([ω / env.c1, ω / env.cb])
	# p = SA[env.c1, env.cb, env.ρ1, env.ρb, env.depth]
	func(kr, p) =
		@. nm.tan(nm.sqrt((ω / p[1])^2 - kr^2) * p[5]) + (p[4] / p[3]) * nm.sqrt((ω / p[1])^2 - kr^2) / nm.sqrt(kr^2 - (ω / p[2])^2)

	kr_try = range(kr_min + eps(kr_min), kr_max - eps(kr_min); length = n_points)
	kr0 = find_spans(func, kr_try, p)
	# kr0 = [0.39380573744004627, 0.4028554311353545, 0.40996684021864305, 0.4149566134428289, 0.41790332991380486]
	# @show kr0
	if isempty(kr0)
		return Vector{eltype(kr_min)}()
	end
	prob = NonlinearProblem{false}(func, kr0, p)
	sol = solve(prob, method)
	return reverse(sol.u)
end

"""
    find_spans(func, x, p)

    Find the spans where the function changes sign.

    # Arguments
    - `func::Function`: The function to find the spans for.
    - `x::Vector{Real}`: The x values.
    - `p::Vector{Real}`: The parameters.

    # Returns
    - `Vector{eltype(x)}`: The spans.
"""
function find_spans(func, x, p)
    idxes = [i for i in 1:(length(x) - 1) if (func(x[i], p) > 0 && func(x[i + 1], p) < 0)]
    if isempty(idxes)
        return Vector{eltype(p)}()
    end
    spans = [mean((x[i], x[i + 1])) for i in idxes]
    return spans
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
function solve_pekeris_equation(span, func_to_min, p; method = ITP())
	kr0 = mean(span)
	prob = NonlinearProblem(func_to_min, kr0, p)
	sol = solve(prob, TrustRegion())
	return sol.u
end

"""
    get_modal_function(env, krs, freq, zr, zs)

    Get the modal function for the Pekeris model.

    # Arguments
    - `env::PekerisUnderwaterEnv`: The underwater environment.
    - `krs::Vector{Real}`: The horizontal wavenumbers.
    - `freq::Real`: The frequency of the signal.
    - `zr::Real`: The receiver depth.
    - `zs::Real`: The source depth.

    # Returns
    - `Vector{Real}`: The modal function at the receiver depth.
    - `Vector{Real}`: The modal function at the source depth.
"""
function get_modal_function(env, krs, freq, zr, zs)
	@unpack c1, cb, ρ1, ρb, depth = env
	ω = 2π * freq
	kzw = sqrt.((ω / c1)^2 .- krs .^ 2) # vertical wavenumber in the water column
	kzb = sqrt.(krs .^ 2 .- (ω / cb)^2) # vertical wavenumber in the bottom half-space
	# The amplitude `amplitude` of the modes
	amplitude =
		sqrt.(
			(4 .* kzb .* kzw .* ρ1 .* ρb) ./
			(2 .* kzw .* ρ1 .* sin.(depth .* kzw) .^ 2 .- kzb .* ρb .* (-2 .* depth .* kzw .+ sin.(2 .* depth .* kzw)))
		)
	# The mode function
	function Ψ(z, mode)
		amplitude[mode] * sin(kzw[mode] * z) * (z .<= depth) +
		amplitude[mode] * sin(kzw[mode] * depth) * exp(-kzb[mode] * (z - depth)) * (z .> depth)
	end
	modes_zr = [Ψ(zr, mode) for mode in eachindex(krs)]
	modes_zs = [Ψ(zs, mode) for mode in eachindex(krs)]
	return modes_zr, modes_zs
end

"""
    pressure_f(env, krs, freq, r, zs, zr; t0 = 0.1, max_modes = Inf)

    Compute the pressure field for the Pekeris model.

    # Arguments
    - `env::PekerisUnderwaterEnv`: The underwater environment.
    - `krs::Vector{Real}`: The horizontal wavenumbers.
    - `freq::Real`: The frequency of the signal.
    - `r::Real`: The range.
    - `zs::Real`: The source depth.
    - `zr::Real`: The receiver depth.
    - `t0::Real=0.1`: The time offset.
    - `max_modes::Int=Inf`: The maximum number of modes to use.

    # Returns
    - `Complex{Real}`: The pressure field.
"""
function pressure_f(env::PekerisUnderwaterEnv, krs, freq, r, zs, zr; t0 = 0.1, max_modes = Inf)
	if freq == 0
		return 0.0 + 0.0im
	end
	if length(krs) == 0
		return 0.0 + 0.0im
	end
	@unpack c1, cb, ρ1, ρb, depth = env
	ψ_zr, ψ_zs = get_modal_function(env, krs, freq, zr, zs)
	t_offset = r / min(c1, cb) - t0  # align window correctly in time
	n_modes = min(max_modes, length(krs)) |> Int
	pf =
		-1im * exp(-1im * π / 4) / (ρ1 * sqrt(8π * r)) *
		sum([ψ_zs[mode] * ψ_zr[mode] * exp(-1im * krs[mode] * r) / sqrt(krs[mode]) for mode in range(1, n_modes)])
	return pf * exp(2im * pi * freq * t_offset)
end
