using NonlinearSolve
using UnPack
using StaticArrays
using IntervalArithmetic, IntervalArithmetic.Symbols, IntervalRootFinding
using Infiltrator

export PekerisUnderwaterEnv, find_kr, get_modal_function, pressure_f

struct PekerisUnderwaterEnv{T1 <: Real}
    c1::T1
    cb::T1
    ρ1::T1
    ρb::T1
    depth::T1

    function PekerisUnderwaterEnv(c1::T1, cb, ρ1, ρb, depth) where {T1}
        @assert cb>c1 "Bottom sound speed must be greater than water sound speed"
        new{T1}(c1, cb, ρ1, ρb, depth)
    end
end

function find_kr(env::PekerisUnderwaterEnv, freq, p; n_points = 5_000, method = NewtonRaphson())
    ω = 2π * freq
    kr_min, kr_max = extrema([ω / env.c1, ω / env.cb])
    # p = SA[env.c1, env.cb, env.ρ1, env.ρb, env.depth]
	func(kr, p) = @. tan(sqrt((ω / p[1])^2 - kr^2) * p[5]) +
                   (p[4] / p[3]) * sqrt((ω / p[1])^2 - kr^2) / sqrt(kr^2 - (ω / p[2])^2)

    kr_try = range(kr_min + eps(kr_min), kr_max - eps(kr_min); length = n_points)
    kr0 = find_spans(func, kr_try, p)
    prob = NonlinearProblem(func, kr0, p)
    sol = solve(prob, method)
    return reverse(sol.u)
end

function find_spans(func, x, p)
    idxes = [i for i in 1:(length(x) - 1) if (func(x[i], p) > 0 && func(x[i + 1], p) < 0)]
    if isempty(idxes)
        return Vector{eltype(x)}()
    end
    spans = [mean((x[i], x[i + 1])) for i in idxes]
    return spans
end

function solve_pekeris_equation(spans, func_to_min, p; method = ITP())
    prob = IntervalNonlinearProblem(func_to_min, spans, p)
    sol = solve(prob, method)
    return sol.u
end

function solve_pekeris_equation_v2(span, func_to_min, p; method = ITP())
    kr0 = mean(span)
    prob = NonlinearProblem(func_to_min, kr0, p)
    sol = solve(prob, TrustRegion())
    return sol.u
end

function get_modal_function(env, krs, freq, zr, zs)
    @unpack c1, cb, ρ1, ρb, depth = env
    ω = 2π * freq
    kzw = sqrt.((ω / c1)^2 .- krs .^ 2) # vertical wavenumber in the water column
    kzb = sqrt.(krs .^ 2 .- (ω / cb)^2) # vertical wavenumber in the bottom half-space
    # The amplitude `amplitude` of the modes
    amplitude = sqrt.(
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
    pf = -1im * exp(-1im * π / 4) / (ρ1 * sqrt(8π * r)) *
         sum([ψ_zs[mode] * ψ_zr[mode] * exp(-1im * krs[mode] * r) / sqrt(krs[mode]) for mode in range(1, n_modes)])
    return pf * exp(2im * pi * freq * t_offset)
end
