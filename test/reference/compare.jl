### Comparing the Julia solver against Fortran KRAKEN on the same environment.

using LinearAlgebra: dot, norm

# ForwardDiff is only needed for the opt-in group-speed comparison, but it is already a test
# dependency and the AD suite loads it anyway, so there is nothing to gain from deferring it.
using ForwardDiff

using Kraken: kraken_jl

"""
    FortranComparison

Result of [`compare_with_fortran`](@ref). `show` prints a per-mode table plus the two numbers that
usually decide whether a run is acceptable: the largest relative wavenumber difference and the
smallest mode-shape correlation.

Fields: `freq`, `n_julia`, `n_fortran`, `n_compared`, `kr_julia`, `kr_fortran`, `kr_absdiff`,
`kr_reldiff`, `mode_corr`, `group_speed_julia`, `group_speed_fortran`, `group_speed_reldiff`,
`depths`, `warnings`, `binary`.
"""
struct FortranComparison
    freq::Float64
    n_julia::Int
    n_fortran::Int
    n_compared::Int
    kr_julia::Vector{Float64}
    kr_fortran::Vector{Float64}
    kr_absdiff::Vector{Float64}
    kr_reldiff::Vector{Float64}
    mode_corr::Vector{Float64}
    group_speed_julia::Union{Nothing,Vector{Float64}}
    group_speed_fortran::Union{Nothing,Vector{Float64}}
    group_speed_reldiff::Union{Nothing,Vector{Float64}}
    depths::Vector{Float64}
    warnings::Vector{String}
    binary::String
end

max_kr_reldiff(c::FortranComparison) = isempty(c.kr_reldiff) ? NaN : maximum(c.kr_reldiff)
min_mode_corr(c::FortranComparison) = isempty(c.mode_corr) ? NaN : minimum(c.mode_corr)
function max_group_speed_reldiff(c::FortranComparison)
    return c.group_speed_reldiff === nothing || isempty(c.group_speed_reldiff) ? NaN : maximum(c.group_speed_reldiff)
end

function Base.show(io::IO, ::MIME"text/plain", c::FortranComparison)
    println(
        io,
        "FortranComparison at $(c.freq) Hz — $(c.n_julia) Julia modes vs $(c.n_fortran) Fortran " *
        "modes ($(basename(c.binary)))",
    )
    if c.n_julia != c.n_fortran
        println(io, "  ! mode-count mismatch: comparing the leading $(c.n_compared)")
    end
    println(io, "  mode        kr (Julia)      kr (Fortran)       |Δkr|      rel Δkr      corr")
    for i in 1:(c.n_compared)
        @printf(
            io,
            "  %4d  %16.10f  %16.10f  %10.3e  %10.3e  %8.6f\n",
            i,
            c.kr_julia[i],
            c.kr_fortran[i],
            c.kr_absdiff[i],
            c.kr_reldiff[i],
            c.mode_corr[i]
        )
    end
    @printf(io, "  max relative wavenumber difference : %.3e\n", max_kr_reldiff(c))
    @printf(io, "  min mode-shape correlation         : %.8f\n", min_mode_corr(c))
    if c.group_speed_reldiff !== nothing
        @printf(io, "  max relative group-speed difference: %.3e\n", max_group_speed_reldiff(c))
    end
    for w in c.warnings
        println(io, "  fortran warning: ", w)
    end
    return nothing
end

function Base.show(io::IO, c::FortranComparison)
    return print(io, "FortranComparison($(c.freq) Hz, $(c.n_julia) vs $(c.n_fortran) modes)")
end

# Piecewise-linear resampling. Deliberately hand-rolled rather than pulling DataInterpolations into
# the test environment: it is ten lines, and the alternative is a dependency that exists only to
# resample two vectors.
function _interp_linear(xs::AbstractVector, ys::AbstractVector, xq::Real)
    xq <= first(xs) && return float(first(ys))
    xq >= last(xs) && return float(last(ys))
    j = searchsortedlast(xs, xq)
    j >= length(xs) && return float(last(ys))
    x0, x1 = xs[j], xs[j + 1]
    x1 == x0 && return float(ys[j])
    t = (xq - x0) / (x1 - x0)
    return float(ys[j]) + t * (float(ys[j + 1]) - float(ys[j]))
end

"""
    _julia_mode_grid(sol) -> (depths, modes)

The Julia mode shapes with the surface point put back.

`get_z_vec` starts each layer's mesh at `z0 + Δz`, so the solver's grid *excludes* `z = 0` and the
returned `modes` matrix has no row for it. The Fortran `zTab` does include 0. Resampling without
that point would clamp to the first interior sample — a visibly wrong value right where the mode is
steepest — so the pressure-release condition `φ(0) = 0` is prepended, which is exact rather than an
approximation.
"""
function _julia_mode_grid(sol)
    z = Float64.(vcat(sol.props.zn_vec...))
    modes = Float64.(real.(sol.modes))
    if isempty(z) || first(z) > 0
        z = vcat(0.0, z)
        modes = vcat(zeros(1, size(modes, 2)), modes)
    end
    return z, modes
end

"""
    mode_correlation(zj, φj, zf, φf) -> Float64

Normalized inner product of two mode shapes sampled on different depth grids, in `[0, 1]`.

The two solvers do not share a mesh, so `φj` is resampled onto `zf` first. They also do not share a
*sign* or a normalization: KRAKEN flips each mode so it is positive at its near-surface turning
point, and normalizes against `∫φ²/ρ` including half-space terms, while Kraken.jl makes its own
choice. Taking `|⟨a,b⟩| / (‖a‖‖b‖)` removes both, leaving shape agreement alone — which is the thing
worth asserting.
"""
function mode_correlation(zj::AbstractVector, φj::AbstractVector, zf::AbstractVector, φf::AbstractVector)
    resampled = [_interp_linear(zj, φj, z) for z in zf]
    na, nb = norm(resampled), norm(φf)
    (na == 0 || nb == 0) && return 0.0
    return abs(dot(resampled, φf)) / (na * nb)
end

# The .mod holds every wavenumber but only in single precision; the .prt prints ten digits of the
# real part, for the subset of modes it lists. Overlaying the second on the first gives the most
# precise reference available for each mode. See MOD_WAVENUMBER_DIGITS.
function _best_fortran_kr(ref)
    kr = Float64.(real.(ref.kᵣ))
    if ref.grp !== nothing
        for (row, mode) in enumerate(ref.grp.m)
            1 <= mode <= length(kr) && (kr[mode] = real(ref.grp.kᵣ[row]))
        end
    end
    return kr
end

"""
    compare_with_fortran(env, freq; kwargs...) -> FortranComparison

Solve `env` at `freq` with both Kraken.jl and Fortran KRAKEN and report the differences.

# Keywords
- `nmodes` — cap on how many modes to compare. Defaults to every mode both solvers found.
- `atol`, `rtol` — root-finder tolerances forwarded to `kraken_jl`.
- `group_speeds` — also compare group speeds. **Off by default**, for two reasons: Kraken.jl gets
  them by differentiating `kraken_jl` with ForwardDiff, which costs seconds rather than
  milliseconds, and `AcousticsToolbox_jll`'s `kraken.exe` does not report them at all (see
  [`has_group_speeds`](@ref)). When the primary run turns out to have none, a second run with
  `complex=true` is made purely to obtain them.
- `complex`, `bindir`, `keep_files` and any `write_env_file` keyword are forwarded to
  [`run_fortran_kraken`](@ref).

# What is compared, and how

Mode counts are reported from each side separately — a mismatch is a finding, not an error, so the
comparison proceeds over the leading `min(n_julia, n_fortran)` modes and `show` flags it.

Wavenumbers are compared directly. Mode shapes cannot be: the two solvers use different depth
meshes, different sign conventions and different normalizations, so they are compared through
[`mode_correlation`](@ref).
"""
function compare_with_fortran(
    env, freq; nmodes=nothing, atol=1e-6, rtol=1e-6, group_speeds::Bool=false, complex::Bool=false, kwargs...
)
    freq = Float64(freq)
    sol = kraken_jl(env, freq; abstol=atol, reltol=rtol)
    ref = run_fortran_kraken(env, freq; complex=complex, kwargs...)

    n_julia = length(sol.kr)
    n_fortran = ref.nmodes
    n = min(n_julia, n_fortran)
    nmodes === nothing || (n = min(n, Int(nmodes)))

    kr_julia = Float64.(real.(sol.kr))[1:n]
    kr_fortran = _best_fortran_kr(ref)[1:n]
    kr_absdiff = abs.(kr_julia .- kr_fortran)
    kr_reldiff = kr_absdiff ./ abs.(kr_fortran)

    zj, φj = _julia_mode_grid(sol)
    mode_corr = [mode_correlation(zj, view(φj, :, m), ref.depths, real.(view(ref.ϕ, :, m))) for m in 1:n]

    vg_julia = nothing
    vg_fortran = nothing
    vg_reldiff = nothing
    if group_speeds
        grp = ref.grp
        binary = ref.binary
        if grp === nothing || !has_group_speeds(grp)
            # The jll's kraken.exe zeroes the column; krakenc.exe on the same input does not.
            second = run_fortran_kraken(env, freq; complex=true, kwargs...)
            grp = second.grp
            binary = second.binary
        end
        if grp !== nothing && has_group_speeds(grp)
            rows = [i for (i, mode) in enumerate(grp.m) if 1 <= mode <= n]
            modes = grp.m[rows]
            # Group speed is dω/dkᵣ; forward-mode AD straight through the solver gives it directly.
            dkr_dω = ForwardDiff.derivative(ω -> kraken_jl(env, ω / (2π); abstol=atol, reltol=rtol).kr, 2π * freq)
            vg_julia = [1 / dkr_dω[m] for m in modes]
            vg_fortran = grp.v[rows]
            vg_reldiff = abs.(vg_julia .- vg_fortran) ./ abs.(vg_fortran)
        end
    end

    return FortranComparison(
        freq,
        n_julia,
        n_fortran,
        n,
        kr_julia,
        kr_fortran,
        kr_absdiff,
        kr_reldiff,
        mode_corr,
        vg_julia,
        vg_fortran,
        vg_reldiff,
        ref.depths,
        ref.warnings,
        ref.binary,
    )
end
