### Comparing the Julia solver against Fortran KRAKEN on the same environment.

using LinearAlgebra: dot, norm

# ForwardDiff is only needed for the opt-in group-speed comparison, but it is already a test
# dependency and the AD suite loads it anyway, so there is nothing to gain from deferring it.
using ForwardDiff

using Kraken: kraken_jl, PressureRelease, transmission_loss

"""
    FortranComparison

Result of [`compare_with_fortran`](@ref). `show` prints a per-mode table plus the two numbers that
usually decide whether a run is acceptable: the largest relative wavenumber difference and the
smallest mode-shape correlation.

Fields: `freq`, `n_julia`, `n_fortran`, `n_compared`, `kr_julia`, `kr_fortran`, `kr_absdiff`,
`kr_reldiff`, `alpha_julia`, `alpha_fortran`, `alpha_absdiff`, `alpha_reldiff`, `mode_corr`,
`group_speed_julia`, `group_speed_fortran`, `group_speed_reldiff`, `depths`, `warnings`, `binary`,
`group_speed_detail`.

`group_speed_detail` is, when group speeds were compared, one named tuple per compared row —
`(julia_mode, fortran_mode, kr_julia, kr_fortran, vg_julia, vg_fortran, reldiff)` — plus the path of the
binary that printed them in `group_speed_binary`. It exists so a group-speed disagreement can be traced
to a mode from a CI log alone, where the Fortran binary cannot be rerun by hand.

`alpha_*` are the **imaginary** parts of the wavenumbers — the modal attenuation. They are zero on
both sides for a lossless environment; see [`max_alpha_reldiff`](@ref) for why they get a tolerance
of their own rather than sharing the wavenumber one.
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
    alpha_julia::Vector{Float64}
    alpha_fortran::Vector{Float64}
    alpha_absdiff::Vector{Float64}
    alpha_reldiff::Vector{Float64}
    mode_corr::Vector{Float64}
    group_speed_julia::Union{Nothing,Vector{Float64}}
    group_speed_fortran::Union{Nothing,Vector{Float64}}
    group_speed_reldiff::Union{Nothing,Vector{Float64}}
    depths::Vector{Float64}
    warnings::Vector{String}
    binary::String
    group_speed_detail::Union{Nothing,Vector{NamedTuple}}
    group_speed_binary::Union{Nothing,String}
end

max_kr_reldiff(c::FortranComparison) = isempty(c.kr_reldiff) ? NaN : maximum(c.kr_reldiff)
min_mode_corr(c::FortranComparison) = isempty(c.mode_corr) ? NaN : minimum(c.mode_corr)

"""
    max_alpha_reldiff(c) -> Float64
    max_alpha_absdiff(c) -> Float64

Worst disagreement on the modal attenuation `Im(kᵣ)` across the compared modes. `NaN` when nothing
was compared; `0.0` when both solvers report a lossless waveguide.

**These warrant a looser tolerance than the wavenumbers, for three compounding reasons**, which is
why they are reported separately rather than folded into `max_kr_reldiff`:

  * `Im(kᵣ)` is three to five orders of magnitude smaller than `Re(kᵣ)`, so it is the *last* digits
    of the complex number and inherits the whole absolute error of the solve.
  * The only usable Fortran source for it is the `.mod` file, which is single precision — about 7
    significant digits (see [`MOD_WAVENUMBER_DIGITS`](@ref)). The `.prt` prints it with `G10.2`, i.e.
    **two**, which is useless here.
  * Both solvers compute it by first-order perturbation on a mode shape taken from their own
    coarsest mesh, and those meshes differ. The wavenumbers get Richardson extrapolation; the
    attenuation does not.

`alpha_reldiff` is relative where the Fortran value is nonzero and falls back to the absolute
difference where it is zero, so a lossless case reports `0.0` rather than `NaN`.
"""
max_alpha_reldiff(c::FortranComparison) = isempty(c.alpha_reldiff) ? NaN : maximum(c.alpha_reldiff)
max_alpha_absdiff(c::FortranComparison) = isempty(c.alpha_absdiff) ? NaN : maximum(c.alpha_absdiff)
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
    if any(!iszero, c.alpha_fortran) || any(!iszero, c.alpha_julia)
        println(io, "  mode     Im kr (Julia)    Im kr (Fortran)      |Δα|      rel Δα")
        for i in 1:(c.n_compared)
            @printf(
                io,
                "  %4d  %16.9e  %16.9e  %10.3e  %10.3e\n",
                i,
                c.alpha_julia[i],
                c.alpha_fortran[i],
                c.alpha_absdiff[i],
                c.alpha_reldiff[i]
            )
        end
        @printf(io, "  max relative attenuation difference: %.3e\n", max_alpha_reldiff(c))
    end
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

The Julia mode shapes with the pressure-release end points put back.

A pressure-release boundary is not a mesh point, because `ψ` is already known there (plan task 6.2).
So a vacuum surface's grid starts at `Δz` and a vacuum bottom's stops at `D − Δz`, while the Fortran
`zTab` runs from 0 to `D`. Resampling without those points would clamp to the nearest interior
sample, which is a visibly wrong value right where the mode is steepest. So `φ = 0` is added at each
missing end, which is exact rather than an approximation. A rigid or half-space end is already a
mesh point and gets nothing.
"""
function _julia_mode_grid(sol)
    z = Float64.(vcat(sol.props.zn_vec...))
    modes = Float64.(real.(sol.modes))
    if isempty(z) || first(z) > 0
        z = vcat(0.0, z)
        modes = vcat(zeros(1, size(modes, 2)), modes)
    end
    depth = Float64(sol.env.depth)
    if sol.env.bottom_bc isa PressureRelease && last(z) < depth
        z = vcat(z, depth)
        modes = vcat(modes, zeros(1, size(modes, 2)))
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
  milliseconds, and `AcousticsToolbox_jll`'s `kraken.exe` cannot be trusted for them (see
  [`has_group_speeds`](@ref)). They always come from `krakenc.exe`: unless the primary run already
  was one (`complex=true`), a second, complex run is made purely to obtain them.
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

    # The modal attenuation. Unlike `Re(kᵣ)` this comes from the `.mod` and only the `.mod` -- the
    # `.prt`'s two printed digits cannot resolve it. See `max_alpha_reldiff`.
    alpha_julia = Float64.(imag.(sol.kr))[1:n]
    alpha_fortran = Float64.(imag.(ref.kᵣ))[1:n]
    alpha_absdiff = abs.(alpha_julia .- alpha_fortran)
    alpha_reldiff = [
        iszero(alpha_fortran[i]) ? alpha_absdiff[i] : alpha_absdiff[i] / abs(alpha_fortran[i]) for i in 1:n
    ]

    zj, φj = _julia_mode_grid(sol)
    mode_corr = [mode_correlation(zj, view(φj, :, m), ref.depths, real.(view(ref.ϕ, :, m))) for m in 1:n]

    vg_julia = nothing
    vg_fortran = nothing
    vg_reldiff = nothing
    vg_detail = nothing
    vg_binary = nothing
    if group_speeds
        grp = ref.grp
        binary = ref.binary
        if !complex
            # Group speeds always come from krakenc.exe, never from AcousticsToolbox_jll v2025.9's
            # kraken.exe, which gets them wrong differently on each platform: the macOS build prints
            # 0.00000 for every mode, and the Linux build prints plausible values that are not all
            # right -- Munk at 100 Hz, mode 199, wavenumber matching krakenc's and Kraken.jl's to 8
            # digits: krakenc 1452.59 m/s, Kraken.jl 1452.47, Linux kraken.exe 1500.15 (CI,
            # 2026-10-13). Falling back only when the column is all zeros let the Linux values
            # through and made the comparison platform-dependent.
            second = run_fortran_kraken(env, freq; complex=true, kwargs...)
            grp = second.grp
            binary = second.binary
        end
        if grp !== nothing && has_group_speeds(grp)
            # Pair each printed row with a Julia mode by *wavenumber*, not by mode number. The rows
            # usually come from a second, krakenc.exe run, a different program from the one `n`
            # counts, and nothing guarantees it numbers the modes the same way: near cutoff it can
            # find a mode the others do not, which shifts every number after it. Index pairing
            # then compares neighbouring modes. A row only counts if its kᵣ matches the Julia
            # mode's to 1e-4, an order looser than any wavenumber tolerance asserted on these
            # environments and orders tighter than the gap between adjacent modes.
            kr_all = Float64.(real.(sol.kr))
            matched = map(eachindex(grp.m)) do i
                k = real(grp.kᵣ[i])
                m = argmin(abs.(kr_all .- k))
                return abs(kr_all[m] - k) <= 1e-4 * abs(k) ? m : 0
            end
            rows = [i for i in eachindex(matched) if matched[i] > 0]
            modes = matched[rows]
            # Group speed is dω/dkᵣ; forward-mode AD straight through the solver gives it directly.
            dkr_dω = ForwardDiff.derivative(ω -> kraken_jl(env, ω / (2π); abstol=atol, reltol=rtol).kr, 2π * freq)
            vg_julia = [1 / dkr_dω[m] for m in modes]
            vg_fortran = grp.v[rows]
            vg_reldiff = abs.(vg_julia .- vg_fortran) ./ abs.(vg_fortran)
            vg_detail = NamedTuple[
                (
                    julia_mode=modes[j],
                    fortran_mode=grp.m[rows[j]],
                    kr_julia=kr_all[modes[j]],
                    kr_fortran=real(grp.kᵣ[rows[j]]),
                    vg_julia=vg_julia[j],
                    vg_fortran=vg_fortran[j],
                    reldiff=vg_reldiff[j],
                ) for j in eachindex(modes)
            ]
            vg_binary = binary
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
        alpha_julia,
        alpha_fortran,
        alpha_absdiff,
        alpha_reldiff,
        mode_corr,
        vg_julia,
        vg_fortran,
        vg_reldiff,
        ref.depths,
        ref.warnings,
        ref.binary,
        vg_detail,
        vg_binary,
    )
end

### Comparing the pressure field against `field.exe`
#
# What can be compared and what cannot. `Kraken.acoustic_field` was built in plan task 7.1 to
# reproduce `Evaluate`'s own expression — the `i·√(2π)·e^{iπ/4}` prefactor, the 1 m reference that
# supplies the 4π, and KRAKEN's `e^{i(ωt-kᵣr)}` phase convention — so the *complex* fields should
# agree, not merely their magnitudes. The one factor Kraken.jl carries and KRAKEN does not is
# `1/ρ(zₛ)`, which cancels exactly for a source in water of 1 g/cm³ and is the physically correct
# generalization otherwise (see the header of `src/kraken_field.jl`).
#
# Two things are deliberately *not* compared:
#
#   * **The range r = 0**, where the mode sum is singular. `Evaluate` guards it with a `WHERE` and
#     leaves the pressure unscaled; Kraken.jl rejects it outright. Neither answer means anything.
#   * **Depths at a pressure-release boundary**, where both fields are zero and the dB difference
#     is `Inf - Inf`. Those are dropped rather than given a fudge factor.

"""
    FieldComparison

Result of [`compare_field_with_fortran`](@ref): the two transmission-loss fields on a common grid,
with the summary statistics plan task 7.4 asks for.

`max_db` and `rms_db` are taken over the compared points only — see `n_compared` and `n_dropped`.
"""
struct FieldComparison
    freq::Float64
    mode::Symbol
    ranges::Vector{Float64}
    zr::Vector{Float64}
    zs::Float64
    tl_julia::Matrix{Float64}
    tl_fortran::Matrix{Float64}
    absdiff::Matrix{Float64}
    max_db::Float64
    rms_db::Float64
    n_compared::Int
    n_dropped::Int
    warnings::Vector{String}
end

function Base.show(io::IO, ::MIME"text/plain", c::FieldComparison)
    println(io, "FieldComparison at ", c.freq, " Hz, :", c.mode)
    println(io, "  grid       : ", length(c.zr), " depths × ", length(c.ranges), " ranges")
    println(io, "  compared   : ", c.n_compared, " points (", c.n_dropped, " dropped)")
    println(io, "  max |ΔTL|  : ", round(c.max_db; digits=4), " dB")
    print(io, "  RMS ΔTL    : ", round(c.rms_db; digits=4), " dB")
    isempty(c.warnings) || print(io, "\n  warnings   : ", length(c.warnings))
    return nothing
end

"""
    compare_field_with_fortran(env, freq, ranges, zs, zr; mode=:coherent, nmodes=9999, kwargs...)

Run both field solvers over the same range-depth grid and return a printable [`FieldComparison`](@ref).

`ranges` are metres, `zs` and `zr` metres. `nmodes` caps both solvers at the same mode count so a
mode-count disagreement does not masquerade as a field disagreement. Extra keyword arguments reach
[`write_env_file`](@ref) and, through it, the `.env` the Fortran side reads.
"""
function compare_field_with_fortran(env, freq, ranges, zs, zr; mode::Symbol=:coherent, nmodes::Integer=9999, kwargs...)
    rs = Float64.(collect(ranges))
    zrv = Float64.(collect(zr))
    ref = run_fortran_field(env, freq, rs, zs, zrv; mode=mode, nmodes=nmodes, kwargs...)

    # field.exe sorts and may subtabulate the grid it was given, so take the grid back from the
    # .shd rather than assuming it survived unchanged.
    sol = kraken_jl(env, freq)
    tl_julia = transmission_loss(sol, ref.ranges, zs, ref.zr; nmodes=nmodes, mode=mode)
    tl_fortran = ref.tl

    size(tl_julia) == size(tl_fortran) ||
        error("grid mismatch: Kraken.jl gave $(size(tl_julia)), field.exe gave $(size(tl_fortran))")

    finite = isfinite.(tl_julia) .& isfinite.(tl_fortran)
    absdiff = abs.(tl_julia .- tl_fortran)
    compared = absdiff[finite]
    return FieldComparison(
        Float64(freq),
        mode,
        collect(ref.ranges),
        collect(ref.zr),
        Float64(zs),
        tl_julia,
        tl_fortran,
        absdiff,
        isempty(compared) ? NaN : maximum(compared),
        isempty(compared) ? NaN : sqrt(sum(abs2, compared) / length(compared)),
        length(compared),
        length(absdiff) - length(compared),
        ref.warnings,
    )
end
