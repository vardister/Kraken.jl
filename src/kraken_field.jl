### Pressure field and transmission loss
#
# Modes and wavenumbers are intermediates; this file turns them into the acoustic field. Everything
# here is a mode sum of the standard far-field form, so the only things worth pinning down precisely
# are the conventions — and all three of them are transcribed from the Acoustics Toolbox rather than
# from a textbook, because `field.exe` is what Milestone 7 validates against.
#
# **Phase convention.** KRAKEN works in `exp(i(ωt - kᵣr))`, so a mode contributes `exp(-i kᵣ r)`.
# `Evaluate` in `KrakenField/EvaluateMod.f90` builds
#
#     p(z, r) = i √(2π) e^{iπ/4} Σ_m φ_m(zₛ) φ_m(z) e^{-i kᵣ r} / √kᵣ / √r
#
# which is the *conjugate* of the `exp(-iωt)` form in Jensen et al. (5.13),
# `i e^{-iπ/4} / (ρ(zₛ)√(8πr)) Σ φφ e^{+ikᵣr}/√kᵣ`, times `4π`. We follow KRAKEN, sign and all, so
# [`acoustic_field`](@ref) can be compared with `field.exe` directly. Anything that expects the
# `exp(-iωt)` convention — `AcousticsToolbox.jl` conjugates for exactly this reason — needs `conj`.
#
# **The 4π.** It is not a slip: `Evaluate`'s header says "Normalized to pressure of point source at
# 1 meter", and `4π/√(8πr) = √(2π/r)` is precisely that normalization. So `-20log10|p|` is
# transmission loss in the usual sense, with no further scaling, which is also what the toolbox's own
# `plotshd` plots.
#
# **The 1/ρ(zₛ).** `Evaluate` has no density factor at all, because KRAKEN carries densities in
# g/cm³, where water is 1. Kraken.jl carries them in kg/m³, and its modes are normalized by
# `∫ φ²/ρ dz = 1` in those units, so they are √1000 larger than KRAKEN's. Dividing by ρ(zₛ) in kg/m³
# cancels both factors exactly and reproduces `field.exe` for a source in the water column — while
# also being the physically correct expression (Jensen's `1/ρ(zₛ)`) for a source in a sediment, where
# KRAKEN's missing factor would be wrong by ρ.
#
# **Depth interpolation** is linear, and receivers below the last mesh point ride the half-space
# exponential tail. That is `ReadOneMode` in `KrakenField/ReadModes.f90` line for line: `Weight`
# returns linear weights, and a receiver past `DepthB` gets `Phi(NTot)·exp(-γ_B(z - DepthB))`.
# Quadratic interpolation would be a defensible choice on its own, but not one that can be checked
# against anything.
#
# **Summation conventions: there are two, not three.** `field.exe` reads its task code from
# `Opt(4:4)` of the `.flp` file and accepts `'C'`, a blank, or `'I'`; anything else is a fatal
# `'Unknown option for coherent vs. incoherent mode addition'` (`KrakenField/field.f90`). The Matlab
# port agrees — `Matlab/Kraken/evalri.m` branches on `Opt(4:4) ~= 'I'` and nothing else.
# *Semi-coherent* is a BELLHOP run type, not a KRAKEN one: it multiplies an incoherent TL by the
# Lloyd-mirror source pattern `√2·sin(ω·zₛ·sinθ/c)` indexed by ray take-off angle
# (`Bellhop/influence.f90`), and a mode sum has no take-off angle — the surface image is already
# inside `φ_m(zₛ)`, so weighting by it again would double-count. So `mode` takes `:coherent` and
# `:incoherent`, and `:semicoherent` raises rather than inventing a convention nothing can check.
#
# The incoherent sum is `√(Σ_m |term_m|²)` — magnitudes, no phase. `Evaluate` writes it as
# `SQRT(SUM((Cmat*Hank)**2))` after forcing `ik = REAL(ik)` (which strips the oscillation and keeps
# the attenuation decay), and carries `1/√k`'s phase into the square; `evalri.m`, and `Evaluate`'s
# own commented-out line right above, take `ABS` of the whole term instead. The two agree exactly for
# a lossless solve and differ only in a phase for a lossy one, so we follow the `ABS` form, which is
# the one that is actually a magnitude. The result has no phase left to speak of, so the prefactor
# enters as `|Q|`.

export acoustic_field, transmission_loss, mode_amplitudes

"""
    tabulated_modes(sol::NormalModeSolution) -> (z, Φ)

The mode shapes together with the depths they are tabulated at, extended with the boundary samples
the mesh leaves out.

`sol.modes` is sampled on `reduce(vcat, sol.props.zn_vec)` and nowhere else, and that mesh omits a
boundary whose value the boundary condition already fixes: a pressure-release surface has no `z = 0`
sample and a vacuum bottom has no `z = D` sample, because in both cases `φ = 0` there and the
finite-difference system has no unknown to solve for. Interpolating without them would extrapolate
across the very region where the mode goes to zero, so they are put back here.

Returns the depth vector and a mode matrix with one column per mode, both with the same number of
rows.
"""
function tabulated_modes(sol::NormalModeSolution)
    z = reduce(vcat, sol.props.zn_vec)
    Φ = sol.modes
    T = eltype(Φ)
    if sol.env.top_bc isa PressureRelease
        z = vcat(zero(eltype(z)), z)
        Φ = vcat(zeros(T, 1, size(Φ, 2)), Φ)
    end
    if sol.env.bottom_bc isa PressureRelease
        z = vcat(z, sol.env.depth)
        Φ = vcat(Φ, zeros(T, 1, size(Φ, 2)))
    end
    return z, Φ
end

"""
    halfspace_decay(sol::NormalModeSolution, m)

The vertical decay rate `γ = √(kᵣ² - (ω/c_b)²)` of mode `m` in the bottom half-space.

Only meaningful over an [`AcousticHalfspace`](@ref) bottom; a rigid or vacuum bottom has no medium
below `D` for a receiver to sit in.
"""
function halfspace_decay(sol::NormalModeSolution, m)
    kb2 = (2π * sol.props.freq / sol.env.cb)^2
    return sqrt(sol.kr[m]^2 - kb2)
end

# Linear interpolation of one mode, with the half-space tail below the deepest sample. `zt` is
# sorted and starts at 0 (the only two supported surfaces both put a sample there, one of them via
# `tabulated_modes`), so a receiver above the first sample can only mean a negative depth.
function mode_amplitude(zt, Φ, m, z, sol)
    z < zero(z) &&
        throw(ArgumentError("receiver depth $z is above the surface; depths are measured positive downwards"))
    if z > zt[end]
        sol.env.bottom_bc isa AcousticHalfspace || throw(
            ArgumentError(
                "receiver depth $z is below the bottom at $(zt[end]) m, and a $(sol.env.bottom_bc) " *
                "bottom has no medium below it to propagate into",
            ),
        )
        return Φ[end, m] * exp(-halfspace_decay(sol, m) * (z - zt[end]))
    end
    j = min(searchsortedlast(zt, z), length(zt) - 1)
    j = max(j, 1)
    w = (z - zt[j]) / (zt[j + 1] - zt[j])
    return Φ[j, m] + w * (Φ[j + 1, m] - Φ[j, m])
end

"""
    mode_amplitudes(sol::NormalModeSolution, z; nmodes=length(sol.kr))

Evaluate the mode shapes of `sol` at the depths `z` (metres, positive downwards).

Returns a `length(z) × nmodes` matrix — row per depth, column per mode — for a vector `z`, and a
`1 × nmodes` matrix for a scalar one. Interpolation is linear between mesh points; below the last
mesh point the mode continues as `φ(D)·exp(-γ(z - D))` into the bottom half-space, which is what
KRAKEN's `ReadOneMode` does. Depths below a rigid or vacuum bottom are an error, as is a negative
depth.

# Example
```julia
sol = kraken_jl(UnderwaterEnv(pekeris_env()...), 100.0)
mode_amplitudes(sol, 0:10:100)     # 11 × M
```
"""
function mode_amplitudes(sol::NormalModeSolution, z; nmodes=length(sol.kr))
    depths = as_vector(z)
    M = min(nmodes, length(sol.kr))
    M <= 0 && return zeros(eltype(sol.modes), length(depths), 0)
    zt, Φ = tabulated_modes(sol)
    return [mode_amplitude(zt, Φ, m, depths[i], sol) for i in eachindex(depths), m in 1:M]
end

as_vector(x::Real) = [x]
as_vector(x) = collect(x)

"""
    SUMMATION_MODES

The mode-summation conventions [`acoustic_field`](@ref) accepts: `:coherent` and `:incoherent`,
matching `field.exe`'s `.flp` task codes `C` (or blank) and `I`. There is no third one — see the
note at the top of `kraken_field.jl` for why `:semicoherent` is not among them.
"""
const SUMMATION_MODES = (:coherent, :incoherent)

# Resolved once per call, outside the field loop, so the summation stays type stable.
function summation_mode(mode::Symbol)
    mode in SUMMATION_MODES && return Val(mode)
    mode === :semicoherent && throw(
        ArgumentError(
            "`mode = :semicoherent` has no normal-mode definition. It is a BELLHOP run type: an " *
            "incoherent TL weighted by the Lloyd-mirror source pattern, which is indexed by ray " *
            "take-off angle. A mode sum has no take-off angle, and the surface image is already " *
            "inside φ(zₛ), so the weight would double-count. `field.exe` has no such option " *
            "either — it accepts only C/blank and I. Use :coherent or :incoherent.",
        ),
    )
    return throw(ArgumentError("`mode = :$mode` is not a summation mode; expected one of $SUMMATION_MODES"))
end

# The mode sum proper. `terms` is a generator of the per-mode contributions, consumed once.
mode_sum(::Val{:coherent}, Q, terms) = Q * sum(terms)
mode_sum(::Val{:incoherent}, Q, terms) = complex(abs(Q) * sqrt(sum(abs2, terms)))

"""
    acoustic_field(sol, ranges, zs, zr; nmodes=length(sol.kr), mode=:coherent)

The complex acoustic pressure field of the normal-mode solution `sol` for a point source at depth
`zs`, on the grid of `ranges` (metres) and receiver depths `zr` (metres).

Returns a `length(zr) × length(ranges)` matrix — depth down the rows, range across the columns —
always, including when `ranges` or `zr` is a single number. `transmission_loss` is the same field in
dB.

# Keyword arguments
- `nmodes`: how many modes to sum, from the first (largest `kᵣ`). Defaults to all of them; the
  solution has already discarded everything that is not propagating.
- `mode`: `:coherent` or `:incoherent`, matching `field.exe`'s `.flp` task codes `C` (or blank) and
  `I`. A coherent sum adds the modes as complex amplitudes and so shows the interference structure —
  the nulls a real receiver sees at one frequency. An incoherent sum adds their intensities,
  `√(Σ|term|²)`, which discards the phase and leaves the smooth mean level: what you want when the
  interference pattern is unresolved, unknown, or about to be averaged away anyway. `:semicoherent`
  is deliberately absent; ask for it and the error says why.

# Conventions

The phase convention is KRAKEN's `exp(i(ωt - kᵣr))`, and the amplitude is normalized to the pressure
of a point source at 1 m, so `field.exe`'s output can be compared without rescaling. Take `conj` for
the `exp(-iωt)` convention used in Jensen et al. See the comments at the top of `kraken_field.jl` for
where each factor comes from.

An incoherent field has thrown its phase away by construction, so what comes back is a magnitude
carried in a complex number for a uniform return type — real, non-negative, and not a pressure you
can interfere with anything. Only `:coherent` returns a genuine complex pressure.

The field is the far-field mode sum, which is singular at `r = 0` and inaccurate within a wavelength
or so of the source; zero and negative ranges are rejected rather than returned as `Inf`.

# Example
```julia
sol = kraken_jl(UnderwaterEnv(pekeris_env()...), 100.0)
p = acoustic_field(sol, 1_000:100:10_000, 25.0, 0:5:100)
```
"""
function acoustic_field(sol::NormalModeSolution, ranges, zs, zr; nmodes=length(sol.kr), mode=:coherent)
    conv = summation_mode(mode)
    r = as_vector(ranges)
    zrv = as_vector(zr)
    any(<=(0), r) && throw(ArgumentError("ranges must be positive; the mode sum is singular at r = 0"))

    M = min(nmodes, length(sol.kr))
    T = promote_type(eltype(sol.kr), eltype(sol.modes), eltype(r), typeof(zs), eltype(zrv), ComplexF64)
    M <= 0 && return zeros(T, length(zrv), length(r))

    krs = sol.kr[1:M]
    ϕs = vec(mode_amplitudes(sol, zs; nmodes=M))
    ϕr = mode_amplitudes(sol, zrv; nmodes=M)
    # `Evaluate`'s `factor`, with the density of the source medium folded in — see the file header.
    Q = im * sqrt(2π) * exp(im * π / 4) / density(sol.env.ρ, zs)

    return [
        T(mode_sum(conv, Q, (ϕs[m] * ϕr[i, m] * exp(-im * krs[m] * r[j]) / sqrt(krs[m]) for m in 1:M)) / sqrt(r[j])) for
        i in eachindex(zrv), j in eachindex(r)
    ]
end

"""
    transmission_loss(sol, ranges, zs, zr; nmodes=length(sol.kr), mode=:coherent)

Transmission loss in dB, `-20log10|p|`, for the field described by [`acoustic_field`](@ref) — same
arguments, same `length(zr) × length(ranges)` layout, same conventions.

The reference is the pressure a point source of the same strength would produce at 1 m in free
space, which is the normalization built into the field itself, so nothing further is subtracted
here. Larger numbers mean a weaker signal.

# Example
```julia
sol = kraken_jl(UnderwaterEnv(pekeris_env()...), 100.0)
tl = transmission_loss(sol, 1_000:100:10_000, 25.0, 50.0)
smooth = transmission_loss(sol, 1_000:100:10_000, 25.0, 50.0; mode=:incoherent)
```
"""
function transmission_loss(sol::NormalModeSolution, ranges, zs, zr; kwargs...)
    return -20 .* log10.(abs.(acoustic_field(sol, ranges, zs, zr; kwargs...)))
end

### Broadband synthesis
#
# A pulse is the inverse Fourier transform of the field over frequency. Two things decide how this
# is written.
#
# **No FFT dependency.** The transform is the Riemann sum of the continuous inverse transform over
# the positive frequencies. In KRAKEN's `exp(i(ωt - kᵣr))` convention a mode carries `e^{+i2πft}` in
# time, and for a real signal the negative frequencies are the conjugates of the positive ones, so
#
#     p(t) = 2 Δf Re[ Σ_k S(f_k) P(f_k) e^{i2πf_k t} ]
#
# is the whole of it. That costs `O(N_f · N_t)` against an FFT's `O(N log N)`, which would matter if
# it were the expensive part — but every one of those `N_f` coefficients costs a full `kraken_jl`
# eigen-solve to produce, and that dominates by orders of magnitude. Writing the sum out buys
# freedom from FFTW (a binary artifact this package pruned in 0.3.0 and would now be re-adding for
# one line), a transform that is differentiable with no extra rules for task 7.5, and a formula that
# states its own phase convention instead of inheriting one from a library.
#
# **The reduction time is the caller's.** The staged `dev/` code folded `exp(2πif·t₀)` with
# `t₀ = r/c_max - t0_offset` into the field itself. That is a windowing convenience — it slides the
# arrival to a known place in the time axis so a short window can hold it — not physics, and baking
# it into the field is what made the staged code's two methods disagree about what they returned. It
# lives here as `t_reduce` on the synthesis, where it is visible.

export broadband_field, synthesize_pulse

"""
    broadband_field(env, freqs, ranges, zs, zr; nmodes=typemax(Int), mode=:coherent, kwargs...)

Solve `env` at every frequency in `freqs` and evaluate the field at each, for pulse synthesis.

Returns a `length(freqs) × length(zr) × length(ranges)` array of complex pressure — frequency down
the first axis, then the same depth × range layout [`acoustic_field`](@ref) returns. Extra `kwargs`
go to [`kraken_jl`](@ref).

A frequency that supports no trapped modes contributes zeros rather than an error, and `f = 0` is
zero by definition; both are normal at the low end of a band, and both are what the transform wants
there anyway.

Each frequency is an independent solve, which is the whole reason this is worth having as its own
function: the loop below is a `map` over `freqs` precisely so that turning it into a threaded map is
the entire change needed to parallelize it. That is deliberately not done yet — a parallel version
of a wrong answer is still wrong, and task 7.4 is what confirms the serial one.

# Example
```julia
env = UnderwaterEnv(pekeris_env()...)
P = broadband_field(env, 10.0:1.0:200.0, 5_000.0, 36.0, 55.0)   # 191 × 1 × 1
```
"""
function broadband_field(env::UnderwaterEnv, freqs, ranges, zs, zr; nmodes=typemax(Int), mode=:coherent, kwargs...)
    r = as_vector(ranges)
    zrv = as_vector(zr)
    conv = summation_mode(mode)   # reject a bad `mode` once, not once per frequency

    # `map`, not a preallocated buffer written in a loop: the point of this function is that the
    # frequencies are independent, and `map` is the form that says so (and that `ThreadsX.map` or an
    # `@threads` fill can replace wholesale when 7.4 has signed off on the serial answer).
    planes = map(freqs) do f
        f <= 0 && return zeros(ComplexF64, length(zrv), length(r))
        sol = kraken_jl(env, f; kwargs...)
        isempty(sol.kr) && return zeros(ComplexF64, length(zrv), length(r))
        return ComplexF64.(acoustic_field(sol, r, zs, zrv; nmodes=nmodes, mode=conv_symbol(conv)))
    end

    P = Array{ComplexF64,3}(undef, length(planes), length(zrv), length(r))
    for (k, plane) in enumerate(planes)
        P[k, :, :] = plane
    end
    return P
end

conv_symbol(::Val{S}) where {S} = S

"""
    synthesize_pulse(P, freqs, times; spectrum=nothing, t_reduce=0)

Transform a complex spectrum `P` sampled at `freqs` into a real time series at `times`.

`P` is indexed by frequency along its first axis — a vector for one receiver, or the array
[`broadband_field`](@ref) returns — and the result replaces that axis with one of length
`length(times)`. `freqs` must be evenly spaced; the spacing is the `Δf` of the transform.

# Keyword arguments
- `spectrum`: the source spectrum, as a vector over `freqs` or a function of frequency. The default
  is flat, which makes the result the band-limited impulse response of the waveguide.
- `t_reduce`: subtracted from the travel time, so `times` are measured from `t_reduce` rather than
  from the source instant. Passing `r / maximum(env.c.c)` puts the earliest possible arrival at
  `t = 0` and lets a short window hold a pulse that travelled a long way.

The phase convention is KRAKEN's, the same one [`acoustic_field`](@ref) documents: a mode carries
`e^{i(ωt - kᵣr)}`, so this sums `e^{+i2πft}` and takes twice the real part, the positive-frequency
half of a real signal's inverse transform.

# Example
```julia
env = UnderwaterEnv(pekeris_env()...)
freqs = 2.0:2.0:300.0
P = broadband_field(env, freqs, 10_000.0, 36.0, 55.0)
t = 0:0.002:2.0
p = synthesize_pulse(P[:, 1, 1], freqs, t; t_reduce=10_000.0 / maximum(env.c.c))
```
"""
function synthesize_pulse(P::AbstractVector, freqs, times; spectrum=nothing, t_reduce=0)
    f = as_vector(freqs)
    length(f) == length(P) ||
        throw(DimensionMismatch("`P` has $(length(P)) frequency samples but `freqs` has $(length(f))"))
    length(f) >= 2 || throw(ArgumentError("need at least two frequencies to set the transform's Δf"))
    Δf = f[2] - f[1]
    all(k -> isapprox(f[k + 1] - f[k], Δf; rtol=1e-8), 1:(length(f) - 1)) ||
        throw(ArgumentError("`freqs` must be evenly spaced; the spacing is the transform's Δf"))

    S = source_spectrum(spectrum, f)
    return [2 * Δf * real(sum(S[k] * P[k] * cis(2π * f[k] * (t + t_reduce)) for k in eachindex(f))) for t in times]
end

function synthesize_pulse(P::AbstractArray{<:Any,3}, freqs, times; kwargs...)
    nt = length(as_vector(times))
    out = Array{Float64,3}(undef, nt, size(P, 2), size(P, 3))
    for j in axes(P, 3), i in axes(P, 2)
        out[:, i, j] = synthesize_pulse(view(P, :, i, j), freqs, times; kwargs...)
    end
    return out
end

source_spectrum(::Nothing, f) = ones(length(f))
function source_spectrum(s::AbstractVector, f)
    return if length(s) == length(f)
        s
    else
        throw(DimensionMismatch("`spectrum` has $(length(s)) entries but `freqs` has $(length(f))"))
    end
end
source_spectrum(s, f) = [s(fk) for fk in f]
