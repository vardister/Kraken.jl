### Serialization of a Kraken.jl environment into an Acoustics Toolbox `.env` file.
#
# The grammar implemented here was read off the Fortran that consumes it — `ReadEnvironment`,
# `ReadTopOpt`, `TopBot` (`misc/ReadEnvironmentMod.f90`), `ReadSSP` (`misc/sspMod.f90`),
# `ReadVector`/`ReadfreqVec` (`misc/SourceReceiverPositions.f90`) and `SubTab`
# (`misc/subtabulate.f90`) — cross-checked against the files in `test/standard_envs/`.
#
# Record order (single-frequency, the case this writer emits):
#
#     'TITLE'
#     FREQ
#     NMEDIA
#     'TOPOPT'
#     NMESH  SIGMA  Z(NSSP)          \  once per medium
#     Z CP CS RHO AP AS  /           /  NSSP rows, last one at Z(NSSP)
#     'BOTOPT'  SIGMA
#     ZB CPB CSB RHOB APB ASB /
#     CLOW  CHIGH
#     RMAX
#     NSD
#     SD(1) ... /
#     NRD
#     RD(1) ... /
#
# With a frequency *vector* two extra records follow (`NFREQ`, then the frequencies) and `'B'` is
# placed in column 6 of the top-option string — that is what `ReadfreqVec` keys off.

using Printf

using Kraken: UnderwaterEnv, UnderwaterEnvFORTRAN

"""
Densities in a KRAKEN `.env` file are in g/cm³, while Kraken.jl's standard environments carry them
in kg/m³ (`pekeris_env` uses `ρ0 = 1000.0`, and the checked-in `Pekeris_AV.env` writes `1.0` for the
same layer). This is the factor applied on the way out; override it with the `density_scale` keyword
if your environment is already in g/cm³.
"""
const DEFAULT_DENSITY_SCALE = 1e-3

"""
Top-option string written by default.

- `C` — C-linear SSP interpolation. **Deliberately not `S`** (cubic spline), which is what the
  checked-in `.env` files use: `Kraken.jl` interpolates its sound-speed profile linearly
  (`SampledSSP` wraps `DataInterpolations.LinearInterpolation`), so `C` is the option that makes the
  Fortran run solve the *same* problem. The two agree for any medium given by two points, which is
  every standard environment except `munk_env`. Milestone 6 adds the other interpolators.
- `V` — vacuum (pressure-release) above the top interface, which is the only surface condition the
  finite-difference scheme currently implements.
- `W` — attenuation in dB/wavelength. Irrelevant while every αp/αs is zero, but it has to be a
  *valid* letter or `ReadTopOpt` calls `ERROUT`.
"""
const DEFAULT_TOPOPT = "CVW"

# ---------------------------------------------------------------------------------------------
# Number formatting
# ---------------------------------------------------------------------------------------------

# Fortran list-directed input is happy with either fixed or exponent form; fixed reads better and
# diffs better against the checked-in files, so it is used wherever it does not lose precision.
function _fmt(x::Real)
    v = Float64(x)
    isfinite(v) || error("Cannot write non-finite value $v to a .env file")
    v == 0 && return "0.0"
    (abs(v) < 1e-3 || abs(v) >= 1e7) && return @sprintf("%.8e", v)
    s = @sprintf("%.6f", v)
    s = rstrip(s, '0')
    endswith(s, '.') && (s *= "0")
    return s
end

# `Title` is read into columns 9:80 of an 80-character buffer, so 72 characters survive. Quotes are
# doubled because the record is read list-directed as a single-quoted string.
function _fmt_title(title::AbstractString)
    t = replace(String(title), '\n' => ' ', '\r' => ' ')
    t = replace(t, '\'' => "''")
    length(t) > 72 && (t = t[1:72])
    return "'" * t * "'"
end

# ---------------------------------------------------------------------------------------------
# Depth / frequency vectors
# ---------------------------------------------------------------------------------------------

# `ReadVector` presets x(2) and x(3) to -999.9, reads list-directed, then calls `SubTab`, which
# expands the array into `n` uniform points on [x(1), x(2)] when x(3) survived untouched. So a
# uniformly spaced vector of three or more points can be written as its two endpoints plus `/`,
# exactly as the checked-in files do. Anything else is written out in full.
function _is_uniform(v::AbstractVector{<:Real})
    length(v) < 3 && return false
    lo, hi = Float64(first(v)), Float64(last(v))
    δ = (hi - lo) / (length(v) - 1)
    tol = 1e-6 * max(1.0, abs(hi), abs(lo))
    return all(abs(Float64(v[i]) - (lo + (i - 1) * δ)) <= tol for i in eachindex(v))
end

# ...but `SubTab` reconstructs in **single** precision (`ReadVector_sngl` declares `REAL :: x(Nx)`),
# as `x(1) + (i-1) * deltax`. That can land the last point past the endpoint: a 0-120 m grid of 101
# points gives deltax = 1.2f0 and a final value of 120.00001, which is *below the seabed*. KRAKEN
# then prints "Warning in ReadSzRz : Receiver depth exceeds bottom bdry has been moved up" and
# clamps it. Harmless, but it is our grid that provoked it, so the terse form is used only when the
# reconstruction is exact and the vector is written out in full otherwise.
function _subtab_is_exact(v::AbstractVector{<:Real})
    n = length(v)
    lo, hi = Float32(first(v)), Float32(last(v))
    δ = (hi - lo) / (n - 1)
    return lo + Float32(n - 1) * δ <= hi
end

function _write_vector(io::IO, v::AbstractVector{<:Real}, label::AbstractString)
    n = length(v)
    n >= 1 || error("$label must contain at least one value")
    println(io, rpad(string(n), 40), "! N", uppercase(label))
    if _is_uniform(v) && _subtab_is_exact(v)
        # Two endpoints + `/`; SubTab fills in the other n-2 points.
        println(io, rpad(_fmt(first(v)) * "  " * _fmt(last(v)) * " /", 40), "! $label(1) ... $label($n), subtabulated")
    else
        for chunk in Iterators.partition(eachindex(v), 8)
            print(io, join((_fmt(v[i]) for i in chunk), "  "))
            println(io, last(chunk) == lastindex(v) ? " /" : "")
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------
# Splitting an SSP table into media
# ---------------------------------------------------------------------------------------------

"""
    _split_media(ssp, layer_depth) -> Vector{Matrix}

Partition the rows of an `[z cp cs ρ αp αs]` table into one sub-table per medium, given the vector
of medium *bottom* depths.

The rows are consumed in order rather than selected by depth range, because a table straddling an
interface carries two rows at (numerically) the same depth — `one_layer_env` puts the water bottom at
`z0` and the sediment top at `z0 + eps(z0)`. A range test would put both rows in both media; walking
the table in order puts the water row in medium 1 and the sediment row in medium 2, which is what
they mean.

Each sub-table's first and last depths are then snapped exactly onto the medium's interfaces.
`ReadSSP` stops reading a medium when the row depth matches `Z(NSSP)` to within `100 * eps(1f0)`
(≈1.2e-5), so the `eps`-sized offsets above would otherwise be fine — but snapping also removes them
from the *first* row, where medium 2 would otherwise start a hair below the interface that medium 1
ended on.
"""
function _split_media(ssp::AbstractMatrix, layer_depth::AbstractVector)
    size(ssp, 2) == 6 || error("SSP table must have 6 columns [z cp cs ρ αp αs], got $(size(ssp, 2))")
    z = @view ssp[:, 1]
    nrows = length(z)
    media = Vector{Matrix{Float64}}(undef, length(layer_depth))
    first_row = 1
    top = 0.0
    for (m, bottom) in enumerate(layer_depth)
        bottom = Float64(bottom)
        tol = 1e-6 * max(1.0, abs(bottom))
        last_row = first_row
        while last_row < nrows && Float64(z[last_row]) < bottom - tol
            last_row += 1
        end
        if last_row - first_row < 1
            error(
                "Medium $m (depths $top to $bottom m) has $(last_row - first_row + 1) SSP point(s); " *
                "KRAKEN requires at least 2 per medium",
            )
        end
        abs(Float64(z[last_row]) - bottom) <= tol || error(
            "The SSP table ends at $(Float64(z[last_row])) m but medium $m ends at $bottom m " *
            "(from `layers`); the two must agree",
        )
        # The snapping below only removes eps-sized offsets; a profile that genuinely starts
        # somewhere other than the interface is an error, not something to quietly relocate.
        # (Medium 1's interface is 0 m because `get_thickness` measures every layer from the
        # surface, so Kraken.jl cannot represent a first medium starting deeper than that anyway.)
        abs(Float64(z[first_row]) - top) <= tol || error(
            "Medium $m starts at $(Float64(z[first_row])) m but its interface is at $top m; " *
            "the SSP table and `layers` must agree",
        )
        rows = Float64.(ssp[first_row:last_row, :])
        rows[1, 1] = top          # snap onto the interface above ...
        rows[end, 1] = bottom     # ... and the one below
        issorted(@view rows[:, 1]) && allunique(@view rows[:, 1]) ||
            error("Depths within medium $m must be strictly increasing, got $(rows[:, 1])")
        media[m] = rows
        first_row = last_row + 1
        top = bottom
    end
    first_row == nrows + 1 ||
        error("$(nrows - first_row + 1) SSP row(s) below the last medium ($(Float64(top)) m) would be dropped")
    return media
end

# ---------------------------------------------------------------------------------------------
# Environment -> the six-column tables the writer needs
# ---------------------------------------------------------------------------------------------

"""
    _env_tables(env) -> (ssp, layer_depth, halfspace)

Normalize either environment type to the `[z cp cs ρ αp αs]` table, the vector of medium bottom
depths, and the bottom half-space row.

`UnderwaterEnvFORTRAN` keeps the original matrices, so it round-trips shear speeds and attenuations
untouched. `UnderwaterEnv` has already thrown those away — it stores only the sound-speed and
density interpolants — so they are reconstructed as zeros, which is exactly the feature set the
solver currently supports.
"""
_env_tables(env::UnderwaterEnvFORTRAN) = (env.ssp, env.layers[:, 3], env.sspHS[2, :])

function _env_tables(env::UnderwaterEnv)
    # SampledSSP1D/SampledDensity1D store *negated* depths (see their inner constructors).
    z = -env.c.z
    cp = env.c.c
    ρ = env.ρ.ρ
    length(z) == length(cp) == length(ρ) || error(
        "Sound-speed and density profiles must be sampled at the same depths " *
        "(got $(length(z)), $(length(cp)), $(length(ρ)) points)",
    )
    ssp = hcat(z, cp, zero(cp), ρ, zero(cp), zero(cp))
    halfspace = [env.depth, env.cb, zero(env.cb), env.ρb, zero(env.cb), zero(env.cb)]
    return ssp, env.layer_depth, halfspace
end

# ---------------------------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------------------------

"""
    _default_nrd(depth, freq, cmin) -> Int

Number of receiver depths to request by default.

This is not cosmetic: `kraken.exe` tabulates its mode shapes in the `.mod` file on `zTab`, the
*merge of the source and receiver depth vectors* (`MergeVectors` in `kraken.f90`). The receiver grid
is therefore the resolution at which mode shapes come back, and the comparison in task 3.5
interpolates the Julia modes onto it. 20 points per wavelength matches the density KRAKEN picks for
its own finite-difference mesh; the bounds keep short/deep environments sane and cap the `.mod`
record length.
"""
function _default_nrd(depth::Real, freq::Real, cmin::Real)
    per_wavelength = 20 * Float64(depth) * Float64(freq) / Float64(cmin)
    return clamp(ceil(Int, per_wavelength) + 1, 101, 2001)
end

# ---------------------------------------------------------------------------------------------
# The writer
# ---------------------------------------------------------------------------------------------

"""
    env_file_string(env, freq; kwargs...) -> String

Render `env` as the text of a KRAKEN `.env` file. See [`write_env_file`](@ref) for the keywords;
this is the same function without the file, which makes the format easy to test and to eyeball.
"""
function env_file_string(
    env::Union{UnderwaterEnv,UnderwaterEnvFORTRAN},
    freq;
    title::AbstractString="Kraken.jl reference environment",
    topopt::AbstractString=DEFAULT_TOPOPT,
    botopt::AbstractString="A",
    clow=nothing,
    chigh=nothing,
    rmax::Real=10.0,
    sd=nothing,
    rd=nothing,
    nmesh=0,
    sigma::Real=0.0,
    freq0=nothing,
    density_scale::Real=DEFAULT_DENSITY_SCALE,
)
    freqs = freq isa Real ? [Float64(freq)] : Float64.(collect(freq))
    isempty(freqs) && error("At least one frequency is required")
    all(>(0), freqs) || error("Frequencies must be positive, got $freqs")
    # The FREQ record is KRAKEN's *nominal* frequency: it sizes the automatic mesh, which is then
    # rescaled per frequency as `NG * freq / freq0` (`kraken.f90`), and it is the reference frequency
    # for attenuation-unit conversion. For a single-frequency run it is simply the frequency.
    nominal = freq0 === nothing ? first(freqs) : Float64(freq0)

    ssp, layer_depth, halfspace = _env_tables(env)
    media = _split_media(ssp, layer_depth)
    nmedia = length(media)
    depth = Float64(last(layer_depth))

    # NMESH may be given once for every medium or once for all of them. 0 means "let KRAKEN size its
    # own mesh" — anything smaller than half of what it wants is a fatal `Mesh is too coarse`.
    nmesh_vec = nmesh isa Integer ? fill(Int(nmesh), nmedia) : Int.(collect(nmesh))
    length(nmesh_vec) == nmedia || error("`nmesh` must be a scalar or one value per medium ($nmedia)")

    cmin = minimum(minimum(@view m[:, 2]) for m in media)
    cb = Float64(halfspace[2])
    # `ReadEnvironment` rejects `clow >= chigh`. The defaults bracket the band in which modes are
    # trapped — phase speeds from the slowest sound speed up to the half-space speed — which is the
    # same band `bisection` searches in Kraken.jl (`kr ∈ [ω/(0.9999 cb), max(ω/c)]`). `clow` is
    # nudged 1% below the minimum rather than sitting on it, so a mode grazing the slowest layer is
    # strictly inside the interval; no mode exists below `cmin`, so the widening costs nothing.
    clow_v = clow === nothing ? 0.99 * cmin : Float64(clow)
    chigh_v = chigh === nothing ? cb : Float64(chigh)
    clow_v < chigh_v || error("Need clow < chigh, got clow=$clow_v, chigh=$chigh_v")

    sd_v = sd === nothing ? [depth / 4] : Float64.(collect(sd))
    rd_v = if rd === nothing
        collect(range(0.0, depth; length=_default_nrd(depth, maximum(freqs), cmin)))
    else
        Float64.(collect(rd))
    end

    broadband = length(freqs) > 1
    topopt_str = String(topopt)
    if broadband
        # `ReadfreqVec` keys the broadband block off column 6 of the top-option record.
        length(topopt_str) >= 6 &&
            topopt_str[6] != ' ' &&
            topopt_str[6] != 'B' &&
            error("Cannot request a broadband run: column 6 of topopt is already $(topopt_str[6])")
        topopt_str = rpad(topopt_str[1:min(end, 5)], 5) * "B"
    end

    io = IOBuffer()
    println(io, rpad(_fmt_title(title), 40), "! TITLE")
    println(io, rpad(_fmt(nominal), 40), "! FREQ (Hz)")
    println(io, rpad(string(nmedia), 40), "! NMEDIA")
    println(io, rpad("'" * topopt_str * "'", 40), "! TOP OPTIONS")

    for (m, rows) in enumerate(media)
        bottom = rows[end, 1]
        mesh_line = "$(nmesh_vec[m])  $(_fmt(sigma))  $(_fmt(bottom))"
        println(io, rpad(mesh_line, 40), "! NMESH  SIGMA (m)  Z(NSSP)  -- medium $m")
        for i in axes(rows, 1)
            fields = (
                _fmt(rows[i, 1]),
                _fmt(rows[i, 2]),
                _fmt(rows[i, 3]),
                _fmt(rows[i, 4] * density_scale),
                _fmt(rows[i, 5]),
                _fmt(rows[i, 6]),
            )
            line = join(fields, "  ") * " /"
            println(io, i == 1 ? rpad(line, 40) * "! Z  CP  CS  RHO  AP  AS" : line)
        end
    end

    println(io, rpad("'" * String(botopt) * "' " * _fmt(sigma), 40), "! BOTOPT  SIGMA")
    if String(botopt)[1] == 'A'
        hs = join(
            (
                _fmt(depth),
                _fmt(halfspace[2]),
                _fmt(halfspace[3]),
                _fmt(halfspace[4] * density_scale),
                _fmt(halfspace[5]),
                _fmt(halfspace[6]),
            ),
            "  ",
        )
        println(io, rpad(hs * " /", 40), "! ZB  CPB  CSB  RHOB  APB  ASB")
    end

    println(io, rpad(_fmt(clow_v) * "  " * _fmt(chigh_v), 40), "! CLOW  CHIGH (m/s)")
    println(io, rpad(_fmt(rmax), 40), "! RMAX (km)")
    _write_vector(io, sd_v, "SD")
    _write_vector(io, rd_v, "RD")
    broadband && _write_vector(io, freqs, "FREQ")

    return String(take!(io))
end

"""
    write_env_file(path, env, freq; kwargs...) -> String

Serialize `env` into a KRAKEN `.env` file and return the path written.

`path` may be given with or without the `.env` extension — `kraken.exe` takes the *file root* as its
single argument, so it is convenient to build both from the same string. `env` is either an
[`UnderwaterEnv`](@ref) or an [`UnderwaterEnvFORTRAN`](@ref); `freq` is a frequency in Hz, or a
vector of them for a broadband run.

# Keywords
- `title` — the `.env` title record. Truncated to the 72 characters KRAKEN keeps.
- `topopt` — top-option string, default `$(DEFAULT_TOPOPT)` (see [`DEFAULT_TOPOPT`](@ref) for why it
  is C-linear rather than the spline the checked-in files use).
- `botopt` — bottom-option string, default `"A"` (acousto-elastic half-space).
- `clow`, `chigh` — phase-speed limits in m/s. Default to `0.99 * min(c)` and the half-space speed,
  which is the band in which modes are trapped.
- `rmax` — maximum range in km, default `10.0`.
- `sd`, `rd` — source and receiver depths in m. Their union is the depth grid on which `kraken.exe`
  tabulates mode shapes in the `.mod` file, so `rd` controls the resolution of the modes you get
  back; it defaults to a uniform grid over the whole water column at ~20 points per wavelength.
- `nmesh` — finite-difference mesh points, scalar or one per medium. `0` (the default) asks KRAKEN
  to size its own mesh.
- `sigma` — interfacial RMS roughness in m, default `0.0`. Roughness is out of scope for the solver.
- `freq0` — nominal frequency for the `FREQ` record, defaulting to the first entry of `freq`. KRAKEN
  sizes its automatic mesh at `freq0` and rescales it per frequency, so this only matters for a
  broadband run where you want the mesh anchored somewhere other than the first frequency.
- `density_scale` — factor converting the environment's densities to the g/cm³ the format wants;
  defaults to `1e-3` for the kg/m³ the standard environments use.
"""
function write_env_file(path::AbstractString, env, freq; kwargs...)
    file = endswith(path, ".env") ? String(path) : String(path) * ".env"
    open(f -> write(f, env_file_string(env, freq; kwargs...)), file, "w")
    return file
end
