### The reverse of `env_writer.jl`: parse a KRAKEN `.env` file into an `UnderwaterEnv`.
#
# This exists so the environments shipped with the Acoustics Toolbox become test cases. Most of them
# use features Kraken.jl does not implement yet; the point of the reader is that those are *detected
# and named* rather than silently mis-parsed into a plausible-looking environment that then
# "disagrees" with Fortran for reasons nobody can find.
#
# Two pieces of Fortran list-directed input behaviour have to be reproduced exactly, because the
# checked-in Acoustics Toolbox files rely on both:
#
#   * A record may stop early with `/`, leaving the remaining variables at their PREVIOUS values.
#     `ReadSSP` reads into the module-level `alphaR, betaR, rhoR, alphaI, betaI`, which persist from
#     row to row and from medium to medium. So `200.0  1530.29 /` means "depth 200, sound speed
#     1530.29, and everything else as it was on the last row" — not "…and zeros".
#   * Commas are separators just like spaces, and `!` starts a comment.
#
# See `ReadEnvironment`/`ReadTopOpt`/`TopBot` in `misc/ReadEnvironmentMod.f90` and `ReadSSP` in
# `misc/sspMod.f90`.

using Kraken:
    UnderwaterEnv, ATTENUATION_UNIT_CHARS, SSP_INTERPOLATION_CHARS, PressureRelease, RigidBoundary, AcousticHalfspace

"""
    UnsupportedEnvFeature

Raised by [`read_env_file`](@ref) when a `.env` file is valid KRAKEN but uses something Kraken.jl
does not model. `feature` is a short category (suitable for grouping) and `detail` says what was
actually found.

The categories are the backlog. Milestones 5 and 6 cleared attenuation, rigid and vacuum boundaries,
and the n²-linear and cubic-spline interpolators; what remains is either out of scope (elastic layers,
roughness, a top half-space, reflection-coefficient files) or named here so it can be picked up later.
"""
struct UnsupportedEnvFeature <: Exception
    path::String
    feature::String
    detail::String
end

function Base.showerror(io::IO, e::UnsupportedEnvFeature)
    print(io, "UnsupportedEnvFeature: ", e.feature)
    isempty(e.detail) || print(io, " (", e.detail, ")")
    isempty(e.path) || print(io, "\n  in ", e.path)
    return nothing
end

"""
    MalformedEnvFile

Raised when a `.env` file cannot be parsed as a KRAKEN environment at all — as opposed to parsing
fine but using an unsupported feature. Most commonly this is a BELLHOP or SCOOTER input deck, which
shares the first few records and then diverges.
"""
struct MalformedEnvFile <: Exception
    path::String
    detail::String
end

Base.showerror(io::IO, e::MalformedEnvFile) = print(io, "MalformedEnvFile: ", e.detail, "\n  in ", e.path)

# ---------------------------------------------------------------------------------------------
# List-directed tokenizer
# ---------------------------------------------------------------------------------------------

struct _Token
    text::String
    quoted::Bool
end

# Tokens of one line: whitespace and commas separate, `!` starts a comment, `/` ends the record,
# and single quotes group (with '' as an escaped quote).
function _tokenize(line::AbstractString)
    tokens = _Token[]
    i, n = firstindex(line), lastindex(line)
    while i <= n
        c = line[i]
        if isspace(c) || c == ','
            i = nextind(line, i)
        elseif c == '!'
            break
        elseif c == '/'
            push!(tokens, _Token("/", false))
            break
        elseif c == '\''
            j = nextind(line, i)
            buf = IOBuffer()
            while j <= n
                if line[j] == '\''
                    k = nextind(line, j)
                    if k <= n && line[k] == '\''
                        write(buf, '\'')
                        j = nextind(line, k)
                    else
                        j = k
                        break
                    end
                else
                    write(buf, line[j])
                    j = nextind(line, j)
                end
            end
            push!(tokens, _Token(String(take!(buf)), true))
            i = j
        else
            j = i
            while j <= n && !isspace(line[j]) && line[j] ∉ (',', '!', '/')
                j = nextind(line, j)
            end
            push!(tokens, _Token(String(line[i:prevind(line, j)]), false))
            i = j
        end
    end
    return tokens
end

mutable struct _Cursor
    path::String
    lines::Vector{String}
    i::Int
end

_Cursor(path) = _Cursor(path, collect(eachline(path)), 1)

# Next line that carries any token at all. Blank lines and comment-only lines are skipped, exactly
# as a Fortran list-directed READ skips them.
function _next_tokens!(cur::_Cursor)
    while cur.i <= length(cur.lines)
        tokens = _tokenize(cur.lines[cur.i])
        cur.i += 1
        isempty(tokens) || return tokens
    end
    return nothing
end

_exhausted(cur::_Cursor) = cur.i > length(cur.lines)

"""
    _read_values!(cur, n) -> Union{Vector{Float64},Nothing}

Read up to `n` numbers, continuing onto further lines the way a list-directed `READ` does, and
stopping early at `/`. Returns fewer than `n` values when the record was terminated by `/`; the
caller decides what the missing ones mean.
"""
function _read_values!(cur::_Cursor, n::Int)
    values = Float64[]
    while length(values) < n
        tokens = _next_tokens!(cur)
        tokens === nothing && return isempty(values) ? nothing : values
        for token in tokens
            token.text == "/" && return values
            token.quoted && continue
            # Fortran accepts D exponents; Julia does not.
            parsed = tryparse(Float64, replace(token.text, 'D' => 'E', 'd' => 'e'))
            parsed === nothing || push!(values, parsed)
            length(values) >= n && break
        end
    end
    return values
end

function _read_string!(cur::_Cursor)
    tokens = _next_tokens!(cur)
    tokens === nothing && return nothing
    for token in tokens
        token.text == "/" && break
        return token.text
    end
    return ""
end

# ---------------------------------------------------------------------------------------------
# Feature support
# ---------------------------------------------------------------------------------------------

const _SSP_TYPE_NAMES = Dict(
    'C' => "C-linear", 'N' => "n²-linear", 'S' => "cubic spline", 'P' => "PCHIP", 'A' => "analytic"
)
const _BC_NAMES = Dict(
    'V' => "vacuum",
    'A' => "acousto-elastic halfspace",
    'R' => "rigid",
    'F' => "reflection-coefficient file",
    'W' => "write IRC",
    'P' => "precomputed IRC",
)

_describe(table, c) = get(table, c, "'$c'")

# The boundary conditions Kraken.jl solves, by option character. A top 'A' is deliberately absent
# (plan task 6.1). 'F', 'W' and 'P' are reflection-coefficient files, which `kraken.f90` rejects.
const _TOP_BCS = Dict('V' => PressureRelease(), 'R' => RigidBoundary())
const _BOTTOM_BCS = Dict('A' => AcousticHalfspace(), 'V' => PressureRelease(), 'R' => RigidBoundary())

"""
    _ssp_interpolation(ssp, medium_ranges, ssp_type) -> Union{Symbol,Nothing}

The `SampledSSP` interpolation that solves the same problem as KRAKEN's `ssp_type` over this profile,
or `nothing` if there is none.

`C` and `N` map straight across. `S` and `P` are only available as a single whole-column spline, so:

  * where every medium is isovelocity, or every medium is two points, any of KRAKEN's per-medium
    interpolants *is* the straight line, and C-linear is exact;
  * otherwise a single-medium `S` is a spline (up to the end condition, see plan task 6.3), and
    anything else has no equivalent.

Checked in that order, so a two-point single-medium `S` becomes `:c_linear` rather than a spline
through two knots.
"""
function _ssp_interpolation(ssp::AbstractMatrix, medium_ranges, ssp_type::Char)
    ssp_type in ('S', 'P') || return get(SSP_INTERPOLATION_CHARS, ssp_type, nothing)
    two_point = all(length(r) == 2 for r in medium_ranges)
    isovelocity = all(allequal(@view ssp[r, 2]) for r in medium_ranges)
    (two_point || isovelocity) && return :c_linear
    return ssp_type == 'S' && length(medium_ranges) == 1 ? :cubic_spline : nothing
end

function _interpolation_refusal(ssp_type::Char, nmedia::Int)
    ssp_type == 'S' && return "cubic spline over $nmedia media; Kraken.jl splines single-medium profiles only"
    return "$(_describe(_SSP_TYPE_NAMES, ssp_type)) over a varying profile; Kraken.jl has no such interpolator"
end

# ---------------------------------------------------------------------------------------------
# The reader
# ---------------------------------------------------------------------------------------------

"""
    read_env_file(path; strict=true) -> NamedTuple

Parse a KRAKEN `.env` file.

Returns `(; env, ssp, layers, sspHS, freqs, title, topopt, botopt, clow, chigh, rmax, sigmas, nmesh,
atten_units)`, where `env` is an [`UnderwaterEnv`](@ref) ready for `kraken_jl` and `freqs` is the
frequency vector (one entry unless the file declares a broadband run). `nmesh` is the per-medium
`NMESH` record, `0` meaning "let KRAKEN choose"; it does not affect Kraken.jl, but writing a file back
with a different mesh makes `kraken.exe` solve on a different grid.

Throws [`UnsupportedEnvFeature`](@ref) when the file is valid KRAKEN but uses something the solver
does not model, and [`MalformedEnvFile`](@ref) when it is not a KRAKEN environment at all. Pass
`strict=false` to skip the feature checks and get whatever environment the file describes — useful
for inspecting a case, not for cross-validating one.

What `strict` requires, and when each requirement lifts:

| requirement | why |
|---|---|
| `TopOpt(1:1)` is `C`, `N` or `S` (`S` over one medium, see below) | the interpolators `SampledSSP` implements |
| `TopOpt(2:2)` is `V` or `R` | a top half-space is out of scope (plan task 6.1) |
| `TopOpt(4:4)` blank | no THORP / Francois-Garrison / biological volume attenuation |
| `BotOpt(1:1)` is `A`, `V` or `R` | `F` is a reflection-coefficient file, which `kraken.f90` itself refuses |
| `cs == 0` everywhere | the scheme is a scalar pressure equation, not a 4-field elastic one |
| `αs == 0` | shear attenuation needs an elastic layer to live in |
| `sigma == 0` | no interfacial roughness |

The option characters map onto the environment as

| column | character | becomes |
|---|---|---|
| `TopOpt(1:1)` | `C` / `N` / `S` | `env.c.mode` = `:c_linear` / `:n2_linear` / `:cubic_spline` |
| `TopOpt(2:2)` | `V` / `R` | `env.top_bc` = `PressureRelease()` / `RigidBoundary()` |
| `TopOpt(3:3)` | see `ATTENUATION_UNIT_CHARS` | `env.atten_units` |
| `BotOpt(1:1)` | `A` / `V` / `R` | `env.bottom_bc` = `AcousticHalfspace()` / `PressureRelease()` / `RigidBoundary()` |

A perfect (`V` or `R`) bottom has no half-space record, so `sspHS`'s bottom row is zeros, which is
what `TopBot` leaves in `HSBot` too.

**Cubic spline is read as C-linear where the two are the same problem.** `SampledSSP` builds one
spline through the whole column and so refuses a profile with a layer interface in it. KRAKEN splines
each medium separately. Where every medium is isovelocity, or given by exactly two points, KRAKEN's
spline *is* the straight line, and the file is read as `:c_linear`. `Pekeris_AV.env` is such a file:
it declares `S` and matches the C-linear file this repo generates to all ten printed digits. A
single-medium profile is read as `:cubic_spline`. A multi-medium profile that genuinely curves inside
a medium is refused. PCHIP (`P`) gets the same equivalence and is otherwise refused.

Compressional attenuation (`αp`, SSP column 5) is **read and used** — Milestone 5 lifted that
restriction. `TopOpt(3:3)` selects the units it is read in and is parsed even when `strict` is off,
because getting it wrong silently rescales every attenuation in the file; an unknown character is
rejected, exactly as `ReadTopOpt` rejects it.

Densities are converted from the file's g/cm³ to the kg/m³ the standard environments use, so a
`write_env_file` / `read_env_file` round trip is the identity.
"""
function read_env_file(path::AbstractString; strict::Bool=true)
    isfile(path) || error("No .env file at $path")
    cur = _Cursor(String(path))
    malformed(detail) = throw(MalformedEnvFile(String(path), detail))
    unsupported(feature, detail="") = throw(UnsupportedEnvFeature(String(path), feature, detail))

    title = _read_string!(cur)
    title === nothing && malformed("empty file")

    freq0 = _read_values!(cur, 1)
    (freq0 === nothing || isempty(freq0)) && malformed("no frequency record")
    nmedia_v = _read_values!(cur, 1)
    (nmedia_v === nothing || isempty(nmedia_v)) && malformed("no NMEDIA record")
    # A non-integral value here means the record was not NMEDIA at all — BELLHOP3D decks put a
    # quoted option string second, so the reader slides onto an angle fan and lands on -14.66.
    isinteger(nmedia_v[1]) || malformed("NMEDIA is not an integer ($(nmedia_v[1])); not a KRAKEN deck")
    nmedia = Int(nmedia_v[1])
    (1 <= nmedia <= 100) || malformed("implausible NMEDIA = $nmedia")

    topopt = _read_string!(cur)
    topopt === nothing && malformed("no top-option record")
    topopt = rpad(topopt, 6)

    ssp_type = topopt[1]
    # The attenuation units live in column 3 and govern how column 6 of every SSP record is read, so
    # this is parsed whether or not `strict` is on — misreading dB/wavelength as nepers/m is exactly
    # the "plausible-looking environment that then disagrees with Fortran" this reader exists to
    # prevent. `'m'` (dB/m with a power law) needs two extra per-medium parameters that this reader
    # does not model; it is named as unsupported rather than quietly folded into `'M'`.
    atten_char = topopt[3]
    atten_char == 'm' &&
        unsupported("power-law attenuation", "top option 3 = 'm' needs the per-medium beta and fT parameters")
    atten_units = get(ATTENUATION_UNIT_CHARS, atten_char, nothing)
    atten_units === nothing && malformed("unknown attenuation-units option '$atten_char' in top options '$topopt'")

    if strict
        haskey(_SSP_TYPE_NAMES, ssp_type) || malformed("unknown SSP interpolation option '$ssp_type'")
        # 'A' has no SSP records at all — KRAKEN generates the profile from `ANALYT`. That is an
        # unsupported feature, not a malformed file, and it has to be caught before the parse loop
        # walks into the bottom-option record looking for SSP rows.
        ssp_type == 'A' && unsupported("analytic SSP", "the profile is generated by ANALYT, not tabulated")
        # The interpolation check needs the SSP itself and so happens after the media are read.
        haskey(_TOP_BCS, topopt[2]) || unsupported("top boundary", _describe(_BC_NAMES, topopt[2]))
        topopt[4] in (' ', '\0') || unsupported("added volume attenuation", "top option 4 = '$(topopt[4])'")
    end
    # An acoustic halfspace above the first medium adds a record here that we would otherwise
    # mis-read as a mesh line, so this check is not optional even when `strict` is off.
    topopt[2] == 'A' && unsupported("top boundary", "acousto-elastic halfspace above the surface")
    top_bc = get(_TOP_BCS, topopt[2], PressureRelease())

    # Per-medium: the mesh line, then SSP rows until the depth matches Z(NSSP).
    # `alphaR, betaR, rhoR, alphaI, betaI` start at KRAKEN's own defaults and carry forward.
    carried = [1500.0, 0.0, 1.0, 0.0, 0.0]     # cp, cs, ρ, αp, αs
    rows = Vector{Vector{Float64}}()
    layer_depths = Float64[]
    sigmas = Float64[]
    medium_ranges = UnitRange{Int}[]
    nmesh = Int[]
    for medium in 1:nmedia
        mesh = _read_values!(cur, 3)
        (mesh === nothing || length(mesh) < 3) && malformed("truncated mesh record for medium $medium")
        isinteger(mesh[1]) || malformed("NMESH for medium $medium is not an integer ($(mesh[1]))")
        push!(nmesh, Int(mesh[1]))
        push!(sigmas, mesh[2])
        bottom = mesh[3]
        push!(layer_depths, bottom)
        first_row = length(rows) + 1
        while true
            values = _read_values!(cur, 6)
            (values === nothing || isempty(values)) && malformed("SSP of medium $medium never reaches $bottom m")
            z = values[1]
            for k in 1:5
                length(values) >= k + 1 && (carried[k] = values[k + 1])
            end
            push!(rows, [z; carried...])
            # ReadSSP's own stopping rule.
            abs(z - bottom) < 100 * eps(1.0f0) && break
            length(rows) - first_row > 10_000 && malformed("runaway SSP in medium $medium")
        end
        length(rows) - first_row < 1 && malformed("medium $medium has fewer than 2 SSP points")
        push!(medium_ranges, first_row:length(rows))
    end

    bot_tokens = _next_tokens!(cur)
    bot_tokens === nothing && malformed("no bottom-option record")
    botopt = ""
    bot_sigma = 0.0
    # `READ( ENVFile, * ) BotOpt, Sigma` is a list-directed read of exactly two items: it consumes the
    # option string and the *first* number, and **discards the rest of the record**. Several toolbox
    # files rely on that — `SedAtten/calibS_0.6dB.env` writes `'A'  0.0 2.5 2000`, where the trailing
    # two are power-law parameters that only a `'m'` attenuation unit would go back for. Taking the
    # last number instead of the first read that file's roughness as 2000 m and rejected it.
    got_sigma = false
    for token in bot_tokens
        token.text == "/" && break
        if token.quoted && isempty(botopt)
            botopt = token.text
        elseif !got_sigma
            parsed = tryparse(Float64, token.text)
            if parsed !== nothing
                bot_sigma = parsed
                got_sigma = true
            end
        end
    end
    isempty(botopt) && malformed("bottom-option record has no option string")
    botopt = rpad(botopt, 2)
    push!(sigmas, bot_sigma)

    halfspace = zeros(6)
    bottom_bc = get(_BOTTOM_BCS, botopt[1], AcousticHalfspace())
    if botopt[1] == 'A'
        values = _read_values!(cur, 6)
        # Only the depth is required: `TopBot` reads into the same carried-forward variables as the
        # SSP rows, so `70.0 /` legitimately means "the half-space is whatever the last row was".
        (values === nothing || isempty(values)) && malformed("empty bottom half-space record")
        halfspace[1] = values[1]
        for k in 1:5
            length(values) >= k + 1 && (carried[k] = values[k + 1])
        end
        halfspace[2:6] = carried
    elseif strict && !haskey(_BOTTOM_BCS, botopt[1])
        unsupported("bottom boundary", _describe(_BC_NAMES, botopt[1]))
    end

    ssp = reduce(vcat, permutedims.(rows))
    ssp_interp = _ssp_interpolation(ssp, medium_ranges, ssp_type)
    if ssp_interp === nothing
        strict && unsupported("SSP interpolation", _interpolation_refusal(ssp_type, nmedia))
        # Not strict: build *something* inspectable. The linear profile through the same samples is
        # the nearest thing Kraken.jl can represent.
        ssp_interp = :c_linear
    end
    if strict
        maximum(@view ssp[:, 3]) > 0 && unsupported("elastic layer", "shear speed cs > 0 in the water column")
        halfspace[3] > 0 && unsupported("elastic layer", "shear speed cs > 0 in the bottom half-space")
        # Compressional attenuation (column 5) is Milestone 5 and is now read and used. Shear
        # attenuation (column 6) is not: it only means anything in an elastic layer, and those are
        # already rejected above by the `cs > 0` check — so a nonzero αs here is loss the solver
        # would silently drop rather than a feature it supports.
        shear_atten = max(maximum(@view ssp[:, 6]), halfspace[6])
        shear_atten > 0 && unsupported("elastic layer", "shear attenuation αs = $shear_atten with no shear speed")
        roughness = maximum(abs, sigmas)
        roughness > 0 && unsupported("interfacial roughness", "max sigma = $roughness m")

        # `AcousticProblemProperties` asserts `maxsoundspeed(env.c) < env.cb`: the finite-difference
        # solver only looks for modes trapped *above* the half-space, so a bottom that is not the
        # fastest medium has no trapped spectrum to find and the assertion fires. Common in the
        # toolbox's arctic and surface-duct cases, where the energy is ducted near the surface and
        # the interesting modes are leaky ones. Milestone 5's stretch task (krakenc parity) is what
        # would make these comparable.
        # A perfect bottom has no speed of its own and traps every mode, so it is exempt.
        fastest = maximum(@view ssp[:, 2])
        bottom_bc isa AcousticHalfspace &&
            halfspace[2] <= fastest &&
            unsupported(
                "bottom half-space is not the fastest medium",
                "cb = $(halfspace[2]) m/s vs max water speed $fastest m/s; no trapped modes",
            )

        # `get_thickness` measures every layer from the surface, so the profile has to start there.
        abs(ssp[1, 1]) > 1e-6 &&
            unsupported("profile does not start at the surface", "first SSP depth is $(ssp[1, 1]) m")
    end

    # Everything past here is optional: BELLHOP and SCOOTER decks share the records above and then
    # diverge, so a file that stops making sense here still yields a usable environment.
    clow, chigh, rmax = nothing, nothing, nothing
    speeds = _read_values!(cur, 2)
    if speeds !== nothing && length(speeds) == 2 && speeds[1] < speeds[2]
        clow, chigh = speeds
        r = _read_values!(cur, 1)
        r === nothing || isempty(r) || (rmax = r[1])
    end

    # The broadband frequency vector trails the source/receiver depth records, which this reader
    # does not model — the frequency to run at is the caller's choice anyway, and `run_fortran_kraken`
    # takes it as an argument. The nominal frequency is what the file leads with.
    freqs = [freq0[1]]

    # Convert g/cm³ back to the kg/m³ the standard environments use.
    ssp[:, 4] .*= 1 / DEFAULT_DENSITY_SCALE
    halfspace[4] /= DEFAULT_DENSITY_SCALE

    # `UnderwaterEnv` interpolates the SSP with DataInterpolations, which needs strictly increasing
    # depths — but an interface contributes two rows at the same depth. The standard environment
    # builders offset the second by `eps`; do the same here.
    for i in 2:size(ssp, 1)
        ssp[i, 1] <= ssp[i - 1, 1] && (ssp[i, 1] = nextfloat(ssp[i - 1, 1]))
    end

    layers = hcat(zeros(nmedia), zeros(nmedia), layer_depths)
    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0;
        halfspace[1] halfspace[2] halfspace[3] halfspace[4] halfspace[5] halfspace[6]
    ]
    ssp[end, 1] = layer_depths[end]

    env = UnderwaterEnv(ssp, layers, sspHS; atten_units, ssp_interp, top_bc, bottom_bc)
    return (; env, ssp, layers, sspHS, freqs, title, topopt, botopt, clow, chigh, rmax, sigmas, nmesh, atten_units)
end

"""
    categorize_env_tree(dir) -> NamedTuple

Run [`read_env_file`](@ref) over every `.env` file under `dir` and group the results.

Returns `(; supported, unsupported, malformed, total)`, where `supported` is a vector of paths and
the other two are `Dict{String,Vector{String}}` keyed by the reason. Printing the returned value
gives the categorized report — the concrete backlog for Milestones 5 and 6.
"""
function categorize_env_tree(dir::AbstractString)
    supported = String[]
    unsupported = Dict{String,Vector{String}}()
    malformed = Dict{String,Vector{String}}()
    total = 0
    for (root, _, files) in walkdir(dir), file in files
        endswith(file, ".env") || continue
        total += 1
        path = joinpath(root, file)
        rel = relpath(path, dir)
        try
            read_env_file(path)
            push!(supported, rel)
        catch e
            if e isa UnsupportedEnvFeature
                push!(get!(unsupported, e.feature, String[]), rel)
            elseif e isa MalformedEnvFile
                push!(get!(malformed, e.detail, String[]), rel)
            else
                push!(get!(malformed, "reader error: $(typeof(e))", String[]), rel)
            end
        end
    end
    return (; supported, unsupported, malformed, total)
end

"""
    print_env_tree_report(io, report)

Print the summary produced by [`categorize_env_tree`](@ref), most common blocker first.
"""
function print_env_tree_report(io::IO, report)
    println(io, "Scanned $(report.total) .env files")
    println(io, "  supported by Kraken.jl today: $(length(report.supported))")
    for (label, group) in (("unsupported feature", report.unsupported), ("not a KRAKEN environment", report.malformed))
        total = sum(length, values(group); init=0)
        println(io, "  $label: $total")
        for (reason, files) in sort(collect(group); by=p -> -length(p[2]))
            println(io, "    ", rpad(reason, 46), lpad(length(files), 4), "   e.g. ", first(files))
        end
    end
    return nothing
end
