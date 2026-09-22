### `UnderwaterEnv` geometry → a KRAKEN `.flp` field-parameter file.
#
# `field.exe` reads two files: the `.mod` that `kraken.exe` wrote, and a `.flp` describing where to
# evaluate the field. The `.flp` is list-directed throughout and its records are, in order
# (`KrakenField/field.f90`):
#
#     SHDTitle          quoted string
#     Opt               quoted string; only characters 1-4 are read
#     MLimit            how many modes to include
#     NProf, RProf      profile ranges in km — 1 and 0.0 for a range-independent case
#     NRr,   Rr         receiver ranges, km
#     NSz,   Sz         source depths, m
#     NRz,   Rz         receiver depths, m
#     NRro,  Rro        receiver range offsets (array tilt), m — must have NRro == NRz
#
# Every vector goes through `ReadVector`, so `_write_vector` from `env_writer.jl` applies unchanged,
# including its care about `SubTab`'s single-precision reconstruction.
#
# `Opt` is four positions and each one is checked:
#
#   1  X / R / S      line source / point source / point source without the √r spreading. `R`.
#   2  C / A          coupled or adiabatic modes. Only read when NProf > 1, but the checked-in
#                     files all write `A`, so we do too.
#   3  O / ' ' / *    source beam pattern. Blank is omnidirectional; `*` reads a `.sbp` file.
#   4  C / ' ' / I    coherent or incoherent mode addition — and *nothing else*, which is the whole
#                     content of plan task 7.2's finding. There is no `S`.

"""
    FIELD_OPT_POSITIONS

What each character of a `.flp` `Opt` string selects, as `field.f90` reads it. Position 4 is the
mode-summation code, and `'C'`, `' '` and `'I'` are the only values it accepts.
"""
const FIELD_OPT_POSITIONS = (
    coordinates=('X', 'R', 'S'), modes=('C', 'A'), beam_pattern=('O', ' ', '*'), summation=('C', ' ', 'I')
)

# The `Kraken` summation modes, mapped onto the `.flp` characters. `:semicoherent` is deliberately
# absent from both — see the note at the top of `src/kraken_field.jl`.
const SUMMATION_CODE = Dict(:coherent => 'C', :incoherent => 'I')

"""
    flp_file_string(ranges, zs, zr; mode=:coherent, nmodes=9999, title="Kraken.jl field", tilt=nothing)

The text of a `.flp` field-parameter file for a range-independent case.

`ranges` are in **metres** here and written out in kilometres, because that is what every other
depth and range in this package is measured in and converting in one place beats remembering which
unit a given argument wants. `zs` is a single source depth, `zr` a vector of receiver depths, both
in metres.

`tilt` gives the per-receiver range offsets in metres (array tilt); the default is no tilt. It must
have one entry per receiver depth, which is a constraint `field.exe` itself enforces.
"""
function flp_file_string(
    ranges, zs, zr; mode::Symbol=:coherent, nmodes::Integer=9999, title::AbstractString="Kraken.jl field", tilt=nothing
)
    rs = Float64.(collect(ranges))
    zrv = Float64.(collect(zr))
    isempty(rs) && error("At least one receiver range is required")
    isempty(zrv) && error("At least one receiver depth is required")
    all(>=(0), rs) || error("Receiver ranges must be non-negative")

    haskey(SUMMATION_CODE, mode) || error(
        "`mode = :$mode` has no .flp code. field.exe accepts only coherent and incoherent " *
        "addition; see src/kraken_field.jl on why there is no third.",
    )
    offsets = tilt === nothing ? zeros(length(zrv)) : Float64.(collect(tilt))
    length(offsets) == length(zrv) || error(
        "field.exe requires one range offset per receiver depth (NRro == NRz); got " *
        "$(length(offsets)) offsets for $(length(zrv)) depths",
    )

    io = IOBuffer()
    println(io, _fmt_title(title))
    # Positions 1-4: point source, adiabatic, omnidirectional, and the summation code.
    println(io, rpad("'RA " * SUMMATION_CODE[mode] * "'", 40), "! OPT: coords, couple/adiab, beam pattern, C/I")
    println(io, rpad(string(nmodes), 40), "! MLimit — modes to include")
    # Range-independent: one profile, and field.f90 insists the first sits at 0 km.
    _write_vector(io, [0.0], "RProf")
    _write_vector(io, rs ./ 1000, "Rr")
    _write_vector(io, [Float64(zs)], "Sz")
    _write_vector(io, zrv, "Rz")
    _write_vector(io, offsets, "Rro")
    return String(take!(io))
end

"""
    write_flp_file(path, ranges, zs, zr; kwargs...) -> String

Write [`flp_file_string`](@ref) to `path` (`.flp` is appended when missing) and return the path.
"""
function write_flp_file(path::AbstractString, ranges, zs, zr; kwargs...)
    file = endswith(path, ".flp") ? String(path) : String(path) * ".flp"
    open(f -> write(f, flp_file_string(ranges, zs, zr; kwargs...)), file, "w")
    return file
end
