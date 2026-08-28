### Readers for the two things `kraken.exe` produces: the binary `.mod` mode file and the `Group
### Speed` table in the `.prt` print file.
#
# The `.mod` layout is a Fortran DIRECT-access file: fixed-length records of `LRecordLength`
# *longwords* (4 bytes), addressed by number. Record 1 opens with that length, so the first four
# bytes of the file tell you how to seek everything else. Written by `Vector`/`Solve` in
# `Kraken/kraken.f90`:
#
#     rec 1   LRecordLength, Title(80), Nfreq, NMedia, NzTab, NzTab
#     rec 2   ( N(m), Material(m) )            for each acoustic medium   Int32 + 8 chars
#     rec 3   ( Depth(m), rho(m) )             for each acoustic medium   2 x Float32
#     rec 4   freqVec( 1 : Nfreq )                                        Float64
#     rec 5   zTab( 1 : NzTab )                                           Float32
#
# then, once per frequency, starting at record 6:
#
#     rec r         M                          number of modes           Int32
#     rec r+1       top and bottom halfspace descriptors
#     rec r+1+mode  PhiTab( 1 : Nmat )         mode shape                ComplexF32
#     rec r+1+M+j   k( ... )                   wavenumbers, j = 1, 2, …   ComplexF32
#     next          r += 3 + M + (2M-1) / LRecordLength
#
# Adapted from `AcousticsToolbox.jl`'s `_read_mod` (the cleaner of the two existing readers) with
# the multi-frequency record stepping from `KrakenFortran.jl/src/io.jl`. Two things both of those
# get wrong and this one does not: the wavenumbers spill across *several* records once `2M` exceeds
# `LRecordLength`, and `Nfreq > 1` needs `M` re-read per frequency because the stride depends on it.

using Printf

"""
Precision note, and it matters for choosing tolerances — neither file dominates the other.

The `.mod` stores wavenumbers as `CMPLX(k)`: default `KIND=4`, i.e. **single** precision throughout,
about 7 significant digits on both the real and imaginary parts.

The `.prt` prints the same complex number across two columns with different formats —
`'( I5, G18.10, G10.2, … )'`, where a Fortran `COMPLEX` consumes two edit descriptors. So `k` (the
real part) gets **10** significant digits, and `alpha` (the imaginary part) gets **2**.

Hence:

| quantity | better source | why |
|---|---|---|
| `Re(kᵣ)` | `.prt` | 10 printed digits beat single precision |
| `Im(kᵣ)` — attenuation | `.mod` | 7 digits beat 2 printed ones |
| mode shapes | `.mod` | the `.prt` has none |

Measured on `Pekeris_AV_BroadBand`: the two agree on `Re(kᵣ)` to 6e-8 (the single-precision floor)
and on `Im(kᵣ)` to only ~4e-3 relative, which is exactly what two printed digits buy. Milestone 5
tests attenuation, so it wants the `.mod`.
"""
const MOD_WAVENUMBER_DIGITS = 7

# ---------------------------------------------------------------------------------------------
# .mod
# ---------------------------------------------------------------------------------------------

# Read `n` items of type `T` starting at a record boundary. Records are addressed 0-based here
# (record 1 of the Fortran file is offset 0), which is what makes the arithmetic above line up.
function _read_record(io::IO, reclen::Int, rec::Int, ::Type{T}, n::Integer) where {T}
    seek(io, rec * reclen)
    return read!(io, Vector{T}(undef, n))
end

"""
    read_mod_file(path; freq=nothing) -> NamedTuple

Parse a KRAKEN binary `.mod` file.

Returns `(; ϕ, kᵣ, depths, nmodes, freq, freqs, title, nmedia)`:

- `ϕ` — `nmat x nmodes` matrix of mode shapes, tabulated at `depths`.
- `kᵣ` — horizontal wavenumbers, one per mode. Single precision in the file; see
  [`MOD_WAVENUMBER_DIGITS`](@ref).
- `depths` — the depth grid the modes are sampled on. This is `zTab`, the merge of the source and
  receiver depth vectors from the `.env` file, *not* the solver's internal mesh.
- `nmodes` — number of modes at the selected frequency.
- `freq` — the frequency actually read; `freqs` — every frequency in the file.

`freq` selects which frequency block to read from a broadband file (the nearest one). It defaults to
the first, which is the only one for a single-frequency run.
"""
function read_mod_file(path::AbstractString; freq=nothing)
    isfile(path) || error("No .mod file at $path")
    return open(path, "r") do io
        lrecl = Int(read(io, Int32))                      # record length in 4-byte longwords
        lrecl > 0 || error("Bad .mod file $path: record length $lrecl")
        reclen = 4 * lrecl
        title = strip(String(read(io, 80)))
        nfreq = Int(read(io, Int32))
        nmedia = Int(read(io, Int32))
        ntot = Int(read(io, Int32))
        nmat = Int(read(io, Int32))
        (nfreq > 0 && nmedia > 0 && ntot > 0) ||
            error("Bad .mod header in $path: nfreq=$nfreq, nmedia=$nmedia, ntot=$ntot, nmat=$nmat")

        freqs = Float64.(_read_record(io, reclen, 3, Float64, nfreq))
        depths = Float64.(_read_record(io, reclen, 4, Float32, ntot))

        ifreq = freq === nothing ? 1 : argmin(abs.(freqs .- Float64(freq)))

        # Walk the per-frequency blocks. The stride depends on M, so every block up to the one we
        # want has to be visited -- there is no index.
        rec = 5
        nmodes = 0
        for i in 1:ifreq
            seek(io, rec * reclen)
            nmodes = Int(read(io, Int32))
            nmodes >= 0 || error("Bad mode count $nmodes in $path at frequency $(freqs[i])")
            i == ifreq && break
            rec += 3 + nmodes + (2 * nmodes - 1) ÷ lrecl
        end
        nmodes == 0 && error(
            "$path contains no modes at $(freqs[ifreq]) Hz. kraken.exe writes a stub .mod and " *
            "reports 'No modes for given phase speed interval' -- widen CLOW/CHIGH.",
        )

        ϕ = Matrix{ComplexF64}(undef, nmat, nmodes)
        for m in 1:nmodes
            ϕ[:, m] = _read_record(io, reclen, rec + 1 + m, ComplexF32, nmat)
        end

        # The wavenumbers spill across ceil(2M / LRecordLength) records, LRecordLength/2 complex
        # values at a time. Reading them all from one record -- as both existing readers do -- is
        # correct only while 2M <= LRecordLength.
        kᵣ = Vector{ComplexF64}(undef, nmodes)
        per_record = lrecl ÷ 2
        per_record > 0 || error("Bad .mod file $path: record length $lrecl cannot hold a wavenumber")
        first_mode = 1
        block = 1
        while first_mode <= nmodes
            last_mode = min(nmodes, first_mode + per_record - 1)
            chunk = _read_record(io, reclen, rec + 1 + nmodes + block, ComplexF32, last_mode - first_mode + 1)
            kᵣ[first_mode:last_mode] = chunk
            first_mode = last_mode + 1
            block += 1
        end

        return (; ϕ, kᵣ, depths, nmodes, freq=freqs[ifreq], freqs, title, nmedia)
    end
end

# ---------------------------------------------------------------------------------------------
# .prt group-speed table
# ---------------------------------------------------------------------------------------------

"""
The `Group Speed` table is **subsampled** when there are many modes.

`kraken.f90` prints it with `DO mode = 1, M, MAX( 1, M / 30 )` — "print every mode unless there are
an awful lot" — so a run with 102 modes reports every third one, 34 rows. That is why [`read_grp`](@ref)
returns the mode *indices* alongside the values rather than assuming row `i` is mode `i`, and why
`length(grp.m) < nmodes` is normal rather than a parse failure.
"""
const GRP_MAX_ROWS_PRINTED = 30

# One table row: `WRITE( PRTFile, '( I5, G18.10, G10.2, G18.10, G14.6 )' ) mode, k, omega/k, VG`.
# `k` is COMPLEX, so it consumes two edit descriptors -- hence five numbers per row, with the
# second and third being the real and imaginary parts of the wavenumber.
function _parse_grp_row(line::AbstractString)
    fields = split(line)
    length(fields) == 5 || return nothing
    mode = tryparse(Int, fields[1])
    mode === nothing && return nothing
    values = map(f -> tryparse(Float64, f), fields[2:5])
    any(isnothing, values) && return nothing
    return (mode, values[1], values[2], values[3], values[4])
end

"""
    read_grp_blocks(path) -> Vector{NamedTuple}

Every `Group Speed` table in a `.prt` file, in file order — one per frequency for a broadband run.
Each entry is `(; freq, m, kᵣ, phase_speed, v)`. See [`read_grp`](@ref) for the usual entry point.
"""
function read_grp_blocks(path::AbstractString)
    isfile(path) || error("No .prt file at $path")
    blocks = NamedTuple[]
    current = NaN                     # frequency the next table belongs to
    in_table = false
    m, kr, cp, vg = Int[], ComplexF64[], Float64[], Float64[]

    flush_block!() =
        if !isempty(m)
            push!(blocks, (; freq=current, m=copy(m), kᵣ=copy(kr), phase_speed=copy(cp), v=copy(vg)))
            empty!(m)
            empty!(kr)
            empty!(cp)
            empty!(vg)
        end

    for line in eachline(path)
        # Both the nominal frequency (single-frequency runs) and the per-frequency banner
        # (broadband runs) look like `... Frequency = <value> Hz`.
        freq_match = match(r"Frequency\s*=\s*([-+0-9.EeDd]+)", line)
        if freq_match !== nothing
            flush_block!()
            in_table = false
            parsed = tryparse(Float64, replace(freq_match.captures[1], 'D' => 'E', 'd' => 'e'))
            parsed === nothing || (current = parsed)
            continue
        end
        if occursin("Group Speed", line)
            flush_block!()
            in_table = true
            continue
        end
        in_table || continue
        row = _parse_grp_row(line)
        if row === nothing
            # The units line and the blank line after the table both land here; the table has ended
            # only once we have actually collected something.
            isempty(m) || (in_table = false)
            continue
        end
        push!(m, row[1])
        push!(kr, complex(row[2], row[3]))
        push!(cp, row[4])
        push!(vg, row[5])
    end
    flush_block!()
    return blocks
end

"""
    has_group_speeds(grp) -> Bool

Whether a table returned by [`read_grp`](@ref) actually carries group speeds.

**`kraken.exe` as shipped by `AcousticsToolbox_jll` v2025.9 writes `0.00000` in the Group Speed
column for every mode.** Its wavenumbers and mode shapes are correct and agree digit for digit with
a local Acoustics Toolbox build; only `VG` is lost. The same jll's `krakenc.exe` is fine, and so is
the 2023 OALIB build of both binaries — so this is a regression in that one binary, not a property
of the format or of the environment.

A group-speed comparison therefore has to run `complex=true`, or point `KRAKEN_FORTRAN_BIN` at a
local build. Checking this predicate is what keeps that from surfacing as a spurious "group speeds
differ by 100%".
"""
has_group_speeds(grp) = !isempty(grp.v) && any(!iszero, grp.v)

"""
    read_grp(path; freq=nothing) -> NamedTuple

Parse the `Group Speed` table out of a `.prt` file, returning `(; m, kᵣ, v, phase_speed, freq)`.

`m` is the vector of **mode indices** the rows describe. It is not necessarily `1:nmodes` — see
[`GRP_MAX_ROWS_PRINTED`](@ref). `kᵣ` here is more precise than the `.mod` file's copy of the same
quantity (10 printed digits versus single precision), so this is the better reference for a
wavenumber comparison; `v` is the group speed in m/s, which is what Milestone 4 validates
AD-computed `dω/dkᵣ` against.

`freq` picks a frequency block from a broadband run; it defaults to the first table in the file.
"""
function read_grp(path::AbstractString; freq=nothing)
    blocks = read_grp_blocks(path)
    isempty(blocks) && error("No Group Speed table found in $path")
    block = if freq === nothing
        first(blocks)
    else
        blocks[argmin(abs.([b.freq for b in blocks] .- Float64(freq)))]
    end
    return (; m=block.m, kᵣ=block.kᵣ, v=block.v, phase_speed=block.phase_speed, freq=block.freq)
end
