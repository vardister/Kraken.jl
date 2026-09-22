### Reading the `.shd` shade file `field.exe` writes.
#
# A direct-access binary file whose record length is its own first word. Everything is
# little-endian on the platforms the jll builds for, and the layout below is transcribed from
# `Matlab/ReadWrite/read_shd_bin.m`, which is the Acoustics Toolbox's own reader:
#
#   record 1   int32 recl, then 80 characters of title          (each record is 4*recl bytes)
#   record 2   10 characters of PlotType
#   record 3   int32 × 7: Nfreq Ntheta Nsx Nsy Nsz Nrz Nrr, then float64 freq0, float64 atten
#   record 4   float64 × Nfreq   frequencies
#   record 5   float64 × Ntheta  bearings
#   record 6   float64 × Nsx     source x
#   record 7   float64 × Nsy     source y
#   record 8   float32 × Nsz     source depths      ← note: single, where x and y are double
#   record 9   float32 × Nrz     receiver depths    ← likewise
#   record 10  float64 × Nrr     receiver ranges
#   record 11+ float32 × 2Nrr    interleaved re/im, one record per (freq, θ, sz, rz)
#
# The mixed precision in records 8-10 is not a transcription slip: `Pos%Sz` and `Pos%Rz` are
# `REAL` in the Fortran while the ranges are `REAL(KIND=8)`, and the file mirrors the types.

"""
    SHD_FIRST_DATA_RECORD

Index (0-based, as `read_shd_bin.m` counts) of the first record holding pressure. Records 0 through
9 are the header; data starts at 10.
"""
const SHD_FIRST_DATA_RECORD = 10

_shd_seek(io, recl, record) = seek(io, record * 4 * recl)

"""
    read_shd_file(path; freq=nothing) -> NamedTuple

Read the complex pressure field `field.exe` wrote.

Returns `(; p, ranges, zs, zr, freqs, theta, title, plot_type)`, where `p` is indexed
`[source depth, receiver depth, range]` for the selected frequency — the source-depth axis is kept
because the file has one, even though this harness only ever writes a single source.

`freq` picks one of a broadband file's frequencies (nearest match); the default is the first.

Ranges are in metres and depths in metres, as they are in the file.
"""
function read_shd_file(path::AbstractString; freq=nothing)
    isfile(path) || error("No shade file at $path")
    open(path, "r") do io
        recl = read(io, Int32)
        recl > 0 || error("$path: nonsensical record length $recl; is this a .shd file?")
        title = strip(String(read(io, 80)))

        _shd_seek(io, recl, 1)
        plot_type = strip(String(read(io, 10)))

        _shd_seek(io, recl, 2)
        nfreq, ntheta, nsx, nsy, nsz, nrz, nrr = (Int(read(io, Int32)) for _ in 1:7)
        freq0 = read(io, Float64)
        read(io, Float64)   # atten, unused here

        _shd_seek(io, recl, 3)
        freqs = [read(io, Float64) for _ in 1:nfreq]
        _shd_seek(io, recl, 4)
        theta = [read(io, Float64) for _ in 1:ntheta]
        _shd_seek(io, recl, 7)
        zs = [Float64(read(io, Float32)) for _ in 1:nsz]
        _shd_seek(io, recl, 8)
        zr = [Float64(read(io, Float32)) for _ in 1:nrz]
        _shd_seek(io, recl, 9)
        ranges = [read(io, Float64) for _ in 1:nrr]

        # 'irregular' collapses the receiver-depth axis: the depths pair off with the ranges instead
        # of forming a grid, so there is one record per source depth rather than per (sz, rz).
        #
        # `field.exe` writes **ten blanks** here and nothing else: `KrakenField/field.f90` declares
        # `PlotType = '          '` and never assigns it. So the blank is the normal case, not a
        # corrupt file, and a reader that insisted on seeing "rectilin" would reject every shade
        # file field.exe produces. `read_shd_bin.m` handles it the same way, in its `otherwise`
        # branch. Only `field3d` and Scooter fill this field in.
        per_range = startswith(plot_type, "irregular") ? 1 : nrz
        ifreq = freq === nothing ? 1 : argmin(abs.(freqs .- Float64(freq)))

        p = Array{ComplexF64,3}(undef, nsz, per_range, nrr)
        for itheta in 1:ntheta, isz in 1:nsz, irz in 1:per_range
            record =
                SHD_FIRST_DATA_RECORD +
                (ifreq - 1) * ntheta * nsz * per_range +
                (itheta - 1) * nsz * per_range +
                (isz - 1) * per_range +
                (irz - 1)
            _shd_seek(io, recl, record)
            raw = [read(io, Float32) for _ in 1:(2 * nrr)]
            for ir in 1:nrr
                p[isz, irz, ir] = ComplexF64(raw[2ir - 1], raw[2ir])
            end
        end

        return (; p, ranges, zs, zr, freqs, freq0, theta, title, plot_type, nfreq, ntheta)
    end
end

"""
    shd_transmission_loss(shd) -> Array

`-20log10|p|` of a [`read_shd_file`](@ref) result, which is what the toolbox's own `plotshd` plots
and so the quantity to compare a Kraken.jl [`transmission_loss`](@ref) against directly.
"""
shd_transmission_loss(shd) = -20 .* log10.(abs.(shd.p))
