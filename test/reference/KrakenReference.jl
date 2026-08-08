"""
    KrakenReference

Test-only harness for running **unmodified** Fortran KRAKEN as a reference oracle.

This is deliberately not part of the `Kraken` package: validation machinery is not public API. It
lives under `test/reference/` and is `include`d by the test suite.

Binaries come from [`AcousticsToolbox_jll`](https://github.com/JuliaBinaryWrappers/AcousticsToolbox_jll.jl),
which ships prebuilt `kraken.exe`/`krakenc.exe` for every platform — so CI needs no Fortran
toolchain. Set `KRAKEN_FORTRAN_BIN` to a directory containing `kraken.exe` to compare against a
local build of a specific toolbox version instead.

```julia
include("reference/KrakenReference.jl")
KrakenReference.fortran_available()   # false => tests skip rather than error
KrakenReference.describe()            # one line naming the resolved binary and where it came from
```

```julia
KrakenReference.write_env_file(joinpath(dir, "case"), UnderwaterEnv(pekeris_env()...), 100.0)
```

The `.mod`/`.prt` readers, runner and comparison utility are added by plan tasks 3.3-3.5 and
`include`d at the bottom of this module.
"""
module KrakenReference

using AcousticsToolbox_jll

export fortran_available, write_env_file

"""
Name of the environment variable that overrides the jll binaries with a local Acoustics Toolbox
build. Should point at a *directory* containing `kraken.exe`, e.g.
`/Users/you/programs/AcousticsToolboxOALIB/bin`.
"""
const BIN_ENV_VAR = "KRAKEN_FORTRAN_BIN"

# Both binaries are named `.exe` on every platform — that is the Acoustics Toolbox convention, not a
# Windows-ism, and the jll artifact follows it on macOS and Linux too.
_exe_name(complex::Bool) = complex ? "krakenc.exe" : "kraken.exe"

"""
    override_dir() -> Union{String,Nothing}

Directory named by `KRAKEN_FORTRAN_BIN`, or `nothing` when the variable is unset or empty.
Does not check that the directory exists — [`binary_source`](@ref) does that.
"""
function override_dir()
    dir = get(ENV, BIN_ENV_VAR, "")
    return isempty(dir) ? nothing : dir
end

"""
    jll_path(; complex=false) -> Union{String,Nothing}

Path to the jll's `kraken.exe` (or `krakenc.exe`), or `nothing` if the jll did not build for this
platform. Wrapped in a `try` because a jll that failed to build still loads — it just has no
products, and touching them throws.
"""
function jll_path(; complex::Bool=false)
    try
        AcousticsToolbox_jll.is_available() || return nothing
        path = complex ? AcousticsToolbox_jll.krakenc_path : AcousticsToolbox_jll.kraken_path
        return isfile(path) ? path : nothing
    catch
        return nothing
    end
end

"""
    binary_source(; complex=false) -> Symbol

Where the reference binary would come from: `:override` (`KRAKEN_FORTRAN_BIN`), `:jll`, or `:none`.

The override wins when it resolves, so pointing `KRAKEN_FORTRAN_BIN` at a local build is enough to
switch the whole suite over. A set-but-unusable override falls back to the jll rather than failing:
the point of the variable is to *prefer* a local build, not to make the suite fragile.
"""
function binary_source(; complex::Bool=false)
    dir = override_dir()
    if dir !== nothing && isfile(joinpath(dir, _exe_name(complex)))
        return :override
    end
    return jll_path(; complex=complex) === nothing ? :none : :jll
end

"""
    binary_path(; complex=false) -> String

Absolute path to the resolved `kraken.exe`/`krakenc.exe`. Throws if neither source resolves — call
[`fortran_available`](@ref) first if you want to skip instead.
"""
function binary_path(; complex::Bool=false)
    source = binary_source(; complex=complex)
    if source === :override
        return abspath(joinpath(override_dir()::String, _exe_name(complex)))
    elseif source === :jll
        return jll_path(; complex=complex)::String
    end
    dir = override_dir()
    hint = if dir === nothing
        "AcousticsToolbox_jll has no binaries for this platform, and $BIN_ENV_VAR is not set."
    else
        "$BIN_ENV_VAR is set to $(repr(dir)), which contains no $(_exe_name(complex)), and " *
        "AcousticsToolbox_jll has no binaries for this platform."
    end
    return error("No Fortran $(_exe_name(complex)) available. " * hint)
end

"""
    kraken_cmd(; complex=false) -> Cmd

A runnable command for the reference binary, taking the `.env` file root as its single argument.

For the jll this is the generated wrapper rather than a bare path, so the environment the binary
needs (library search paths, `PATH`) is set up for it.
"""
function kraken_cmd(; complex::Bool=false)
    source = binary_source(; complex=complex)
    if source === :override
        return Cmd([binary_path(; complex=complex)])
    elseif source === :jll
        return complex ? AcousticsToolbox_jll.krakenc() : AcousticsToolbox_jll.kraken()
    end
    return binary_path(; complex=complex)  # throws with the descriptive message
end

"""
    fortran_available(; complex=false) -> Bool

Whether a reference binary can be resolved. The cross-validation tests gate on this so a platform
without binaries skips with a message instead of erroring.
"""
fortran_available(; complex::Bool=false) = binary_source(; complex=complex) !== :none

"""
    describe() -> String

One-line summary of what the harness resolved, for test output and for debugging a machine where
the wrong binary is being picked up.
"""
function describe()
    parts = String[]
    for (label, complex) in (("kraken", false), ("krakenc", true))
        source = binary_source(; complex=complex)
        if source === :none
            push!(parts, "$label: unavailable")
        else
            push!(parts, "$label: $(binary_path(; complex=complex)) [$source]")
        end
    end
    return join(parts, ", ")
end

include("env_writer.jl")

# Plan tasks 3.3-3.5 add these:
# include("mod_reader.jl")
# include("runner.jl")
# include("compare.jl")

end # module
