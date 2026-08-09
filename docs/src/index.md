# Kraken.jl

!!! warning "Documentation under construction"
    Most of the guides and tutorials arrive in Milestone 8 of the revival plan; what is here now is
    the README, the [Automatic differentiation](@ref) guide, and a generated API reference.

Kraken.jl is a normal-mode simulation package for underwater acoustic propagation, written entirely
in Julia. It is based on [Michael Porter's KRAKEN Fortran code from the Acoustics
Toolbox](https://oalib-acoustics.org/models-and-software/normal-modes/) and on
[UnderwaterAcoustics.jl](https://github.com/org-arl/UnderwaterAcoustics.jl).

It computes horizontal wavenumbers and mode shapes for a range-independent environment by
finite-difference discretization of the depth-separated wave equation, then uses those modes to
synthesize the acoustic pressure field. Because the solver is written to be differentiable, you can
take derivatives of the results with respect to any environment parameter — group speeds, for
instance, are a derivative with respect to frequency. Both forward and reverse mode work; see
[Automatic differentiation](@ref).

## Installation

```julia
using Pkg
Pkg.add("Kraken")
```

## Usage

```julia
using Kraken

# Load the environment (same structure as Acoustics Toolbox .env files)
ssp, layers, sspHS = pekeris_env()
env = UnderwaterEnv(ssp, layers, sspHS)

# Run the simulation
freq = 100.0
sol = kraken_jl(env, freq)

# Access the results
modes = sol.modes
wavenumbers = sol.kr
zn = vcat(sol.props.zn_vec...)
```

## Calculating group speeds

Group speed is the derivative of angular frequency ``\omega = 2\pi f`` with respect to the
horizontal wavenumber ``k_{r,m}``:

```math
c_g = \frac{\partial \omega}{\partial k_{r,m}}
```

That is a derivative with respect to one parameter, so this is a direct application of
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). For gradients with respect to *many*
parameters — a whole sound-speed profile, say — use reverse mode instead; see
[Automatic differentiation](@ref).

```julia
using ForwardDiff
using Kraken

function calculate_kr_pekeris(freq)
    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    return find_kr(env, props, cache)
end

freq = 100.0
group_speeds = 2pi ./ ForwardDiff.derivative(calculate_kr_pekeris, freq)
```

## Plotting

Plotting is a package extension, so it costs nothing unless you ask for it. Load any Makie backend
and [`plot_modes`](@ref) and [`plot_ssp`](@ref) become available:

```julia
using Kraken, CairoMakie   # or GLMakie

env = UnderwaterEnv(pekeris_env()...)
sol = kraken_jl(env, 100.0)
plot_modes(sol; modes=1:5)
```

## Validation against Fortran KRAKEN

The solver is checked against **unmodified** Fortran KRAKEN on every push. The reference binaries
come from [`AcousticsToolbox_jll`](https://github.com/JuliaBinaryWrappers/AcousticsToolbox_jll.jl),
so CI needs no Fortran toolchain, and the comparison runs over the toolbox's own `.env`/`.mod` file
interface. Kraken.jl itself links against no Fortran code and ships no shared library — the harness
is test-only, under `test/reference/`, and is not part of the package's public surface.

Across the five standard environments at 25–400 Hz the largest relative wavenumber difference is
2.7e-5 and the smallest mode-shape correlation is 0.99995. Of the 402 environment files shipped with
the Acoustics Toolbox, 65 use only features Kraken.jl models today, and all 65 agree; the rest are
reported with the specific feature that blocks them. `test/README.md` in the repository has the
per-environment table.

To compare against a particular Acoustics Toolbox build instead of the packaged one, point
`KRAKEN_FORTRAN_BIN` at a directory containing `kraken.exe`.

!!! note "KrakenFortran.jl is a different thing"
    [KrakenFortran.jl](https://github.com/vardister/KrakenFortran.jl) is a separate, optional package
    that calls Fortran KRAKEN *in process* through `ccall`. Kraken.jl does not depend on it and does
    not use it for validation: its sources are a MEX-adapted fork of an older KRAKEN, so it is a
    performance option for broadband sweeps rather than a statement about correctness. Validating
    against a fork would prove nothing about agreement with upstream KRAKEN.

## Where to go next

* [Automatic differentiation](@ref) — forward and reverse mode, the cost of each, and a worked
  gradient-based sound-speed inversion.
* [API reference](@ref) — the full list of exported functions and types.
