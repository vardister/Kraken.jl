# Kraken.jl

!!! warning "Documentation under construction"
    This is scaffolding. The guides, tutorials, and worked examples arrive in Milestone 8 of the
    revival plan; what is here now is the README plus a generated API reference.

Kraken.jl is a normal-mode simulation package for underwater acoustic propagation, written entirely
in Julia. It is based on [Michael Porter's KRAKEN Fortran code from the Acoustics
Toolbox](https://oalib-acoustics.org/models-and-software/normal-modes/) and on
[UnderwaterAcoustics.jl](https://github.com/org-arl/UnderwaterAcoustics.jl).

It computes horizontal wavenumbers and mode shapes for a range-independent environment by
finite-difference discretization of the depth-separated wave equation, then uses those modes to
synthesize the acoustic pressure field. Because the solver is written to be differentiable, you can
take derivatives of the results with respect to any environment parameter — group speeds, for
instance, are a derivative with respect to frequency.

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

The solver is differentiable, so this is a direct application of
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl):

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

## Where to go next

See [API reference](@ref) for the full list of exported functions and types.
