# Kraken.jl

Kraken.jl is a normal-mode model for underwater acoustic propagation, written entirely in Julia. It
reimplements [Michael Porter's KRAKEN](https://oalib-acoustics.org/models-and-software/normal-modes/)
from the Acoustics Toolbox, with ideas from
[UnderwaterAcoustics.jl](https://github.com/org-arl/UnderwaterAcoustics.jl).

Given a range-independent waveguide, it computes the horizontal wavenumbers and mode shapes by a
finite-difference discretization of the depth-separated wave equation, then sums the modes into the
acoustic pressure field. It follows KRAKEN's numerics closely: Sturm-sequence bracketing, inverse
iteration, and Richardson extrapolation over successively finer meshes. It is checked against
unmodified Fortran KRAKEN on every push.

What it adds is **differentiability**. The whole solve, from the environment through the modes to the
pressure field, can be differentiated in forward or reverse mode. So group speeds, sensitivities, and
gradient-based inversion come from the same code that computes the field.

## Installation

```julia
using Pkg
Pkg.add("Kraken")
```

## Quick start

```@example index
using Kraken

env = UnderwaterEnv(pekeris_env()...)       # 100 m of water over a fluid half-space
sol = kraken_jl(env, 100.0)                  # modes at 100 Hz
tl = transmission_loss(sol, 1_000.0:1_000.0:5_000.0, 36.0, 55.0)   # dB, source at 36 m, receiver at 55 m
```

## What it models

- **Layered fluid waveguides**: any number of media, each with its own sound-speed, density and
  attenuation profile, over a fluid half-space or a perfect boundary.
- **Boundary conditions**: pressure-release, rigid, or acoustic half-space, as in KRAKEN.
- **Sound-speed interpolation**: C-linear, n²-linear, or KRAKEN's not-a-knot cubic spline.
- **Attenuation** in all six of KRAKEN's unit conventions, giving complex wavenumbers.
- **Fields**: coherent and incoherent transmission loss, and broadband pulse synthesis.
- **Automatic differentiation**: ForwardDiff for a few parameters, Zygote or Mooncake for many.

It does not model leaky modes, elastic media, or range dependence. The
[Validation](@ref validation) page lists the gaps against KRAKEN alongside the measured agreement.

## Where to go next

- [Getting started](@ref) — from an environment to modes, a field, a pulse, and group speeds.
- [Environments](@ref environments) — the matrix layout and every option an environment takes.
- [Automatic differentiation](@ref) — forward and reverse mode, what each costs, and two worked
  inversions.
- [Validation](@ref validation) — agreement with Fortran KRAKEN, computed when these docs are built.
- [API reference](@ref) — every exported function and type.
