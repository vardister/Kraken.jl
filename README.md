# Kraken.jl

[![CI](https://github.com/vardister/Kraken.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/vardister/Kraken.jl/actions/workflows/CI.yml)
[![Format](https://github.com/vardister/Kraken.jl/actions/workflows/Format.yml/badge.svg)](https://github.com/vardister/Kraken.jl/actions/workflows/Format.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://vardister.github.io/Kraken.jl/dev/)

Kraken.jl computes underwater sound propagation with normal modes. It is a pure-Julia
reimplementation of Michael Porter's [KRAKEN](https://oalib-acoustics.org/models-and-software/normal-modes/)
from the Acoustics Toolbox, with ideas from
[UnderwaterAcoustics.jl](https://github.com/org-arl/UnderwaterAcoustics.jl). Given a range-independent
ocean (a sound-speed profile, sediment layers and a seabed), it finds the horizontal wavenumbers and
mode shapes, and from those the pressure field and transmission loss.

Every step is differentiable, in forward mode and in reverse mode. That makes group speeds,
sensitivities and gradient-based inversions a matter of calling an AD package on the solver.

- **Validated against Fortran KRAKEN.** Every push compares the solver with the *unmodified*
  `kraken.exe` and `field.exe` from
  [`AcousticsToolbox_jll`](https://github.com/JuliaBinaryWrappers/AcousticsToolbox_jll.jl), on
  standard environments, lossy ones, every boundary condition and every profile interpolation.
- **What it models:** layered fluid media, compressional attenuation in six unit conventions,
  pressure-release, rigid and acoustic half-space boundaries, three sound-speed interpolations,
  coherent and incoherent transmission loss, and broadband pulses.
- **Differentiable:** forward mode through [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl),
  reverse mode through [Zygote.jl](https://github.com/FluxML/Zygote.jl) or
  [Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl). Both use implicit-function rules for the
  root finder and the eigenvectors, so neither differentiates solver iterations.
- **Optional plotting:** `using CairoMakie` (or `GLMakie`) enables `plot_modes`, `plot_ssp`,
  `plot_tl`, `plot_tl_range`, `plot_dispersion` and `plot_arrivals`.

The [documentation](https://vardister.github.io/Kraken.jl/dev/) has a getting-started guide, the
environment format, the AD details, and the validation results measured by the docs build itself.

## Installation

```julia-repl
pkg> add Kraken
```

## Quick start

An environment is three matrices laid out like an Acoustics Toolbox `.env` file. `pekeris_env()`
returns them for a 100 m isovelocity water column over a fast seabed:

```julia
using Kraken

ssp, layers, sspHS = pekeris_env()
env = UnderwaterEnv(ssp, layers, sspHS)

sol = kraken_jl(env, 100.0)   # solve at 100 Hz
sol.kr                        # horizontal wavenumbers (1/m), one per mode
sol.modes                     # mode shapes, one column per mode
```

`kraken_jl` refines the depth mesh and Richardson-extrapolates the wavenumbers until they converge.
From the modes, the field at any ranges and receiver depths is a mode sum:

```julia
ranges = 100.0:100.0:10_000.0   # m
zs, zr = 25.0, [10.0, 50.0]     # source depth, receiver depths (m)

p = acoustic_field(sol, ranges, zs, zr)        # complex pressure, length(zr) × length(ranges)
tl = transmission_loss(sol, ranges, zs, zr)    # dB re 1 m
tl_smooth = transmission_loss(sol, ranges, zs, zr; mode=:incoherent)
```

With a Makie backend loaded, `plot_tl(sol, ranges, zs, 0.0:1.0:100.0)` draws the TL as a
range-depth image and `plot_modes(sol)` draws the modes beside the profile.

## Environments

Each row of `ssp` is `[depth, cp, cs, ρ, αp, αs]`: depth in m (positive downwards), speeds in m/s,
density in **kg/m³** (a `.env` file uses g/cm³), and attenuation in the environment's units. A layer
interface is two rows at the same depth, one for each side. `layers` has one row per medium, with the
medium's lower interface depth in the third column. `sspHS` holds the top and bottom half-spaces.

Here is a 74 m shelf with a soft, slow sediment layer over a hard, lossy basement:

```julia
shelf_ssp = [
    0.0 1468.3 0.0 1030.0 0.0 0.0
    74.4 1468.3 0.0 1030.0 0.0 0.0
    74.4 1445.0 0.0 1600.0 0.05 0.0   # sediment: slower than the water above it
    85.4 1555.0 0.0 1600.0 0.05 0.0
]
shelf_layers = [0.0 0.0 74.4; 0.0 0.0 85.4]
shelf_hs = [0.0 343.0 0.0 1.21 0.0 0.0; 85.4 1830.0 0.0 2000.0 0.25 0.0]

shelf = UnderwaterEnv(shelf_ssp, shelf_layers, shelf_hs)
kraken_jl(shelf, 50.0).kr   # complex: attenuation makes every mode decay with range
```

The keyword arguments cover the rest of what a `.env` file can say:

```julia
UnderwaterEnv(
    shelf_ssp, shelf_layers, shelf_hs;
    atten_units=:dB_per_m,        # or :nepers_per_m, :dB_per_kmHz, :dB_per_wavelength (default), :Q, :loss_parameter
    ssp_interp=:n2_linear,        # or :c_linear (default), :cubic_spline
    top_bc=PressureRelease(),     # or RigidBoundary()
    bottom_bc=AcousticHalfspace(),# or RigidBoundary(), PressureRelease()
)
```

A lossless environment returns real `Float64` wavenumbers. Loss is added by first-order perturbation,
which is what `kraken.exe` does. The full complex solve of `krakenc.exe`, which also finds leaky
modes, is not implemented.

## Automatic differentiation

Kraken.jl depends only on ChainRulesCore for its derivative rules. Add ForwardDiff, Zygote or Mooncake
to your own environment.

**Forward mode** costs one solve per input, so it suits a handful of parameters. Group speed is
`dω/dkᵣ`, a derivative with respect to frequency alone:

```julia
using ForwardDiff

dkr_df = ForwardDiff.derivative(f -> kraken_jl(env, f).kr, 100.0)
group_speeds = 2π ./ dkr_df   # m/s, one per mode
```

**Reverse mode** costs a fixed multiple of one solve however many inputs there are, so it suits a
whole sound-speed profile. Build the depth grid outside the differentiated function:

```julia
using Zygote

z = collect(range(0.0, 100.0, 50))
c = fill(1500.0, 50)

function modal_sum(cvec)
    ssp = hcat(z, cvec, zero(cvec), fill(1000.0, 50), zero(cvec), zero(cvec))
    sspHS = [0.0 343.0 0.0 1.21 0.0 0.0; 100.0 1600.0 0.0 1500.0 0.0 0.0]
    env = UnderwaterEnv(ssp, [0.0 0.0 100.0], sspHS)
    return sum(kraken_jl(env, 100.0).kr)
end

grad = Zygote.gradient(modal_sum, c)[1]   # all 50 derivatives in one pass
```

In the package's benchmark, reverse mode is 4.4× faster than forward mode at 50 profile points and
47× faster at 500. The [AD page](https://vardister.github.io/Kraken.jl/dev/ad/) has the
measurements, a worked inversion that recovers a sound-speed profile, and the main caveat: the
top-level derivative is conditional on the depth mesh, which is chosen from the frequency and
changes in integer steps.

## Validation

Kraken.jl is checked against unmodified Fortran KRAKEN on every push. The binaries come from
`AcousticsToolbox_jll`, so CI needs no Fortran compiler, and the comparison goes through the toolbox's
own `.env`/`.mod` files: nothing here links against Fortran. The test suite compares wavenumbers, mode
shapes, modal attenuation, group speeds, AD gradients and the `field.exe` pressure field. Set
`KRAKEN_OALIB_TESTS` to a toolbox `tests/` directory and it also runs the toolbox's own environment
files. [`test/README.md`](test/README.md) has the tolerances and what each comparison measured.

## Examples

The scripts in [`examples/`](examples/) have their own environment:

```sh
$ julia --project=examples -e 'using Pkg; Pkg.instantiate()'
$ julia --project=examples examples/dispersion.jl
```

| Script | What it shows |
|---|---|
| [`deep_water_modes.jl`](examples/deep_water_modes.jl) | A measured 5 km profile over a sediment layer and a rigid basement, and its modes |
| [`dispersion.jl`](examples/dispersion.jl) | Phase speed, group speed (by ForwardDiff) and attenuation against frequency on a lossy shelf |
| [`broadband_sweep.jl`](examples/broadband_sweep.jl) | A frequency sweep across threads, checked against the serial result |
| [`reverse_ad.jl`](examples/reverse_ad.jl) | A tour of reverse-mode AD, with each gradient checked against ForwardDiff and finite differences |

## Limitations

- No shear waves: `cs` and `αs` are read and ignored.
- No leaky modes: attenuation is perturbational, like `kraken.exe`, not the complex solve of `krakenc.exe`.
- Range-independent only.
- `:cubic_spline` interpolation is limited to a single medium, because the spline runs through the
  whole column where KRAKEN splines each medium separately.

## Contributing

Open an issue before starting a major change. `test/README.md` explains
how to run the suite, including against a local Acoustics Toolbox build.
