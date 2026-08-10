# Kraken.jl

[![CI](https://github.com/vardister/Kraken.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/vardister/Kraken.jl/actions/workflows/CI.yml)
[![Format](https://github.com/vardister/Kraken.jl/actions/workflows/Format.yml/badge.svg)](https://github.com/vardister/Kraken.jl/actions/workflows/Format.yml)

**❗Documentation is currently under construction.**

KRAKEN.jl is a Normal-Mode based simulation package for underwater acoustic propagation. It's heavily based on [Michael Porter's KRAKEN Fortran code located in the Acoustics Toolbox](https://oalib-acoustics.org/models-and-software/normal-modes/) and [UnderwaterAcoustics.jl](https://github.com/org-arl/UnderwaterAcoustics.jl).

This reimplementation is fully written in Julia, and is designed to be more user-friendly, and easier to extend. It is also designed to be more efficient, and to take advantage of Julia's parallelization capabilities.

## Features

- Normal-Mode based simulation for underwater acoustic propagation fully written in Julia
- User-friendly and easy to extend
- Differentiable in **both** modes — forward via [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), reverse via [Zygote.jl](https://github.com/FluxML/Zygote.jl) or [Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl), so the gradient with respect to a whole sound-speed profile costs about one solve instead of one solve per point (see below)
- Re-using existing environmental data files from the Acoustics Toolbox — `.env` files are both read and written
- Cross-validated against the **unmodified** Fortran `kraken.exe` on every push (see below)
- Optional plotting via a package extension — `using CairoMakie` (or `GLMakie`) enables `plot_modes` and `plot_ssp` at no cost to anyone who does not

## Validation

Every push runs the solver against unmodified Fortran KRAKEN. The binaries come from
[`AcousticsToolbox_jll`](https://github.com/JuliaBinaryWrappers/AcousticsToolbox_jll.jl), so CI needs
no Fortran toolchain, and the comparison is driven over the toolbox's own `.env`/`.mod` file
interface — nothing in this package links against Fortran code or ships a shared library.

Across the five standard environments at 25–400 Hz, the largest relative wavenumber difference
against `kraken.exe` is **2.7e-5** and the smallest mode-shape correlation is **0.99995**. Of the 402
environment files shipped with the Acoustics Toolbox, **167** use only features Kraken.jl models
today — 101 of them lossy, which attenuation support unlocked. See [`test/README.md`](test/README.md)
for the per-environment table, the attenuation comparison, and what blocks the rest.

## Missing features
- [ ] Inclusion of shear wave properties in environment
- [ ] Boundary conditions other than a pressure-release surface over a fluid half-space
- [ ] SSP interpolation other than C-linear
- [ ] Leaky modes — the full complex solve, `krakenc.exe`'s job (see the caveat under Attenuation)
- [x] Compressional wave attenuation in environment
- [x] Reverse-mode automatic differentiation

## Attenuation

An environment's `αp` column is read and used, and a lossy waveguide gives **complex** wavenumbers
whose imaginary part is the modal attenuation:

```julia
env = UnderwaterEnv(pekeris_env(; αb=0.5)...)   # 0.5 dB per wavelength in the seabed
sol = kraken_jl(env, 100.0)
sol.kr            # ComplexF64; Im(kᵣ) < 0, so the mode decays with range
imag(sol.kr[1])   # nepers per metre of range for mode 1
```

A lossless environment still returns real `Float64` wavenumbers, unchanged and bit-identical to
before attenuation existed — the complex path is entered only when there is loss to model.

Units follow KRAKEN's own six conventions, selected with `atten_units` and defaulting to
dB/wavelength (`'W'`, KRAKEN's `TopOpt(3:3)`):

| `atten_units` | `.env` char | meaning |
|---|---|---|
| `:nepers_per_m` | `N` | nepers/m |
| `:dB_per_m` | `M` | dB/m |
| `:dB_per_kmHz` | `F` | dB/(m·kHz), i.e. dB/(km·Hz) |
| `:dB_per_wavelength` | `W` | dB/wavelength *(default)* |
| `:Q` | `Q` | quality factor |
| `:loss_parameter` | `L` | loss parameter |

```julia
UnderwaterEnv(pekeris_env(; α0=1e-3)...; atten_units=:nepers_per_m)
```

Four of the six depend on frequency, so an environment stores the value you gave it and converts at
solve time — exactly as the Fortran's `CRCI` does.

### The approximation, and where it runs out

Attenuation is added by **first-order perturbation** of the real eigenvalue, which is what
`kraken.exe` does; `krakenc.exe` is the separate program that solves the complex problem outright.
The two agree to first order in α *exactly* — shrink α and the disagreement falls as α² until it hits
the single precision of the reference file.

Two things limit how accurately `Im(kᵣ)` comes out, and they bite on opposite ends of the mode
spectrum:

- **The perturbation is first order in `v = 2ω α_b/(c_b γ²)`**, not in `α_b`, and `γ → 0` at the
  bottom cutoff — so a strongly attenuating seabed degrades the *least*-trapped mode first.
  Agreement with `kraken.exe` runs from 1.8e-3 on a weakly attenuating waveguide to ~10% for a
  near-cutoff mode over a 0.5 dB/λ seabed. Inherent to the method; the fix is the full complex solve.
- **Everything else is ordinary discretization**, and it is second order in the mesh spacing. Both the
  energy normalization and the perturbation integral are taken medium by medium, so a jump in `ρ` or
  `α` at a layer interface is resolved exactly rather than averaged across — on
  `one_layer_env(; α1=0.4)` that is the difference between 8.1e-2 and 2.9e-3 agreement, and it is
  what keeps a lossy sediment layer usable at the default mesh.

Neither affects `Re(kᵣ)`, which stays within 1e-4 of `kraken.exe` throughout.
[`test/README.md`](test/README.md) has the measured table.

## Installation

```julia
using Pkg
Pkg.add("Kraken")
```

## Usage

```julia
using Kraken

# Load the environment
ssp, layers, sspHS = pekeris_env() # Similar structure to environment files from the Acoustics Toolbox
env = UnderwaterEnv(ssp, layers, sspHS)

# Run the simulation
freq = 100.0
sol = kraken_jl(env, freq)

# Access the results
modes = sol.modes
wavenumbers = sol.kr
zn = vcat(sol.props.zn_vec...)
```

## Automatic differentiation

The solver is differentiable in **both** modes, and which one you want depends on how many parameters
you are differentiating with respect to. Forward mode costs one solve per parameter; reverse mode
costs a fixed multiple of one solve regardless of the parameter count. Kraken.jl itself depends only
on ChainRulesCore.jl — ForwardDiff, Zygote and Mooncake are yours to add.

### Group speeds (forward mode — one parameter)

Group speeds are defined as the derivative of the angular frequency $\omega = 2\pi f$ with respect to the wavenumbers $k_{r,m}$. This is written as

$$ c_g = \frac{\partial \omega}{\partial k_{r,m}} $$

That is a derivative with respect to a single parameter, so forward mode is the right tool:

```julia
using ForwardDiff
using Kraken


function calculate_kr_pekeris(freq)
    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    wavenumbers = find_kr(env, props, cache)
    return wavenumbers
end

freq = 100.0
group_speeds = 2pi ./ ForwardDiff.derivative(calculate_kr_pekeris, freq)
```

These are checked against Fortran KRAKEN's own group-speed table on every push, to within 0.1%.

### Whole-profile gradients (reverse mode — many parameters)

For an inversion the unknown is not one number but a whole sound-speed profile, and reverse mode
returns every derivative in one pass:

```julia
using Kraken, Zygote

z = collect(range(0.0, 100.0, 50))     # build the depth grid OUTSIDE the differentiated function
c = fill(1500.0, 50)

function modal_sum(cvec)
    ssp = hcat(z, cvec, zero(cvec), fill(1000.0, 50), zero(cvec), zero(cvec))
    sspHS = [0.0 343.0 0.0 0.00121 0.0 0.0; 100.0 1600.0 0.0 1500.0 0.0 0.0]
    env = UnderwaterEnv(ssp, [0.0 0.0 100.0], sspHS)
    return sum(kraken_jl(env, 100.0).kr)
end

grad = Zygote.gradient(modal_sum, c)[1]   # all 50 derivatives, one pass
```

Measured on a 2021 M1 — the sum of all modal wavenumbers of that waveguide at 100 Hz, differentiated
with respect to an `M`-point profile. (The measurement uses the single-mesh
`bisection`/`solve_for_kr` path rather than `kraken_jl`, so the ratios isolate the cost of
differentiation from the mesh-refinement loop.)

| `M` | forward / primal | reverse / primal |
|---:|---:|---:|
| 1 | 1.1× | 5.3× |
| 5 | 2.6× | 6.5× |
| 10 | 5.2× | 6.0× |
| 50 | 25.7× | 5.8× |
| 500 | 253× | 5.4× |

Forward mode grows linearly with `M`; reverse mode stays flat. The two cross near a dozen parameters,
and at a 500-point profile reverse mode is 47× faster.
The residual 5.8× is mostly Zygote's fixed ~0.3 ms tape overhead on a very small solve — at 400 Hz,
where the same waveguide carries 18 modes, a 50-parameter gradient costs 1.9× the primal.

Reverse mode goes through wavenumbers, mode shapes, and the top-level `kraken_jl`, and works with
Zygote and Mooncake alike — the rules are written once with ChainRulesCore.

The [documentation's AD page](https://vardister.github.io/Kraken.jl/dev/ad/) has the details,
including a worked gradient-based inversion that recovers a sound-speed profile from synthetic
wavenumbers, and the caveats (the top-level gradient is conditional on the mesh schedule).

## More examples

[`examples/reverse_ad.jl`](examples/reverse_ad.jl) is a runnable tour of the AD support:
`julia --project=test examples/reverse_ad.jl`.

❗The other files in `examples/` do not currently run: they call an `EnvKRAKEN` API that was removed
when the Fortran sources moved out of this repository. They are being rewritten against the current
API.

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes. For major changes, please open an issue first to discuss what you would like to change.
