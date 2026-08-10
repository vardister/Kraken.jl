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

## Attenuation

The `αp` column of an environment is read and used. A lossy waveguide gives **complex** wavenumbers,
whose imaginary part is the modal attenuation in nepers per metre of range:

```julia
env = UnderwaterEnv(pekeris_env(; αb=0.5)...)   # 0.5 dB per wavelength in the seabed
sol = kraken_jl(env, 100.0)
sol.kr             # ComplexF64
imag(sol.kr[1])    # < 0: the mode loses energy with range
```

A **lossless** environment still returns real `Float64` wavenumbers, bit-identical to what it
returned before attenuation was supported — the complex path is entered only when there is loss to
model, so nothing downstream of a lossless solve changes type.

### Units

Attenuation units follow KRAKEN's own six conventions, chosen with `atten_units` and defaulting to
dB per wavelength. The characters are those of `TopOpt(3:3)` in a `.env` file, and a file read with
[`read_env_file`](https://github.com/vardister/Kraken.jl/blob/master/test/reference/env_reader.jl)
carries its own declaration through.

| `atten_units` | `.env` char | meaning | conversion to nepers/m |
|---|---|---|---|
| `:nepers_per_m` | `N` | nepers/m | ``\alpha`` |
| `:dB_per_m` | `M` | dB/m | ``\alpha / 8.6858896`` |
| `:dB_per_kmHz` | `F` | dB/(m·kHz) ≡ dB/(km·Hz) | ``\alpha f / 8685.8896`` |
| `:dB_per_wavelength` | `W` | dB/wavelength *(default)* | ``\alpha f / (8.6858896\,c)`` |
| `:Q` | `Q` | quality factor | ``\omega / (2 c \alpha)`` |
| `:loss_parameter` | `L` | loss parameter | ``\alpha \omega / c`` |

```julia
UnderwaterEnv(pekeris_env(; α0=1e-3)...; atten_units=:nepers_per_m)
```

Four of the six depend on frequency, which is why an environment stores the value you gave it and
converts only once a frequency is known — the same arrangement as the Fortran's `CRCI`.

### How it is computed, and where the approximation runs out

Kraken.jl follows `kraken.exe`: it solves the **real** eigenproblem and adds the loss as a
first-order perturbation of the eigenvalue,

```math
\delta(k_r^2) = -2i\omega \int \frac{a(z)\,\psi(z)^2}{c(z)\,\rho(z)}\,\mathrm{d}z
              \;-\; \frac{i\,\omega\,a_b\,\psi(D)^2}{\gamma\,c_b\,\rho_b},
\qquad k_r = \sqrt{k_r^2 + \delta(k_r^2)}
```

with ``\gamma = \sqrt{k_r^2 - (\omega/c_b)^2}`` and ``a`` the attenuation in nepers/m. `krakenc.exe`
is the separate program that solves the complex problem outright; that path is not implemented here.

Two separate things limit how accurately ``\mathrm{Im}(k_r)`` comes out, and they bite on opposite
modes — worth knowing which one you are in.

**1. The perturbation is first order, and in ``v = 2\omega a_b/(c_b\gamma^2)``, not in ``a_b``.**
Since ``\gamma \to 0`` at the bottom cutoff, `v` grows without bound for the *least*-trapped mode, so
a strongly attenuating half-space degrades the top of the mode spectrum first. Agreement with
`kraken.exe` runs from 1.8e-3 on a weakly attenuating waveguide to ~10% for a near-cutoff mode over a
0.5 dB/λ seabed. This one is inherent to the method — the fix is the full complex solve.

**2. Everything else is ordinary discretization, and it is second order.** Both the energy
normalization and the perturbation integral are taken *medium by medium*, matching what
`kraken.f90`'s `Normalize` does, so a jump in ``\rho`` or ``\alpha`` at a layer interface is resolved
exactly instead of averaged across by a trapezoid that straddles it. That matters more than it
sounds: a single straddled interval per interface is enough to drop the whole quadrature to first
order, and on `one_layer_env(; α1=0.4)` it was the difference between 8.1e-2 and 2.9e-3 agreement
with `kraken.exe`. The observed convergence order on that case is 2.00.

Worth knowing when calibrating against Fortran: at a lossy sediment layer `kraken.exe` itself carries
about 1.6e-3 of discretization error on its automatic mesh, so agreement below that says more about
its mesh than about ours.

Neither affects ``\mathrm{Re}(k_r)``, which stays within 1e-4 of `kraken.exe` on all of these. The
measured table is in
[`test/README.md`](https://github.com/vardister/Kraken.jl/blob/master/test/README.md).

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
