# Getting started

This guide goes from an environment to modes, to a pressure field, to a pulse, using the Pekeris
waveguide throughout: 100 m of water at 1500 m/s over a 1600 m/s seabed. Every block on this page is
run when the documentation is built, and the output shown is what it printed.

```@setup guide
using Kraken
```

## Solve for the modes

An environment comes from three matrices in the Acoustics Toolbox `.env` layout (see
[Environments](@ref environments)). The standard environments return them ready-made:

```@example guide
using Kraken

env = UnderwaterEnv(pekeris_env()...)
sol = kraken_jl(env, 100.0)   # frequency in Hz
```

[`kraken_jl`](@ref) returns a `NormalModeSolution`. Its `kr` field holds the horizontal wavenumbers
in 1/m, largest first. The largest belongs to the mode that travels most nearly horizontally:

```@example guide
sol.kr
```

The phase speed of each mode, `ω/kᵣ`, lies between the water and seabed sound speeds. That is what
makes a mode *trapped*:

```@example guide
2π * 100.0 ./ sol.kr
```

## Mode shapes

`sol.modes` holds one column per mode, sampled on the solver's depth mesh. To evaluate the modes at
depths of your choosing, use [`mode_amplitudes`](@ref), which interpolates between mesh points and
continues each mode into the seabed as the decaying exponential it is there:

```@example guide
mode_amplitudes(sol, [0.0, 25.0, 50.0, 100.0, 120.0])
```

Every mode is zero at the pressure-release surface. Mode `m` crosses zero `m − 1` times in the water
column. The modes are normalized so that ``\int \psi_m^2/\rho \, dz = 1``, including the part of the
integral in the seabed.

## Pressure field and transmission loss

[`transmission_loss`](@ref) sums the modes into the field of a point source and returns it in dB,
with receiver depths down the rows and ranges across the columns. For a source at 36 m and one
receiver at 55 m:

```@example guide
ranges = 1_000.0:1_000.0:10_000.0
round.(transmission_loss(sol, ranges, 36.0, 55.0); digits=1)
```

The level does not fall smoothly with range, because the modes interfere. An *incoherent* sum adds
their intensities instead of their complex amplitudes. That discards the interference pattern and
leaves the mean trend, which is often what you want when the pattern is too fine to resolve:

```@example guide
round.(transmission_loss(sol, ranges, 36.0, 55.0; mode=:incoherent); digits=1)
```

[`acoustic_field`](@ref) gives the complex pressure itself, on any grid of ranges and depths:

```@example guide
p = acoustic_field(sol, 1_000.0:10.0:10_000.0, 36.0, 0.0:1.0:100.0)
size(p)
```

The phase convention is KRAKEN's, `exp(i(ωt − kᵣr))`, and the amplitude is normalized to a source of
unit pressure at 1 m. That is why `field.exe`'s output can be compared with no rescaling, and the
[Validation](@ref validation) page does exactly that.

## A pulse

A single frequency gives interference; a band of frequencies gives arrivals in time.
[`broadband_field`](@ref) solves the environment at each frequency, and [`synthesize_pulse`](@ref)
Fourier-synthesizes a time series from the result. `t_reduce` starts the time axis at the earliest
arrival any mode could make, so a short window can hold a pulse that has travelled 5 km:

```@example guide
freqs = 5.0:5.0:300.0
P = broadband_field(env, freqs, 5_000.0, 36.0, 55.0)
t = 0:0.001:0.6
pulse = synthesize_pulse(P[:, 1, 1], freqs, t; t_reduce=5_000.0 / 1500.0)
t[argmax(abs.(pulse))]   # seconds after r/c
```

Each frequency is an independent solve. A frequency that supports no trapped modes, or zero,
contributes nothing rather than an error, which is what a band that starts near the lowest mode's
cutoff needs.

## Group speeds, by differentiation

The solver can be differentiated, so a derivative of its output is a derivative in the calculus
sense rather than a finite difference. Group speed is ``c_g = \partial\omega/\partial k_r``, which
is one derivative with respect to frequency. That is a job for forward mode:

```@example guide
using ForwardDiff

dkr_df = ForwardDiff.derivative(f -> kraken_jl(env, f).kr, 100.0)
group_speeds = 2π ./ dkr_df
```

The group speeds lie *below* the water sound speed, as they should: in this waveguide, energy
travels more slowly than the sound speed of the water.

For gradients with respect to many parameters at once, a whole sound-speed profile for example,
reverse mode is the tool. [Automatic differentiation](@ref) covers both modes and works through a
full inversion.

## Beyond Pekeris

The same calls work for any environment. Here is the Munk deep-water sound channel at 50 Hz, which
traps about a hundred modes:

```@example guide
munk = kraken_jl(UnderwaterEnv(munk_env()...), 50.0)
```

Layered seabeds, rigid or vacuum boundaries, other sound-speed interpolations, and attenuation are
all options of [`UnderwaterEnv`](@ref); see [Environments](@ref environments). The Pekeris waveguide
also has a closed-form solution, [`PekerisUnderwaterEnv`](@ref), which the finite-difference solver
reproduces:

```@example guide
closed_form = find_kr(PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0), 100.0)
maximum(abs.(closed_form .- sol.kr) ./ closed_form)
```

## Plotting

Plotting is a package extension, so it costs nothing until you load a Makie backend:

```julia
using Kraken, CairoMakie   # or GLMakie

plot_modes(sol; modes=1:5)   # the sound-speed profile beside the mode shapes
plot_ssp(UnderwaterEnv(munk_env()...))
```

Both return a `Makie.Figure`.
