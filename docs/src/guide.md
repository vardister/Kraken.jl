# Getting started

This guide goes from an environment to modes, to a pressure field, to a pulse, using the Pekeris
waveguide throughout: 100 m of water at 1500 m/s over a 1600 m/s seabed. Every block on this page is
run when the documentation is built, and the output shown, figures included, is what it produced.

The figures come from Kraken.jl's plotting functions, which load with any Makie backend:

```@example guide
using Kraken
using CairoMakie   # activates Kraken's plotting extension
nothing # hide
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

[`plot_modes`](@ref) draws them beside the sound-speed profile. Below the dashed seafloor each mode
decays exponentially into the seabed, and the higher modes, nearer their cutoff, reach further:

```@example guide
plot_modes(sol; modes=1:5)
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

Over the whole water column the interference pattern is easier to see than to read off a table.
[`plot_tl`](@ref) computes the TL on a grid and draws it, with the source marked:

```@example guide
plot_tl(sol, 50.0:25.0:10_000.0, 36.0, 0.0:0.5:100.0)
```

[`plot_tl_range`](@ref) takes a single receiver depth and draws the coherent and incoherent sums
together. The incoherent curve is the mean level the interference oscillates about:

```@example guide
plot_tl_range(sol, 50.0:5.0:10_000.0, 36.0, 55.0)
```

The phase convention is KRAKEN's, `exp(i(ωt − kᵣr))`, and the amplitude is normalized to a source of
unit pressure at 1 m, so `field.exe`'s output can be compared with no rescaling. The
[Validation](@ref validation) page makes that comparison.

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

[`plot_arrivals`](@ref) does this at a series of ranges and stacks the results into a record
section. Here the source spectrum is a band around 60 Hz, set with `spectrum`. The pulse stretches
out as it travels, because its modes and frequencies travel at different group speeds, and by 50 km
the arrival is several times as long as it was at 5 km:

```@example guide
band(f) = exp(-((f - 60) / 20)^2)
plot_arrivals(env, 1.0:1.0:200.0, 5_000.0:5_000.0:50_000.0, 36.0, 55.0, 0:0.001:0.5; spectrum=band)
```

## Group speeds, by differentiation

The solver can be differentiated, so a derivative of its output is a derivative in the calculus
sense rather than a finite difference. Group speed is ``c_g = \partial\omega/\partial k_r``, which
is one derivative with respect to frequency. That is a job for forward mode:

```@example guide
using ForwardDiff

dkr_df = ForwardDiff.derivative(f -> kraken_jl(env, f).kr, 100.0)
group_speeds = 2π ./ dkr_df
```

The group speeds lie *below* the water's 1500 m/s: in this waveguide, energy travels more slowly than
sound in the water does.

The dispersion behind these speeds shows across frequency. [`plot_dispersion`](@ref) solves at each
frequency and draws each mode's phase speed. A mode appears at its cutoff, where its phase speed is
the seabed's 1600 m/s, and falls towards the water's 1500 m/s as frequency rises. The steeper the
curve, the more its group speed differs from its phase speed:

```@example guide
plot_dispersion(env, 5.0:1.0:150.0; nmodes=6)
```

For gradients with respect to many parameters at once, a whole sound-speed profile for example,
reverse mode is the tool. [Automatic differentiation](@ref) covers both modes and works through a
full inversion.

## Beyond Pekeris

The same calls work for any environment. Here is the Munk deep-water sound channel at 50 Hz, which
traps about a hundred modes:

```@example guide
munk = kraken_jl(UnderwaterEnv(munk_env()...), 50.0)
```

The low modes are trapped around the sound-channel axis at 1300 m and never touch the surface or
the bottom. The higher modes span more and more of the water column:

```@example guide
plot_modes(munk; modes=[1, 2, 10, 40])
```

The [Validation](@ref validation) page shows this environment's full transmission-loss field, with its
shadow zone and convergence zones, computed by both Kraken.jl and Fortran KRAKEN.

Layered seabeds, rigid or vacuum boundaries, other sound-speed interpolations, and attenuation are
all options of [`UnderwaterEnv`](@ref); see [Environments](@ref environments). The Pekeris waveguide
also has a closed-form solution, [`PekerisUnderwaterEnv`](@ref), which the finite-difference solver
reproduces:

```@example guide
closed_form = find_kr(PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0), 100.0)
maximum(abs.(closed_form .- sol.kr) ./ closed_form)
```

## Plotting

Plotting is a package extension, so it costs nothing until a Makie backend is loaded (`using
CairoMakie` for files, `using GLMakie` for interactive windows). The figure-level functions each
return a `Makie.Figure`:

| Function | Draws |
|---|---|
| [`plot_modes`](@ref) | the sound-speed profile beside the mode shapes |
| [`plot_ssp`](@ref) | one sound-speed profile, or several overlaid |
| [`plot_tl`](@ref) | transmission loss over range and depth |
| [`plot_tl_range`](@ref) | transmission loss against range at one depth |
| [`plot_dispersion`](@ref) | modal phase speed against frequency |
| [`plot_arrivals`](@ref) | broadband arrivals at a series of ranges |

[`plot_ssp!`](@ref), [`plot_tl!`](@ref) and [`plot_tl_range!`](@ref) draw into an axis you already
have and take plain arrays. Use them to build your own layouts, or to draw data from elsewhere with
the same styling; the [Validation](@ref validation) page uses them to put `field.exe`'s output beside
Kraken.jl's.

!!! note "`density` after `using CairoMakie`"
    Makie exports a `density` plot of its own, so with both packages loaded by `using`, a bare
    `density` is ambiguous and Julia refuses it. Write `Kraken.density(env.ρ, z)`, or load the backend
    with `import CairoMakie`, which activates the plotting extension just as well.
