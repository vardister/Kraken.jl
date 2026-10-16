### Plotting
#
# Implemented in ext/KrakenMakieExt.jl, which loads as soon as any Makie backend is present
# (`using CairoMakie` or `using GLMakie` both pull in Makie). The stubs exist here so the names are
# exported, so the docstrings live with the package rather than the extension (the extension's
# methods carry none, or every docstring would be listed twice once a backend is loaded), and so
# calling one without a backend gives an actionable message instead of a bare `MethodError`.
#
# Each figure-level function returns a `Makie.Figure`. The ones ending in `!` draw into an axis you
# already have and take plain arrays, which is what lets output from somewhere else — `field.exe`,
# say — be drawn with exactly the same styling as Kraken.jl's own.

export plot_modes, plot_ssp, plot_ssp!, plot_tl, plot_tl!, plot_tl_range, plot_tl_range!
export plot_dispersion, plot_arrivals

"""
    plot_modes(sol::NormalModeSolution; modes=1:min(10, length(sol.kr)))

Plot the sound-speed profile beside the mode shapes of `sol`, both against depth with the surface
at the top.

Requires a Makie backend: run `using CairoMakie` (or `using GLMakie`) to enable it. Without one this
throws an error saying so.

# Arguments
- `sol`: the solution from [`kraken_jl`](@ref).
- `modes`: which modes to draw, by index.

# Example
```julia
using CairoMakie
sol = kraken_jl(UnderwaterEnv(pekeris_env()...), 100.0)
fig = plot_modes(sol; modes=1:3)
```
"""
function plot_modes end

"""
    plot_ssp(env::UnderwaterEnv; npoints=500)
    plot_ssp(envs::AbstractVector{<:Pair}; npoints=500)

Plot the sound-speed profile of an environment against depth, from the surface to the bottom
boundary. Given `label => env` pairs, overlays their profiles with a legend.

Requires a Makie backend: run `using CairoMakie` (or `using GLMakie`) to enable it.

# Arguments
- `env`: the [`UnderwaterEnv`](@ref) whose profile to draw.
- `envs`: a vector of `label => env` pairs.
- `npoints`: how many depths to evaluate each profile at.

# Example
```julia
using CairoMakie
fig = plot_ssp(UnderwaterEnv(munk_env()...))
```
"""
function plot_ssp end

"""
    plot_ssp!(ax, env::UnderwaterEnv; npoints=500, kwargs...)

Draw the sound-speed profile of `env` into the Makie axis `ax`, as sound speed against depth.
`kwargs` go to `Makie.lines!` (`label`, `color`, `linestyle`, …). Returns the plot object.

The axis is not configured; [`plot_ssp`](@ref) makes one with reversed depth and labels.

# Example
```julia
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1]; yreversed=true)
plot_ssp!(ax, UnderwaterEnv(munk_env()...); label="C-linear")
plot_ssp!(ax, UnderwaterEnv(munk_env()...; ssp_interp=:n2_linear); label="n²-linear")
```
"""
function plot_ssp! end

"""
    plot_tl(sol::NormalModeSolution, ranges, zs, zr; mode=:coherent, nmodes=length(sol.kr),
            colorrange=nothing, colormap=Reverse(:viridis))

Plot transmission loss over a grid of ranges and receiver depths, as a heatmap with the surface at
the top and a colorbar in dB. The source position is marked.

Requires a Makie backend: run `using CairoMakie` (or `using GLMakie`) to enable it.

# Arguments
- `sol`: the solution from [`kraken_jl`](@ref).
- `ranges`: the ranges (m). They are drawn in km.
- `zs`: source depth (m).
- `zr`: receiver depths (m).
- `mode`, `nmodes`: passed to [`transmission_loss`](@ref).
- `colorrange`: the dB span of the color scale. The default runs from the 2nd percentile of the
  field to 50 dB above it, which keeps a deep null from washing out the rest.
- `colormap`: any Makie colormap. The default makes loud (low TL) bright.

# Example
```julia
using CairoMakie
sol = kraken_jl(UnderwaterEnv(munk_env()...), 50.0)
fig = plot_tl(sol, 1_000.0:250.0:100_000.0, 1000.0, 0.0:25.0:5000.0)
```
"""
function plot_tl end

"""
    plot_tl!(ax, ranges, zr, tl::AbstractMatrix; colorrange=nothing, colormap=Reverse(:viridis), kwargs...)

Draw a transmission-loss matrix into the Makie axis `ax` as a heatmap, with ranges in km along the
horizontal axis. `tl` has one row per receiver depth in `zr` and one column per range in `ranges`,
which is the layout [`transmission_loss`](@ref) returns. Non-finite values (the `Inf` of a
pressure-release surface) are left blank. `kwargs` go to `Makie.heatmap!`. Returns the plot object,
so a `Colorbar` can be attached to it.

Because it takes a plain matrix, it draws TL from any source with the same styling, and a
difference of two TL matrices too (pass a diverging `colormap` and a symmetric `colorrange`).

# Example
```julia
using CairoMakie
sol = kraken_jl(UnderwaterEnv(pekeris_env()...), 100.0)
r, z = 1_000.0:50.0:10_000.0, 0.0:1.0:100.0
fig = Figure()
ax = Axis(fig[1, 1]; yreversed=true)
hm = plot_tl!(ax, r, z, transmission_loss(sol, r, 36.0, z))
Colorbar(fig[1, 2], hm; label="TL (dB)")
```
"""
function plot_tl! end

"""
    plot_tl_range(sol::NormalModeSolution, ranges, zs, zr::Real; nmodes=length(sol.kr))

Plot transmission loss against range at a single receiver depth, coherent and incoherent together.
The coherent curve shows the modal interference; the incoherent one is the mean level it oscillates
about. TL increases downwards, as is conventional.

Requires a Makie backend: run `using CairoMakie` (or `using GLMakie`) to enable it.

# Arguments
- `sol`: the solution from [`kraken_jl`](@ref).
- `ranges`: the ranges (m), drawn in km.
- `zs`, `zr`: source and receiver depths (m).
- `nmodes`: passed to [`transmission_loss`](@ref).

# Example
```julia
using CairoMakie
sol = kraken_jl(UnderwaterEnv(pekeris_env()...), 100.0)
fig = plot_tl_range(sol, 100.0:10.0:10_000.0, 36.0, 55.0)
```
"""
function plot_tl_range end

"""
    plot_tl_range!(ax, ranges, tl; kwargs...)

Draw transmission loss against range into the Makie axis `ax`, with ranges in km. `tl` is a vector,
or a one-row matrix as [`transmission_loss`](@ref) returns for a single receiver depth. `kwargs` go
to `Makie.lines!`. Returns the plot object.

# Example
```julia
using CairoMakie
sol = kraken_jl(UnderwaterEnv(pekeris_env()...), 100.0)
r = 100.0:10.0:10_000.0
fig = Figure()
ax = Axis(fig[1, 1]; yreversed=true)
plot_tl_range!(ax, r, transmission_loss(sol, r, 36.0, 55.0); label="coherent")
```
"""
function plot_tl_range! end

"""
    plot_dispersion(env::UnderwaterEnv, freqs; nmodes=10, kwargs...)

Plot the phase speed `ω/kᵣ` of the first `nmodes` modes against frequency. Each mode appears at its
cutoff frequency, entering at the bottom half-space speed, and falls towards the slowest sound speed
in the water as frequency rises. Both limits are drawn as reference lines.

Solves `env` at every frequency in `freqs`; `kwargs` go to [`kraken_jl`](@ref). Only real
wavenumbers are drawn: for a lossy environment the phase speed uses `Re(kᵣ)`.

Requires a Makie backend: run `using CairoMakie` (or `using GLMakie`) to enable it.

# Arguments
- `env`: the [`UnderwaterEnv`](@ref) to solve.
- `freqs`: the frequencies (Hz).
- `nmodes`: how many modes to draw, from the first.

# Example
```julia
using CairoMakie
fig = plot_dispersion(UnderwaterEnv(pekeris_env()...), 5.0:1.0:150.0; nmodes=6)
```
"""
function plot_dispersion end

"""
    plot_arrivals(env::UnderwaterEnv, freqs, ranges, zs, zr, times;
                  spectrum=nothing, reduction_speed=maximum(env.c.c), kwargs...)

Plot a record section: the broadband pulse received at depth `zr` at each of several ranges, one
trace per range, offset by range. It shows how a pulse spreads out as the modes, travelling at
different group speeds, separate.

Each trace is synthesized with [`broadband_field`](@ref) and [`synthesize_pulse`](@ref), in
*reduced time* `t - r / reduction_speed`, so an arrival travelling at `reduction_speed` lines up
vertically and the whole section fits in one short window. Every trace is normalized to its own
peak, so late, weak arrivals stay visible.

Requires a Makie backend: run `using CairoMakie` (or `using GLMakie`) to enable it.

# Arguments
- `env`: the [`UnderwaterEnv`](@ref) to solve.
- `freqs`: evenly spaced frequencies (Hz), as [`synthesize_pulse`](@ref) requires.
- `ranges`: the ranges (m) at which to draw a trace.
- `zs`, `zr`: source and receiver depths (m).
- `times`: the reduced times (s) to evaluate each trace at.
- `spectrum`: the source spectrum, passed to [`synthesize_pulse`](@ref). Flat by default.
- `reduction_speed`: the speed (m/s) that sets the reduced time.
- `kwargs...`: passed to [`broadband_field`](@ref), and through it to [`kraken_jl`](@ref).

# Example
```julia
using CairoMakie
env = UnderwaterEnv(pekeris_env()...)
fig = plot_arrivals(env, 2.0:2.0:200.0, 2_000.0:2_000.0:10_000.0, 36.0, 55.0, 0:0.001:0.25)
```
"""
function plot_arrivals end

for f in (
    :plot_modes,
    :plot_ssp,
    :plot_ssp!,
    :plot_tl,
    :plot_tl!,
    :plot_tl_range,
    :plot_tl_range!,
    :plot_dispersion,
    :plot_arrivals,
)
    # With the extension loaded, reaching this catch-all means the arguments matched none of its
    # methods (a vector where a single depth was wanted, say), and blaming a missing backend would
    # send the user looking in the wrong place. Only without the extension is that the problem.
    @eval function $f(args...; kwargs...)
        Base.get_extension(@__MODULE__, :KrakenMakieExt) === nothing || throw(MethodError($f, args))
        return error(
            $("`" * string(f) * "` needs a Makie backend, which is a weak dependency of Kraken.jl. ") *
            "Run `using CairoMakie` (or `using GLMakie`) first; that loads Makie and activates " *
            "Kraken's plotting extension.",
        )
    end
end
