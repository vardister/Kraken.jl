module Kraken

### New Julia implementation
include("kraken_core.jl")
include("kraken_ad.jl")
include("kraken_pekeris.jl")

include("kraken_standard_environments.jl")

### Plotting
# Implemented in ext/KrakenMakieExt.jl, which loads as soon as any Makie backend is present
# (`using CairoMakie` or `using GLMakie` both pull in Makie). The stubs exist here so the names are
# exported and give an actionable message instead of a bare MethodError when no backend is loaded.
export plot_modes, plot_ssp

"""
    plot_modes(sol::NormalModeSolution; modes=1:min(10, length(sol.kr)))

Plot the sound-speed profile alongside the computed mode shapes.

Requires a Makie backend: run `using CairoMakie` (or `using GLMakie`) to enable it.
"""
function plot_modes end

"""
    plot_ssp(env::UnderwaterEnv; npoints=500)

Plot the sound-speed profile of an environment.

Requires a Makie backend: run `using CairoMakie` (or `using GLMakie`) to enable it.
"""
function plot_ssp end

for f in (:plot_modes, :plot_ssp)
    @eval function $f(args...; kwargs...)
        return error(
            $("`" * string(f) * "` needs a Makie backend, which is a weak dependency of Kraken.jl. ") *
            "Run `using CairoMakie` (or `using GLMakie`) first — that loads Makie and activates " *
            "Kraken's plotting extension.",
        )
    end
end

end
