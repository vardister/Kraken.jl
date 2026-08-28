module KrakenMakieExt

using Kraken
using Kraken: NormalModeSolution, soundspeed
using Makie

"""
    plot_modes(sol::NormalModeSolution; modes=1:min(10, length(sol.kr)))

Two-panel figure: the environment's sound-speed profile on the left, the requested mode shapes on
the right, both against depth with the surface at the top.

Available only once a Makie backend is loaded (`using CairoMakie` or `using GLMakie`).
"""
function Kraken.plot_modes(sol::NormalModeSolution; modes=1:min(10, length(sol.kr)))
    isempty(sol.kr) && throw(ArgumentError("solution has no modes to plot"))
    n_modes = size(sol.modes, 2)
    all(m -> 1 <= m <= n_modes, modes) ||
        throw(ArgumentError("requested modes $(modes) but the solution has $(n_modes); pass `modes=1:$(n_modes)`"))

    # The mode shapes are sampled on the solver's own mesh, which is the concatenation of the
    # per-layer depth vectors — the same order as the rows of `sol.modes`.
    zn = reduce(vcat, sol.props.zn_vec)

    fig = Figure()

    ax1 = Axis(fig[1, 1]; xlabel="Sound speed (m/s)", ylabel="Depth (m)", title="Sound speed profile", yreversed=true)
    lines!(ax1, soundspeed(sol.env.c, zn), zn)

    ax2 = Axis(fig[1, 2]; xlabel="Mode amplitude ϕ", title="Modes", yreversed=true)
    for m in modes
        lines!(ax2, sol.modes[:, m], zn; label="mode $(m)")
    end
    linkyaxes!(ax1, ax2)
    hideydecorations!(ax2; grid=false)
    axislegend(ax2; position=:rb)

    return fig
end

"""
    plot_ssp(env::UnderwaterEnv; npoints=500)

Sound-speed profile of `env`, sampled on a uniform grid from the surface to the bottom.
"""
function Kraken.plot_ssp(env::UnderwaterEnv; npoints=500)
    z = range(0, env.depth; length=npoints)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Sound speed (m/s)", ylabel="Depth (m)", title="Sound speed profile", yreversed=true)
    lines!(ax, soundspeed(env.c, collect(z)), collect(z))
    return fig
end

end # module
