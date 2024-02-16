using CairoMakie

function plot_modes(res::Dict; modes = 1:10)
    f = Figure()

    ax1 = Axis(
        f[1, 1],
        xlabel = "Sound Speed (m/s)",
        ylabel = "Depth (m)",
        title = "Sound Speed Profile",
    )
    lines!(ax1, ssp[:, 2], -ssp[:, 1])

    ax2 = Axis(f[1, 2], xlabel = "Mode Amplitude ϕ", title = "Modes")
    for (i, mode) in enumerate(eachcol(res["modes"][:, modes]))
        lines!(ax2, mode, vec(-res["zm"]), label = "mode $(modes[i])")
    end
    axislegend(ax2; position = :rb)
    return f
end