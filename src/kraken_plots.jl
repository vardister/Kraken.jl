function is_headless()
    if Sys.isapple()
        # Check if we can access the WindowServer process
        try
            output = read(`ps -e`, String)
            return !occursin("WindowServer", output)
        catch
            # If we can't run the command, assume it's headless
            return true
        end
    else
        # For non-Mac systems, use the previous method
        has_display = haskey(ENV, "DISPLAY")
        return !has_display
    end
end


# Usage
# if is_headless()
#     using CairoMakie
# else
#     using GLMakie
# end
using CairoMakie


function plot_modes(res::Dict, ssp; modes = 1:10)
    f = Figure()

    ax1 = Axis(
        f[1, 1],
        xlabel = "Sound Speed (m/s)",
        ylabel = "Depth (m)",
        title = "Sound Speed Profile",
        yreversed = true
    )
    lines!(ax1, ssp[:, 2], ssp[:, 1])

    ax2 = Axis(f[1, 2], xlabel = "Mode Amplitude ϕ", title = "Modes")
    for (i, mode) in enumerate(eachcol(res["modes"][:, modes]))
        lines!(ax2, mode, vec(-res["zm"]), label = "mode $(modes[i])")
    end
    axislegend(ax2; position = :rb)
    hideydecorations!(ax2, grid = false)
    return f
end