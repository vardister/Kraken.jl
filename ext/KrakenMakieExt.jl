module KrakenMakieExt

# The docstrings for everything here are on the stubs in src/kraken_plotting.jl. The methods below
# deliberately carry none: a docstring here would be attached to the same `Kraken.plot_*` binding and
# listed a second time wherever a backend is loaded, the documentation build included.

using Kraken
using Kraken: NormalModeSolution, soundspeed
using Makie
using Statistics: quantile

const TL_COLORMAP = Reverse(:viridis)

# Loud (low TL) is bright. The span starts at the 2nd percentile rather than the minimum so that the
# few metres around the source, which are tens of dB louder than anything else, do not compress the
# rest of the field into one color, and it is 50 dB wide because that is about what a mode sum
# resolves before the interference nulls take over.
function default_tl_range(tl)
    finite = filter(isfinite, tl)
    isempty(finite) && return (0.0, 1.0)
    lo = quantile(finite, 0.02)
    return (lo, lo + 50)
end

# `Inf` is what a pressure-release surface gives (p = 0 there); NaN is what Makie leaves blank.
finite_or_nan(x) = isfinite(x) ? Float64(x) : NaN

### Mode shapes and profiles

function Kraken.plot_modes(sol::NormalModeSolution; modes=1:min(10, length(sol.kr)))
    isempty(sol.kr) && throw(ArgumentError("solution has no modes to plot"))
    n_modes = size(sol.modes, 2)
    all(m -> 1 <= m <= n_modes, modes) ||
        throw(ArgumentError("requested modes $(modes) but the solution has $(n_modes); pass `modes=1:$(n_modes)`"))

    env = sol.env
    # `mode_amplitudes` rather than the raw `sol.modes`: the solver's mesh leaves out the boundary
    # samples a pressure-release surface fixes at zero, and only `mode_amplitudes` continues each mode
    # into a half-space bottom, which is a fifth of the picture for a weakly trapped mode.
    halfspace = env.bottom_bc isa AcousticHalfspace
    zmax = halfspace ? 1.2 * env.depth : env.depth
    z = collect(range(0, zmax; length=600))
    Φ = mode_amplitudes(sol, z; nmodes=maximum(modes))

    fig = Figure(; size=(900, 500))

    ax1 = Axis(fig[1, 1]; xlabel="Sound speed (m/s)", ylabel="Depth (m)", title="Sound speed profile", yreversed=true)
    zw = filter(<=(env.depth), z)
    cw = soundspeed(env.c, zw)
    ssp_color = Makie.wong_colors()[1]
    lines!(ax1, cw, zw; color=ssp_color)
    halfspace && lines!(ax1, [env.cb, env.cb], [env.depth, zmax]; color=ssp_color)
    # An isovelocity column would otherwise be autoscaled to a span of thousands of m/s.
    c_all = halfspace ? vcat(cw, env.cb) : cw
    lo, hi = extrema(c_all)
    pad = max(0.1 * (hi - lo), 10.0)
    xlims!(ax1, lo - pad, hi + pad)

    ax2 = Axis(fig[1, 2]; xlabel="Mode amplitude ϕ", title="Modes", yreversed=true)
    for m in modes
        lines!(ax2, real.(Φ[:, m]), z; label="mode $(m)")
    end
    for ax in (ax1, ax2)
        hlines!(ax, [env.depth]; color=:gray, linestyle=:dash)
    end
    linkyaxes!(ax1, ax2)
    ylims!(ax1, zmax, 0)
    hideydecorations!(ax2; grid=false)
    Legend(fig[1, 3], ax2)

    return fig
end

function Kraken.plot_ssp!(ax, env::UnderwaterEnv; npoints=500, kwargs...)
    z = collect(range(0, env.depth; length=npoints))
    return lines!(ax, soundspeed(env.c, z), z; kwargs...)
end

function ssp_axis(fig)
    return Axis(fig[1, 1]; xlabel="Sound speed (m/s)", ylabel="Depth (m)", title="Sound speed profile", yreversed=true)
end

function Kraken.plot_ssp(env::UnderwaterEnv; npoints=500)
    fig = Figure()
    Kraken.plot_ssp!(ssp_axis(fig), env; npoints=npoints)
    return fig
end

function Kraken.plot_ssp(envs::AbstractVector{<:Pair}; npoints=500)
    fig = Figure()
    ax = ssp_axis(fig)
    for (label, env) in envs
        Kraken.plot_ssp!(ax, env; npoints=npoints, label=string(label))
    end
    Legend(fig[1, 2], ax)   # outside: any corner of a profile plot can hold a curve
    return fig
end

### Transmission loss

function Kraken.plot_tl!(ax, ranges, zr, tl::AbstractMatrix; colorrange=nothing, colormap=TL_COLORMAP, kwargs...)
    size(tl) == (length(zr), length(ranges)) || throw(
        DimensionMismatch(
            "TL is $(size(tl)) but the grid is $(length(zr)) depths × $(length(ranges)) ranges; " *
            "rows are depths, columns ranges, as transmission_loss returns",
        ),
    )
    values = finite_or_nan.(tl)
    cr = colorrange === nothing ? default_tl_range(values) : colorrange
    return heatmap!(
        ax,
        collect(ranges) ./ 1000,
        collect(zr),
        permutedims(values);
        colormap=colormap,
        colorrange=cr,
        nan_color=:transparent,
        kwargs...,
    )
end

function Kraken.plot_tl(
    sol::NormalModeSolution,
    ranges,
    zs,
    zr;
    mode=:coherent,
    nmodes=length(sol.kr),
    colorrange=nothing,
    colormap=TL_COLORMAP,
)
    tl = transmission_loss(sol, ranges, zs, zr; mode=mode, nmodes=nmodes)
    fig = Figure(; size=(800, 400))
    ax = Axis(
        fig[1, 1];
        xlabel="Range (km)",
        ylabel="Depth (m)",
        yreversed=true,
        title="Transmission loss ($(mode)), $(sol.props.freq) Hz, source at $(zs) m",
    )
    hm = Kraken.plot_tl!(ax, ranges, zr, tl; colorrange=colorrange, colormap=colormap)
    scatter!(ax, [0.0], [zs]; marker=:star5, markersize=16, color=:red)
    Colorbar(fig[1, 2], hm; label="TL (dB)")
    return fig
end

function Kraken.plot_tl_range!(ax, ranges, tl; kwargs...)
    length(tl) == length(ranges) ||
        throw(DimensionMismatch("$(length(tl)) TL values for $(length(ranges)) ranges; pass a single receiver depth"))
    return lines!(ax, collect(ranges) ./ 1000, finite_or_nan.(vec(tl)); kwargs...)
end

function Kraken.plot_tl_range(sol::NormalModeSolution, ranges, zs, zr::Real; nmodes=length(sol.kr))
    fig = Figure(; size=(800, 350))
    ax = Axis(
        fig[1, 1];
        xlabel="Range (km)",
        ylabel="TL (dB)",
        yreversed=true,
        title="Transmission loss, $(sol.props.freq) Hz, source at $(zs) m, receiver at $(zr) m",
    )
    Kraken.plot_tl_range!(ax, ranges, transmission_loss(sol, ranges, zs, zr; nmodes=nmodes); label="coherent")
    Kraken.plot_tl_range!(
        ax,
        ranges,
        transmission_loss(sol, ranges, zs, zr; nmodes=nmodes, mode=:incoherent);
        label="incoherent",
        linewidth=3,
    )
    axislegend(ax; position=:lb)
    return fig
end

### Frequency and time

function Kraken.plot_dispersion(env::UnderwaterEnv, freqs; nmodes=10, kwargs...)
    fs = collect(freqs)
    # One column per mode, NaN where the mode is below cutoff, so each line starts where its mode
    # first exists rather than being joined across frequencies it is absent from.
    speeds = fill(NaN, length(fs), nmodes)
    for (i, f) in enumerate(fs)
        f > 0 || continue
        kr = real.(kraken_jl(env, f; kwargs...).kr)
        n = min(nmodes, length(kr))
        speeds[i, 1:n] .= 2π * f ./ kr[1:n]
    end

    fig = Figure(; size=(800, 450))
    ax = Axis(fig[1, 1]; xlabel="Frequency (Hz)", ylabel="Phase speed (m/s)", title="Modal phase speed")
    for m in 1:nmodes
        all(isnan, view(speeds, :, m)) && continue
        lines!(ax, fs, speeds[:, m]; label="mode $(m)")
    end
    hlines!(ax, [minimum(env.c.c)]; color=:gray, linestyle=:dash)
    env.bottom_bc isa AcousticHalfspace && hlines!(ax, [env.cb]; color=:gray, linestyle=:dash)
    Legend(fig[1, 2], ax)
    return fig
end

function Kraken.plot_arrivals(
    env::UnderwaterEnv,
    freqs,
    ranges,
    zs,
    zr::Real,
    times;
    spectrum=nothing,
    reduction_speed=maximum(env.c.c),
    kwargs...,
)
    rs = collect(Float64, ranges)
    P = broadband_field(env, freqs, rs, zs, zr; kwargs...)   # frequency × 1 depth × range
    spacing = length(rs) > 1 ? minimum(diff(sort(rs))) / 1000 : 1.0

    fig = Figure(; size=(800, 500))
    ax = Axis(
        fig[1, 1];
        xlabel="Reduced time t − r / $(round(Int, reduction_speed)) m/s (s)",
        ylabel="Range (km)",
        title="Arrivals at $(zr) m, source at $(zs) m",
    )
    for (j, r) in enumerate(rs)
        p = synthesize_pulse(P[:, 1, j], freqs, times; spectrum=spectrum, t_reduce=r / reduction_speed)
        peak = maximum(abs, p)
        trace = iszero(peak) ? p : p ./ peak
        lines!(ax, collect(times), r / 1000 .+ 0.45 * spacing .* trace; color=:black, linewidth=0.8)
    end
    return fig
end

end # module
