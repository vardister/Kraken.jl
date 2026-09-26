# A broadband frequency sweep on several threads
# ==============================================
#
# WHAT THIS SHOWS
#   How to solve one environment at many frequencies in parallel. Every frequency is an independent
#   solve, so a sweep spreads over threads with nothing but `Threads.@threads`. The script checks
#   that the threaded answers equal the serial ones exactly and reports the speed-up.
#
# HOW TO RUN IT
#   From the repository root, once per clone:
#
#       julia --project=examples -e 'using Pkg; Pkg.instantiate()'
#
#   and then, with as many threads as you have cores:
#
#       julia --threads=auto --project=examples examples/broadband_sweep.jl
#
#   The figure is written to `examples/figures/`. On Julia 1.10, which ignores the `[sources]` entry in
#   `examples/Project.toml`, run `Pkg.develop(path=".")` in that environment first.

using Kraken
using CairoMakie

#%% The environment and the frequencies
env = UnderwaterEnv(pekeris_env()...)   # 100 m of 1500 m/s water over a 1600 m/s half-space
freqs = collect(5.0:1.0:400.0)

#%% Serial, then threaded
#
# `kraken_jl` allocates everything it works on, so separate calls share no mutable state and can
# run concurrently. Each thread writes only its own slot of `threaded`.
solve(f) = kraken_jl(env, f)
solve(freqs[1])   # compile before timing

t_serial = @elapsed serial = map(solve, freqs)

threaded = similar(serial)
t_threaded = @elapsed Threads.@threads for i in eachindex(freqs)
    threaded[i] = solve(freqs[i])
end

@assert all(s.kr == t.kr for (s, t) in zip(serial, threaded))
@info "Sweep of $(length(freqs)) frequencies" Threads.nthreads() t_serial t_threaded speedup = t_serial / t_threaded

#%% What the sweep shows
#
# The number of trapped modes grows about linearly with frequency. A new mode appears each time
# ωh/c crosses another cutoff, where h is the water depth. This is why the solve time grows with
# frequency too, and why high frequencies are where threads pay off most.
nmodes = [length(s.kr) for s in serial]

fig = Figure(; size=(700, 400))
ax = Axis(fig[1, 1]; xlabel="Frequency (Hz)", ylabel="Trapped modes", title="Pekeris waveguide")
stairs!(ax, freqs, nmodes; step=:post)

mkpath(joinpath(@__DIR__, "figures"))
save(joinpath(@__DIR__, "figures", "broadband_sweep.png"), fig)

# For the acoustic field over a band rather than the modes, `broadband_field(env, freqs, ranges, zs,
# zr)` runs a sweep like this one and sums the modes at each frequency. `synthesize_pulse` then turns
# the result into a time series.
