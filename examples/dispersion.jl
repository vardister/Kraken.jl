# Dispersion on a lossy continental shelf
# ========================================
#
# WHAT THIS SHOWS
#   How each mode's phase speed, group speed and attenuation change with frequency in a shallow,
#   lossy waveguide: phase speeds from the wavenumbers, group speeds by differentiating them with
#   ForwardDiff, and attenuation from their imaginary parts.
#
# HOW TO RUN IT
#   From the repository root, once per clone:
#
#       julia --project=examples -e 'using Pkg; Pkg.instantiate()'
#
#   and then:
#
#       julia --project=examples examples/dispersion.jl
#
#   The figure is written to `examples/figures/`. On Julia 1.10, which ignores the `[sources]` entry in
#   `examples/Project.toml`, run `Pkg.develop(path=".")` in that environment first.

using Kraken
using ForwardDiff
using CairoMakie

#%% The environment
#
# 74.4 m of isovelocity water over an 11 m sediment layer and a fast, lossy basement. The sediment
# starts *slower* than the water, which is typical of soft mud, and speeds up with depth.
# Attenuation is in dB per wavelength, the default: 0.05 in the sediment, 0.25 in the basement.
D, h = 74.4, 11.0
ssp = [
    0.0 1468.3 0.0 1030.0 0.0 0.0
    D 1468.3 0.0 1030.0 0.0 0.0
    D 1445.0 0.0 1600.0 0.05 0.0
    D+h 1555.0 0.0 1600.0 0.05 0.0
]
layers = [
    0.0 0.0 D
    0.0 0.0 D+h
]
sspHS = [
    0.0 343.0 0.0 1.21 0.0 0.0
    D+h 1830.0 0.0 2000.0 0.25 0.0
]
shelf = UnderwaterEnv(ssp, layers, sspHS)

#%% Sweep the frequency
#
# A lossy environment returns complex wavenumbers. Their real part sets the phase speed ω/Re(kᵣ).
# Their imaginary part is the modal attenuation in nepers per metre, negative because the mode
# decays with range.
#
# Group speed is dω/dRe(kᵣ). ForwardDiff differentiates straight through the solver, so the
# derivative is exact to rounding, with no finite-difference step to tune. Each mode appears only
# above its cutoff frequency, so the arrays are padded with NaN, which Makie leaves as gaps.
freqs = 5.0:1.0:100.0
nmodes = 6   # all the modes this shelf has below 100 Hz
phase = fill(NaN, length(freqs), nmodes)
group = fill(NaN, length(freqs), nmodes)
atten = fill(NaN, length(freqs), nmodes)   # dB/km

for (i, f) in enumerate(freqs)
    kr = kraken_jl(shelf, f).kr
    dkr_df = ForwardDiff.derivative(f -> real.(kraken_jl(shelf, f).kr), f)
    n = min(nmodes, length(kr))
    phase[i, 1:n] .= 2π * f ./ real.(kr[1:n])
    group[i, 1:n] .= 2π ./ dkr_df[1:n]
    atten[i, 1:n] .= -imag.(kr[1:n]) .* 1000 .* 20 / log(10)
end

# In a waveguide like this one, energy travels slower than the phase fronts do.
@assert all(isnan(g) || g < p for (g, p) in zip(group, phase))

# Each mode's group speed dips to a minimum a few Hz above its cutoff and then recovers towards the
# water sound speed. That minimum, the Airy phase, is the slowest energy in the waveguide, and a
# broadband pulse ends with it.
airy = [freqs[argmin(replace(group[:, m], NaN => Inf))] for m in 1:nmodes]
@info "Airy-phase frequency of each mode (Hz)" airy

#%% Plot
#
# `plot_dispersion(shelf, freqs)` draws the phase-speed panel on its own. Group speed and
# attenuation need the derivative and the imaginary part, so all three are drawn by hand here.
fig = Figure(; size=(1100, 420))
panels = [(phase, "Phase speed (m/s)"), (group, "Group speed (m/s)"), (atten, "Attenuation (dB/km)")]
for (j, (data, label)) in enumerate(panels)
    ax = Axis(fig[1, j]; xlabel="Frequency (Hz)", ylabel=label)
    for m in 1:nmodes
        lines!(ax, freqs, data[:, m]; label="mode $m")
    end
    j == 3 && (ax.yscale = log10)
end
Legend(fig[1, 4], content(fig[1, 1]))

mkpath(joinpath(@__DIR__, "figures"))
save(joinpath(@__DIR__, "figures", "dispersion.png"), fig)
