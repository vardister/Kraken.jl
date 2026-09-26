# Normal modes of a deep-water waveguide
# =======================================
#
# WHAT THIS SHOWS
#   How to describe a layered environment from tabulated data — a measured sound-speed profile over
#   a sediment layer on a rigid basement — solve it, and look at the modes it traps.
#
# HOW TO RUN IT
#   From the repository root, once per clone:
#
#       julia --project=examples -e 'using Pkg; Pkg.instantiate()'
#
#   and then:
#
#       julia --project=examples examples/deep_water_modes.jl
#
#   Figures are written to `examples/figures/`. On Julia 1.10, which ignores the `[sources]` entry in
#   `examples/Project.toml`, run `Pkg.develop(path=".")` in that environment first.
#
#   The file is written in `#%%` cells, so it also runs block by block in VS Code.

using Kraken
using CairoMakie

#%% The environment
#
# The matrices follow the layout of an Acoustics Toolbox `.env` file. Each row of `ssp` is
# `[depth, cp, cs, ρ, αp, αs]`: depth in m, speeds in m/s, density in kg/m³ (not the g/cm³ of a
# `.env` file), and attenuation in the environment's units, dB/wavelength by default.
#
# Medium 1 is a 5000 m water column with a sound channel at 600 m. The profile is a deep-water CTD
# cast thinned to 21 points; Kraken.jl interpolates it linearly in c, like KRAKEN's `'C'` option.
water = [
    0 1509.52
    50 1507.24
    100 1499.39
    150 1492.42
    200 1488.11
    300 1484.21
    400 1481.89
    500 1480.60
    600 1480.20
    700 1480.27
    800 1480.57
    1000 1481.53
    1250 1483.27
    1500 1485.47
    2000 1490.99
    2500 1498.23
    3000 1506.18
    3500 1514.53
    4000 1523.21
    4500 1531.45
    5000 1536.15
]
nw = size(water, 1)
ssp_water = hcat(water, zeros(nw), fill(1030.0, nw), zeros(nw), zeros(nw))

# Medium 2 is 300 m of sediment whose sound speed rises with depth. It starts slower than the
# water just above it, which is common for soft, porous sediment.
ssp_sediment = [
    5000.0 1513.108 0.0 1751.0 0.0 0.0
    5300.0 1613.108 0.0 1751.0 0.0 0.0
]

ssp = vcat(ssp_water, ssp_sediment)

# One row per medium: `[mesh points, roughness, depth of its lower interface]`. The first two
# columns are accepted for compatibility with `.env` files; the mesh is chosen by the solver.
layers = [
    0.0 0.0 5000.0
    0.0 0.0 5300.0
]

# The half-spaces above and below. With a rigid bottom boundary the lower row is never read, but
# the matrix still needs one.
sspHS = [
    0.0 343.0 0.0 1.21 0.0 0.0
    5300.0 1749.0 0.0 1941.0 0.0 0.0
]

env = UnderwaterEnv(ssp, layers, sspHS; bottom_bc=RigidBoundary())

#%% Solve
#
# `kraken_jl` refines the depth mesh and Richardson-extrapolates the wavenumbers until they
# converge. The first mode is the one trapped closest to the sound channel axis.
freq = 15.0
sol = kraken_jl(env, freq)

nmodes = length(sol.kr)
phase_speeds = 2π * freq ./ sol.kr
@info "Solved at $freq Hz" nmodes slowest_phase_speed = minimum(phase_speeds)

# No mode travels slower than the slowest sound speed in the problem, at the sound channel axis.
# There is no upper limit here: a rigid basement reflects everything, so even modes that hit it
# almost head-on are trapped, and their phase speeds run far above any sound speed in the column.
# Over an acoustic half-space they would leak into the seabed, and the solver would stop at its
# sound speed instead.
@assert all(phase_speeds .> 1480.2)
@info "Steepest trapped mode" fastest_phase_speed = maximum(phase_speeds)

#%% Plot the profile and the first five modes
#
# `plot_modes` comes from Kraken's Makie extension, which `using CairoMakie` loaded. It draws the
# sound-speed profile beside the modes, so each mode's turning points can be read off against it:
# mode 1 stays near the channel axis, and each higher mode reaches further up and down.
mkpath(joinpath(@__DIR__, "figures"))
save(joinpath(@__DIR__, "figures", "deep_water_modes.png"), plot_modes(sol; modes=1:5))

# The same pieces are available as plain arrays for your own plots: `sol.modes` holds one mode per
# column, sampled at the mesh depths.
z = vcat(sol.props.zn_vec...)
@assert size(sol.modes) == (length(z), nmodes)
