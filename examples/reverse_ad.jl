# Reverse-mode automatic differentiation in Kraken.jl
# =====================================================
#
# WHAT THIS SHOWS
#   How to get the gradient of a horizontal wavenumber with respect to many environment parameters
#   at once, using Zygote (reverse mode), and how to convince yourself the answer is right.
#
# HOW TO RUN IT
#   From the repository root:
#
#       julia --project=test examples/reverse_ad.jl
#
#   `--project=test` is used because Zygote is a *test* dependency of this repo, not a dependency of
#   the package — Kraken.jl only depends on ChainRulesCore, which is what carries the rules. In your
#   own project you would just `Pkg.add("Zygote")` alongside Kraken.
#
#   The file is written in `#%%` cells, so it also runs block by block in VS Code or as a notebook.
#
# WHAT WORKS TODAY, AND WHAT DOES NOT YET
#   Reverse mode currently reaches `solve_for_kr` — a single converged wavenumber — and everything
#   feeding into it. It does NOT yet go through `kraken_jl`, the top-level entry point, nor through
#   `inverse_iteration`, which produces the mode shapes. Those are the next two tasks in the plan
#   (4.3 for mode shapes, 4.4 for `kraken_jl`). So every example below goes through the lower-level
#   path, which is only four lines:
#
#       props = AcousticProblemProperties(env, freq)   # pick the depth mesh
#       cache = AcousticProblemCache(env, props)       # build the finite-difference coefficients
#       spans = bisection(env, props, cache)           # bracket each mode
#       kr    = solve_for_kr(spans[m, :], env, props, cache)   # converge mode m
#
#   Forward mode (ForwardDiff) is unaffected by all of this and still works through `kraken_jl`.

using Kraken
using Zygote        # reverse mode: gradient w.r.t. all inputs in one pass
using ForwardDiff   # forward mode: one pass per input — used here as a cross-check
using Printf

# Compare two gradients componentwise. The floor keeps a near-zero entry from being judged against
# its own rounding error: in these environments ∂kr/∂c0 is ~1e-4 while ∂kr/∂ρb is ~1e-7.
relerr(a, b) = maximum(abs.(a .- b) ./ max.(abs.(b), 1e-12))

# Central finite differences: nudge each input up and down and measure the change. Slow and only
# ~7 digits accurate, but it relies on no clever reasoning at all, so it is the honest referee when
# the two AD modes disagree.
function fd_gradient(f, θ; rel=1e-6)
    return map(eachindex(θ)) do i
        h = rel * max(abs(θ[i]), 1.0)
        step = [j == i ? h : zero(eltype(θ)) for j in eachindex(θ)]
        return (f(θ .+ step) - f(θ .- step)) / (2h)
    end
end

#%% ---------------------------------------------------------------------------------------------
# 1. Pekeris waveguide: one wavenumber, five parameters
# -----------------------------------------------------------------------------------------------
# The canonical two-medium case — isovelocity water over a fluid half-space.

"""Mode `mode` at `freq` Hz for a Pekeris waveguide described by `θ = [c0, cb, ρ0, ρb, depth]`."""
function pekeris_kr(θ; freq=100.0, mode=1)
    env = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5])...)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    spans = bisection(env, props, cache)
    return solve_for_kr(spans[mode, :], env, props, cache)
end

θ_pek = [1500.0, 1600.0, 1000.0, 1500.0, 100.0]
names_pek = ["c0 (water speed)", "cb (bottom speed)", "ρ0 (water density)", "ρb (bottom density)", "depth"]

println("\n1. PEKERIS WAVEGUIDE — ∂kr₁/∂θ at 100 Hz")
println("   kr₁ = ", pekeris_kr(θ_pek))

g_rev = Zygote.gradient(pekeris_kr, θ_pek)[1]        # <- the new functionality
g_fwd = ForwardDiff.gradient(pekeris_kr, θ_pek)
g_fd = fd_gradient(pekeris_kr, θ_pek)

# How many points the depth mesh has for a given θ. Printed alongside the finite differences because
# it explains two of the rows — see the note below the table.
function mesh_points(θ; freq=100.0)
    env = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5])...)
    return AcousticProblemProperties(env, freq).Nz_vec[1]
end

@printf("\n   %-22s %14s %14s %14s %10s\n", "parameter", "reverse", "forward", "finite diff", "mesh ±h")
for i in eachindex(θ_pek)
    h = 1e-6 * max(abs(θ_pek[i]), 1.0)
    step = [j == i ? h : 0.0 for j in eachindex(θ_pek)]
    n_lo, n_hi = mesh_points(θ_pek .- step), mesh_points(θ_pek .+ step)
    @printf("   %-22s %14.6e %14.6e %14.6e %5d,%-4d\n", names_pek[i], g_rev[i], g_fwd[i], g_fd[i], n_lo, n_hi)
end
@printf("\n   reverse vs forward: %.2e relative\n", relerr(g_rev, g_fwd))

# Read the signs physically: a faster water column lowers kr, deeper water raises it, and the
# sensitivity to the bottom is small because mode 1 barely penetrates it.
#
# WHY THE FINITE-DIFFERENCE COLUMN IS OFF FOR `cb` AND `depth`, AND WHY THAT IS NOT A BUG
#   Look at the last column. For those two parameters the mesh gains a point (125 -> 126) somewhere
#   inside the finite-difference step, because the number of mesh points is `ceil(...)` of a
#   quantity that depends on both. The discretized wavenumber therefore has a small *jump* there,
#   and a difference quotient taken across the jump measures the jump rather than the slope.
#
#   Both AD modes report the derivative at a FIXED mesh, which is the meaningful quantity and the
#   standard treatment — but it does mean the gradient is discontinuous at the parameter values
#   where the mesh count changes. For `c0`, `ρ0` and `ρb` the mesh does not move and all three
#   columns agree to six digits.

#%% ---------------------------------------------------------------------------------------------
# 2. One sediment layer: three media, eight parameters
# -----------------------------------------------------------------------------------------------
# Water over an isovelocity sediment layer over a fluid half-space. This is the case that exercises
# the interface conditions — the finite-difference coefficients get averaged across the boundary —
# so it is the one worth checking beyond the Pekeris case.

function one_layer_kr(θ; freq=100.0, mode=1)
    c0, c1, cb, ρ0, ρ1, ρb, h0, h1 = θ
    env = UnderwaterEnv(one_layer_env(; c0=c0, c1=c1, cb=cb, ρ0=ρ0, ρ1=ρ1, ρb=ρb, h0=h0, h1=h1)...)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    spans = bisection(env, props, cache)
    return solve_for_kr(spans[mode, :], env, props, cache)
end

θ_1L = [1500.0, 1550.0, 1600.0, 1000.0, 1500.0, 2000.0, 100.0, 20.0]
names_1L = [
    "c0 water",
    "c1 sediment",
    "cb half-space",
    "ρ0 water",
    "ρ1 sediment",
    "ρb half-space",
    "h0 water thickness",
    "h1 sediment thickness",
]

println("\n2. ONE SEDIMENT LAYER — ∂kr₁/∂θ at 100 Hz")
println("   kr₁ = ", one_layer_kr(θ_1L))

g_rev_1L = Zygote.gradient(one_layer_kr, θ_1L)[1]
g_fwd_1L = ForwardDiff.gradient(one_layer_kr, θ_1L)

@printf("\n   %-24s %14s %14s\n", "parameter", "reverse", "forward")
for i in eachindex(θ_1L)
    @printf("   %-24s %14.6e %14.6e\n", names_1L[i], g_rev_1L[i], g_fwd_1L[i])
end
@printf("\n   reverse vs forward: %.2e relative\n", relerr(g_rev_1L, g_fwd_1L))

#%% ---------------------------------------------------------------------------------------------
# 3. A scalar loss over every mode
# -----------------------------------------------------------------------------------------------
# Reverse mode wants a single number out. Any summary of the whole mode set works — this is the
# shape a real misfit function takes.

function modal_sum(θ; freq=100.0)
    env = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5])...)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    spans = bisection(env, props, cache)
    return sum(solve_for_kr(spans[m, :], env, props, cache) for m in axes(spans, 1))
end

println("\n3. SUM OF ALL FIVE WAVENUMBERS — one gradient, all parameters")
g_sum = Zygote.gradient(modal_sum, θ_pek)[1]
@printf("   Σkr = %.10f\n", modal_sum(θ_pek))
@printf("   agreement with forward mode: %.2e relative\n", relerr(g_sum, ForwardDiff.gradient(modal_sum, θ_pek)))

#%% ---------------------------------------------------------------------------------------------
# 4. A whole sound-speed profile as the unknown
# -----------------------------------------------------------------------------------------------
# This is the case reverse mode exists for: the parameter is not a handful of scalars but an
# `M`-point profile, and you want the sensitivity of the wavenumbers to every point of it.
#
# Note the depth grid `z` is built OUTSIDE the differentiated function. It does not depend on the
# unknowns, and `range` is one of the few things reverse mode cannot trace (it stores its endpoints
# in an extended-precision form that has no derivative rule).

function env_from_profile(cvec, z; cb=1600.0, ρ0=1000.0, ρb=1500.0)
    M = length(cvec)
    depth = z[end]
    # Columns are [z, cp, cs, ρ, αp, αs] — the KRAKEN .env SSP record layout.
    ssp = hcat(z, cvec, zero(cvec), fill(ρ0, M), zero(cvec), zero(cvec))
    layers = [0.0 0.0 depth]
    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        depth cb 0.0 ρb 0.0 0.0
    ]
    return UnderwaterEnv(ssp, layers, sspHS)
end

function profile_sum(cvec, z; freq=100.0)
    env = env_from_profile(cvec, z)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    spans = bisection(env, props, cache)
    return sum(solve_for_kr(spans[m, :], env, props, cache) for m in axes(spans, 1))
end

z_prof = collect(range(0.0, 100.0, 25))
c_prof = 1500.0 .+ 5 .* sin.(range(0, 3, 25))

println("\n4. 25-POINT SOUND-SPEED PROFILE — sensitivity of Σkr to every point")
g_prof = Zygote.gradient(c -> profile_sum(c, z_prof), c_prof)[1]
@printf(
    "   agreement with forward mode: %.2e relative\n",
    relerr(g_prof, ForwardDiff.gradient(c -> profile_sum(c, z_prof), c_prof))
)
println("   first five entries: ", round.(g_prof[1:5]; sigdigits=4))
println("   (largest near the surface, where the modes have most of their energy)")

#%% ---------------------------------------------------------------------------------------------
# 5. Cost against parameter count — the honest picture
# -----------------------------------------------------------------------------------------------
# The theory: forward mode costs one solve per parameter, so it grows linearly; reverse mode costs
# roughly one solve regardless. That is the whole reason for this milestone.
#
# The practice, TODAY: reverse mode is CORRECT but not yet FAST. It is still slower than forward
# mode at these sizes. The reason is not the differentiation rules — those are cheap, about 12% on
# top of a wavenumber solve — but Zygote's cost in tracing the construction of the interpolation
# object inside `UnderwaterEnv`. Making that cheap is task 4.7 in the plan; until then, treat this
# table as a measurement of the starting point, not of what reverse mode can do.

bench(f, n=5) = (f(); minimum(@elapsed(f()) for _ in 1:n))

println("\n5. COST vs NUMBER OF PARAMETERS (milliseconds)")
@printf("\n   %5s %10s %10s %10s %10s %10s\n", "M", "primal", "forward", "reverse", "fwd/prim", "rev/prim")
for M in (5, 10, 25, 50)
    local zM = collect(range(0.0, 100.0, M))
    local cM = 1500.0 .+ 5 .* sin.(range(0, 3, M))
    f = cc -> profile_sum(cc, zM)
    t0 = bench(() -> f(cM))
    tf = bench(() -> ForwardDiff.gradient(f, cM))
    tz = bench(() -> Zygote.gradient(f, cM))
    @printf("   %5d %10.3f %10.2f %10.2f %10.1f %10.1f\n", M, t0 * 1e3, tf * 1e3, tz * 1e3, tf / t0, tz / t0)
end
println("""
   Forward mode's column grows linearly with M, exactly as predicted. Reverse mode's does not yet
   flatten, because the fixed overhead above dominates at these sizes.""")

#%% ---------------------------------------------------------------------------------------------
# 6. Where the reverse-mode time actually goes
# -----------------------------------------------------------------------------------------------
# Worth running once so the number above is not mysterious. Almost all of it is spent building the
# environment — specifically the two `DataInterpolations` objects inside `UnderwaterEnv` — and
# almost none in the solver or in the new rules.

let
    M = 50
    z50 = collect(range(0.0, 100.0, M))
    c50 = 1500.0 .+ 5 .* sin.(range(0, 3, M))

    full = bench(() -> Zygote.gradient(cc -> profile_sum(cc, z50), c50))
    envonly = bench(() -> Zygote.gradient(cc -> sum(env_from_profile(cc, z50).c.c), c50))
    interp = bench(() -> Zygote.gradient(cc -> sum(SampledSSP(z50, cc).c), c50))

    println("\n6. WHERE THE TIME GOES (M = 50)")
    @printf("   whole gradient                       %8.2f ms\n", full * 1e3)
    @printf("   ... just building UnderwaterEnv      %8.2f ms\n", envonly * 1e3)
    @printf("   ... just one SampledSSP interpolant  %8.2f ms\n", interp * 1e3)
    @printf("   ... everything else (solve + rules)  %8.2f ms\n", (full - envonly) * 1e3)
end

#%% ---------------------------------------------------------------------------------------------
# 7. A gradient the reference cannot produce: dependence on the root solver's tolerance
# -----------------------------------------------------------------------------------------------
# A property worth seeing, because it is where reverse mode is not merely as good as forward mode
# but better. The rule uses the implicit function theorem: it differentiates the *equation* the
# wavenumber satisfies, evaluated at the converged root, instead of differentiating the arithmetic
# of the root search. So its answer does not depend on how hard the solver worked.
#
# ForwardDiff, which does differentiate the search, inherits the search's truncation error. At
# 250 Hz, mode 1, the two disagree by over 1%, and finite differences confirm the rule is right.

function kr_with_tol(θ; freq=250.0, mode=1, tol=nothing)
    env = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5])...)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    spans = bisection(env, props, cache)
    isnothing(tol) && return solve_for_kr(spans[mode, :], env, props, cache)
    return solve_for_kr(spans[mode, :], env, props, cache; abstol=tol, reltol=tol)
end

println("\n7. TOLERANCE INDEPENDENCE (250 Hz, mode 1, ∂kr/∂cb)")
@printf("   %-12s %18s %18s\n", "solver tol", "reverse", "forward")
for tol in (nothing, 1e-10, 1e-14)
    f = θ -> kr_with_tol(θ; tol=tol)
    @printf("   %-12s %18.10e %18.10e\n", repr(tol), Zygote.gradient(f, θ_pek)[1][2], ForwardDiff.gradient(f, θ_pek)[2])
end
let f = θ -> kr_with_tol(θ)
    @printf("   %-12s %18.10e\n", "finite diff", fd_gradient(f, θ_pek)[2])
end
println("""
   The reverse column is constant to ~13 digits; the forward column drifts and then collapses at
   the tightest setting. Finite differences agree with the reverse column.""")

println("\nDone.")
