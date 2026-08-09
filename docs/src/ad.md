# Automatic differentiation

Kraken.jl is differentiable end to end, in **both** modes. You can take the derivative of a
wavenumber, of a mode shape, or of any scalar built from them, with respect to any environment
parameter — sound speeds, densities, layer thicknesses, frequency.

Which mode you want depends on how many parameters you are differentiating with respect to:

| | Cost | Use it for |
|---|---|---|
| **Forward** ([ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)) | one solve *per parameter* | a handful of parameters — group speeds, a sensitivity to bottom speed |
| **Reverse** ([Zygote.jl](https://github.com/FluxML/Zygote.jl), [Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl)) | a fixed multiple of one solve, *independent of parameter count* | many parameters — fitting a whole sound-speed profile |

Nothing needs to be enabled. Kraken.jl depends on
[ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl), which carries the reverse-mode
rules; Zygote, ForwardDiff and Mooncake are yours to add.

## Forward mode: group speeds

Group speed is the derivative of angular frequency ``\omega = 2\pi f`` with respect to the horizontal
wavenumber ``k_{r,m}``, so it is a derivative with respect to a *single* parameter — forward mode's
home ground.

```julia
using Kraken, ForwardDiff

function kr_vs_freq(freq)
    env = UnderwaterEnv(pekeris_env()...)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    return find_kr(env, props, cache)
end

group_speeds = 2pi ./ ForwardDiff.derivative(kr_vs_freq, 100.0)
```

These are checked against Fortran KRAKEN's own group-speed table on every push, to within 0.1%.

## Reverse mode: many parameters at once

Reverse mode gives the gradient with respect to *every* input in one pass. The example below takes a
25-point sound-speed profile as the unknown and asks how the sum of the wavenumbers responds to each
point of it:

```@example ad
using Kraken, Zygote, ForwardDiff, Printf

depth = 100.0

# Build an environment from an M-point sound-speed profile on a fixed depth grid.
function profile_env(cvec, z)
    M = length(z)
    # Columns are [z, cp, cs, ρ, αp, αs] — the KRAKEN .env SSP record layout.
    ssp = hcat(z, cvec, zero(cvec), fill(1000.0, M), zero(cvec), zero(cvec))
    layers = [0.0 0.0 depth]
    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        depth 1600.0 0.0 1500.0 0.0 0.0
    ]
    return UnderwaterEnv(ssp, layers, sspHS)
end

modal_sum(cvec, z; freq=100.0) = sum(kraken_jl(profile_env(cvec, z), freq).kr)

z = collect(range(0.0, depth, 25))
c = 1500.0 .+ 5 .* sin.(range(0, 3, 25))

g = Zygote.gradient(cc -> modal_sum(cc, z), c)[1]
@printf("largest sensitivity: %.3e at %.1f m\n", g[argmax(abs.(g))], z[argmax(abs.(g))])
```

!!! note "Build the depth grid outside the differentiated function"
    `z` is constructed before the gradient is taken. It does not depend on the unknowns, and `range`
    is one of the few things reverse mode cannot trace — it stores its endpoints in an
    extended-precision form that has no derivative rule. `collect` it once, outside.

Forward mode agrees, which is how the rules are tested:

```@example ad
g_fwd = ForwardDiff.gradient(cc -> modal_sum(cc, z), c)
@printf("reverse vs forward: %.2e relative\n", maximum(abs.(g .- g_fwd) ./ max.(abs.(g_fwd), 1e-12)))
```

## Cost against parameter count

This is the whole reason reverse mode exists here. Timings below are from
`test/performance_tests.jl` on a 2021 M1 laptop — differentiating the sum of the wavenumbers of a
Pekeris-like waveguide at 100 Hz with respect to an `M`-point profile. They use the single-mesh
[`bisection`](@ref)/[`solve_for_kr`](@ref) path rather than [`kraken_jl`](@ref), so the ratios
isolate the cost of differentiation from the mesh-refinement loop:

| `M` | primal | forward | reverse | forward / primal | reverse / primal |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.081 ms | 0.092 ms | 0.431 ms | 1.1× | 5.3× |
| 5 | 0.063 ms | 0.164 ms | 0.410 ms | 2.6× | 6.5× |
| 10 | 0.067 ms | 0.346 ms | 0.403 ms | 5.2× | 6.0× |
| 50 | 0.066 ms | 1.689 ms | 0.381 ms | 25.7× | 5.8× |
| 100 | 0.082 ms | 3.86 ms | 0.463 ms | 47.3× | 5.7× |
| 500 | 0.079 ms | 20.0 ms | 0.424 ms | 253× | 5.4× |

(The test asserts the first four rows; the last two are the same benchmark run out further.)

Forward mode's multiple of the primal grows linearly with `M`, exactly as the theory says; reverse
mode's sits flat at about 6×. The two cross near a dozen parameters, and past that the gap widens
without limit — at `M = 50` reverse mode is already 4.4× faster, at `M = 500` it is **47×** faster,
and nothing about that trend turns around.

That 6× is not a floor on reverse mode's cost, it is a floor on reverse mode's *overhead*: Zygote
spends about 0.3 ms building and replaying the tape whatever the problem is, and the 100 Hz Pekeris
solve above is small enough (5 modes, 0.07 ms) that the overhead dominates the ratio. Hold the
parameter count at 50 and make the physics bigger instead, and the multiple falls away:

| Frequency | Modes | primal | reverse | reverse / primal |
|---:|---:|---:|---:|---:|
| 25 Hz | 1 | 0.008 ms | 0.349 ms | 46× |
| 50 Hz | 2 | 0.019 ms | 0.322 ms | 17× |
| 100 Hz | 5 | 0.064 ms | 0.431 ms | 6.7× |
| 200 Hz | 9 | 0.214 ms | 0.667 ms | 3.1× |
| 400 Hz | 18 | 0.835 ms | 1.60 ms | 1.9× |

So the rule of thumb is the usual reverse-mode one — a gradient costs a small multiple of a solve —
and the multiple gets *better*, not worse, on the problems big enough for it to matter.

The benchmark runs as part of the optional performance suite and asserts the *shape* of this
scaling, so a regression that reintroduced per-parameter cost fails CI:

```bash
KRAKEN_RUN_PERFORMANCE_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
```

## Worked example: recovering a sound-speed profile

This is the use case the reverse-mode work exists for. We generate synthetic wavenumbers from a known
profile, throw the profile away, and recover it from the wavenumbers alone by gradient descent.

The data are the lowest few modal wavenumbers at three frequencies — 14 numbers for 8 unknowns.

```@example ad
z_inv = collect(range(0.0, depth, 8))

# The profile we are pretending not to know: a linear trend with structure on top of it.
c_true = 1500.0 .+ 30.0 .* (z_inv ./ depth) .- 8.0 .* sin.(2π .* z_inv ./ depth)

freqs = (50.0, 100.0, 200.0)
nmodes = (2, 4, 8)   # how many modes are "observed" at each frequency

model_kr(cvec, freq, nm) = kraken_jl(profile_env(cvec, z_inv), freq; abstol=1e-10, reltol=1e-10).kr[1:nm]

data = [model_kr(c_true, f, n) for (f, n) in zip(freqs, nmodes)]
nothing # hide
```

The misfit is the squared wavenumber residual. The `1e-12` divisor is cosmetic — it puts the loss in
units of "how many parts in ``10^{-6}`` of ``k_r`` are we off by", which is easier to read than
``10^{-8}``-sized numbers:

```@example ad
misfit(cvec) = sum(sum(abs2, model_kr(cvec, f, n) .- d) for (f, n, d) in zip(freqs, nmodes, data)) / 1e-12
nothing # hide
```

We start from a prior that has the linear trend right and misses the structure entirely — which is
the realistic situation, since a climatological profile gets you the trend and nothing else. Adam is
written out in full here rather than pulled from a package, because the point is that the *only*
thing Kraken.jl supplies is `Zygote.gradient(misfit, cvec)`:

```@example ad
c_prior = 1500.0 .+ 30.0 .* (z_inv ./ depth)

rmse(a, b) = sqrt(sum(abs2, a .- b) / length(a))

function invert(c0; iters=400, α=0.5, decay=0.995)
    cvec = copy(c0)
    m = zero(cvec)
    v = zero(cvec)
    history = NTuple{3,Float64}[]
    for t in 1:iters
        g = Zygote.gradient(misfit, cvec)[1]     # <- one solve's worth of work, all 8 derivatives
        m = 0.9 .* m .+ 0.1 .* g
        v = 0.999 .* v .+ 0.001 .* g .^ 2
        cvec = cvec .- (α * decay^t) .* (m ./ (1 - 0.9^t)) ./ (sqrt.(v ./ (1 - 0.999^t)) .+ 1e-12)
        iszero(t % 100) && push!(history, (t, misfit(cvec), rmse(cvec, c_true)))
    end
    return cvec, history
end

c_hat, history = invert(c_prior)

@printf("%12s %14s %14s\n", "iteration", "misfit", "RMSE (m/s)")
@printf("%12d %14.3e %14.3f\n", 0, misfit(c_prior), rmse(c_prior, c_true))
for (t, L, e) in history
    @printf("%12d %14.3e %14.3f\n", t, L, e)
end
```

The misfit falls by five orders of magnitude and the profile error by a factor of ~25, from 5.3 m/s
to about 0.2 m/s. Point by point:

```@example ad
@printf("%10s %12s %12s %12s\n", "depth (m)", "prior", "recovered", "truth")
for i in eachindex(z_inv)
    @printf("%10.1f %12.2f %12.2f %12.2f\n", z_inv[i], c_prior[i], c_hat[i], c_true[i])
end
```

Every iteration of that loop cost one gradient — and one gradient costs about what one forward solve
costs, whether the profile has 8 points or 500. With forward mode the same loop would have cost 8
solves per iteration, and 500 solves per iteration for a finely sampled profile. That is the
difference between an inversion you run and one you describe.

!!! warning "Starting point matters — this is a non-convex problem"
    Modal inversion has local minima, and reverse-mode AD does not remove them. Starting the same
    loop from a flat 1500 m/s profile instead of the linear prior gets stuck about 18 m/s away from
    the truth with a misfit that has stopped falling. That is a property of the physics, not of the
    gradients: from a start within a few m/s of the truth the same optimizer converges to
    ``10^{-3}`` m/s. Real inversions handle this with a decent prior, multiple restarts, or a global
    search on top — all of which are still cheaper with reverse-mode gradients driving them.

## What actually carries the derivatives

The solver is not traced. Reverse mode enters through explicit rules attached at four seams, which is
why the internals — in-place cache mutation, an `argmax`, a sign-fixing branch, an iterative
refinement — never had to change:

| Seam | Rule |
|---|---|
| [`solve_for_kr`](@ref) | Implicit function theorem on the Sturm-sequence determinant: ``\partial k_r/\partial\theta = -D_\theta/D_{k_r}`` at the converged root |
| [`inverse_iteration`](@ref) | Eigenvector adjoint, a bordered tridiagonal solve that projects out the ``\psi`` direction. `O(N)`, so a mode-shape gradient costs about what a wavenumber gradient costs |
| [`soundspeed`](@ref), [`density`](@ref) | Linear-interpolant adjoint over the values, the **knots**, and the query depth |
| [`bisection`](@ref) | Non-differentiable: mode counting is integer-valued and piecewise constant, so it correctly contributes zero |

One consequence of the implicit-function rule is worth knowing because it is a genuine improvement
over forward mode, not just parity: **the reverse-mode wavenumber derivative does not depend on how
hard the root solver worked.** It differentiates the equation the root satisfies, not the arithmetic
of the search. Forward mode differentiates the search and inherits its truncation error; at 250 Hz
the two disagree by over 1% at default tolerances, and finite differences side with reverse mode.

Mooncake reaches the same rules through `ext/KrakenMooncakeExt.jl`, a package extension that loads
only if you load Mooncake — it does not read `ChainRulesCore` rules directly, so the extension
declares each one a Mooncake primitive.

## Caveats

* **The top-level gradient is conditional on the mesh schedule.** [`kraken_jl`](@ref) refines the
  depth mesh until the wavenumbers converge, and how many levels that takes is an integer decision in
  the parameters. It is held fixed, which is the standard treatment, but it makes the gradient
  discontinuous at parameter values where the level count changes. Finite differences taken across
  such a point measure the jump rather than the slope — if you see AD and finite differences disagree
  on a thickness or a bottom speed, check whether the mesh moved. `dont_break=true` fixes the
  schedule by construction.
* **Only the wavenumbers are Richardson-extrapolated**; mode shapes come from the coarsest mesh. So a
  mode-shape gradient from `kraken_jl` is the level-1 solve's.
* **Tighten the solver tolerance for top-level gradients.** `kraken_jl`'s default
  `abstol = reltol = 1e-6` leaves ~1e-6 of truncation in the level-1 roots, and both AD modes then
  differentiate a slightly wrong root. This is the one place the implicit rule's tolerance
  independence does *not* save you, because the root it is evaluated at is itself a solver output. At
  `1e-10` the two modes agree to `1e-10`.
* **Attenuation is not modelled yet**, so no derivative with respect to ``\alpha_p`` exists. See the
  README's missing-features list.

## Going further

[`examples/reverse_ad.jl`](https://github.com/vardister/Kraken.jl/blob/master/examples/reverse_ad.jl)
in the repository is a runnable script covering all of the above in more depth, including the
tolerance-independence demonstration, mode-shape functionals, and a case where the layer-thickness
derivatives were validated against central differences of unmodified Fortran `kraken.exe`.
