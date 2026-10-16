# Automatic differentiation

Kraken.jl is differentiable end to end, in **both** modes. You can take the derivative of a
wavenumber, of a mode shape, or of any scalar built from them, with respect to any environment
parameter: sound speeds, densities, layer thicknesses, frequency.

Which mode to use depends on how many parameters you are differentiating with respect to:

| | Cost | Use it for |
|---|---|---|
| **Forward** ([ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)) | one solve *per parameter* | a handful of parameters, such as group speeds or a sensitivity to bottom speed |
| **Reverse** ([Zygote.jl](https://github.com/FluxML/Zygote.jl), [Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl)) | a fixed multiple of one solve, *independent of parameter count* | many parameters, such as a whole sound-speed profile |

Nothing needs to be enabled. Kraken.jl depends on
[ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl), which carries the reverse-mode
rules; Zygote, ForwardDiff and Mooncake are yours to add.

## Forward mode: group speeds

Group speed is the derivative of angular frequency ``\omega = 2\pi f`` with respect to the horizontal
wavenumber ``k_{r,m}``, so it is a derivative with respect to a *single* parameter, which forward mode
computes in one pass.

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

# Build an environment from an N-point sound-speed profile on a fixed depth grid.
function profile_env(cvec, z)
    N = length(z)
    # Columns are [z, cp, cs, ρ, αp, αs], the KRAKEN .env SSP record layout.
    ssp = hcat(z, cvec, zero(cvec), fill(1000.0, N), zero(cvec), zero(cvec))
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
    is one of the few things reverse mode cannot trace: it stores its endpoints in an
    extended-precision form that has no derivative rule. `collect` it once, outside.

Forward mode agrees, which is how the rules are tested:

```@example ad
g_fwd = ForwardDiff.gradient(cc -> modal_sum(cc, z), c)
@printf("reverse vs forward: %.2e relative\n", maximum(abs.(g .- g_fwd) ./ max.(abs.(g_fwd), 1e-12)))
```

## Cost against parameter count

The timings below are from `test/performance_tests.jl` on a 2021 M1 laptop. The benchmark
differentiates the sum of the wavenumbers of a Pekeris-like waveguide at 100 Hz with respect to an
`N`-point profile: `N` inputs and a single output, `M = 1` (throughout, `N` counts inputs and `M`
outputs). It uses the single-mesh [`bisection`](@ref)/[`solve_for_kr`](@ref) path rather than
[`kraken_jl`](@ref), so the ratios isolate the cost of differentiation from the mesh-refinement loop:

| `N` | primal | forward | reverse | forward / primal | reverse / primal |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.081 ms | 0.092 ms | 0.431 ms | 1.1× | 5.3× |
| 5 | 0.063 ms | 0.164 ms | 0.410 ms | 2.6× | 6.5× |
| 10 | 0.067 ms | 0.346 ms | 0.403 ms | 5.2× | 6.0× |
| 50 | 0.066 ms | 1.689 ms | 0.381 ms | 25.7× | 5.8× |
| 100 | 0.082 ms | 3.86 ms | 0.463 ms | 47.3× | 5.7× |
| 500 | 0.079 ms | 20.0 ms | 0.424 ms | 253× | 5.4× |

(The test asserts the first four rows; the last two are the same benchmark run out further.)

Forward mode's multiple of the primal grows linearly with `N`. Reverse mode's stays between 5.3× and
6.5× across a 500-fold range of `N`. The two cross between 10 and 50 parameters. At `N = 50` reverse
mode is 4.4× faster than forward mode, and at `N = 500` it is **47×** faster.

Most of that 6× is fixed overhead. Zygote spends about 0.3 ms building and replaying the tape
regardless of problem size, and the 100 Hz Pekeris solve above takes 0.07 ms for 5 modes, so the
overhead dominates the ratio. Holding the parameter count at 50 and raising the frequency instead,
the multiple drops from 46× to 1.9×:

| Frequency | Modes | primal | reverse | reverse / primal |
|---:|---:|---:|---:|---:|
| 25 Hz | 1 | 0.008 ms | 0.349 ms | 46× |
| 50 Hz | 2 | 0.019 ms | 0.322 ms | 17× |
| 100 Hz | 5 | 0.064 ms | 0.431 ms | 6.7× |
| 200 Hz | 9 | 0.214 ms | 0.667 ms | 3.1× |
| 400 Hz | 18 | 0.835 ms | 1.60 ms | 1.9× |

A gradient costs a small multiple of a solve, and the multiple shrinks as the solve grows: at 400 Hz,
with 18 modes, the gradient of all 50 inputs costs less than two solves.

The benchmark runs as part of the optional performance suite and asserts the *shape* of this
scaling, so a regression that reintroduced per-parameter cost fails CI:

```bash
KRAKEN_RUN_PERFORMANCE_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
```

## Worked example: recovering a sound-speed profile

This example generates synthetic wavenumbers from a known profile, discards the profile, and
recovers it from the wavenumbers alone by gradient descent. The data are the lowest few modal
wavenumbers at three frequencies: 14 numbers for 8 unknowns.

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

The misfit is the squared wavenumber residual. The `1e-12` divisor only rescales: it measures the
residual in units of ``10^{-6}\,\mathrm{m}^{-1}``, so the printed misfit is not a column of
``10^{-8}``-sized numbers:

```@example ad
misfit(cvec) = sum(sum(abs2, model_kr(cvec, f, n) .- d) for (f, n, d) in zip(freqs, nmodes, data)) / 1e-12
nothing # hide
```

The starting prior has the linear trend and none of the sinusoidal structure, as a climatological
profile would. Adam is written out in full rather than taken from a package, to show that the *only*
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

Drawn with [`plot_ssp`](@ref), the recovered profile lies on top of the truth, and the prior's
missing structure is plain:

```@example ad
using CairoMakie
plot_ssp([
    "prior" => profile_env(c_prior, z_inv),
    "recovered" => profile_env(c_hat, z_inv),
    "truth" => profile_env(c_true, z_inv),
])
```

Each iteration of that loop costs one gradient, which the benchmark above puts at 5 to 7 solves
whether the profile has 8 points or 500. Forward mode's cost per iteration grows with the point
count, to 253 solves at 500 points.

!!! warning "The misfit is non-convex, so the starting point matters"
    Modal inversion has local minima, and reverse-mode AD does not remove them. Starting the same
    loop from a flat 1500 m/s profile instead of the linear prior gets stuck about 18 m/s away from
    the truth with a misfit that has stopped falling. The minimum belongs to the misfit, and exact
    gradients do not change it: from a start within a few m/s of the truth the same optimizer
    converges to ``10^{-3}`` m/s. Real inversions handle this with a prior, multiple restarts, or a
    global search on top, and each of those still runs on the same gradients.

## Worked example: a seabed sound speed from transmission loss

The inversion above fits *wavenumbers*, and no hydrophone measures a wavenumber. This one fits
transmission loss, which is what a hydrophone does measure, and solves for the quantity a survey
usually wants out of it: the sound speed of the seabed. What the field adds on top of the solver,
the mode sum and the depth interpolation, is differentiated here without a rule of its own.

The water column is taken as known and the seabed sound speed is the unknown:

```@example ad
# The same environment as above, with the seabed sound speed opened up.
function seabed_env(cb, cvec, z)
    n = length(z)
    ssp = hcat(z, cvec, zero(cvec), fill(1000.0, n), zero(cvec), zero(cvec))
    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        depth cb 0.0 1800.0 0.0 0.0
    ]
    return UnderwaterEnv(ssp, [0.0 0.0 depth], sspHS)
end
nothing # hide
```

The data are transmission loss on a 1–8 km track, at two receiver depths and two frequencies:

```@example ad
tl_ranges = collect(1_000.0:250.0:8_000.0)
tl_depths = [30.0, 70.0]
tl_freqs = (75.0, 150.0)
tl_nmodes = (2, 5)
cb_true = 1650.0

function tl_of(cb, freq, nmodes)
    sol = kraken_jl(seabed_env(cb, c_true, z_inv), freq; abstol=1e-10, reltol=1e-10)
    return transmission_loss(sol, tl_ranges, 25.0, tl_depths; nmodes, mode=:incoherent)
end

tl_data = [tl_of(cb_true, f, n) for (f, n) in zip(tl_freqs, tl_nmodes)]
tl_misfit(cb) = sum(sum(abs2, tl_of(cb, f, n) .- d) for (f, n, d) in zip(tl_freqs, tl_nmodes, tl_data))
nothing # hide
```

Two choices in that code keep the misfit smooth enough for gradient descent to work.

!!! note "Fit incoherent TL, not coherent"
    A coherent mode sum has the interference nulls in it, and a null moves in range as `cb` changes,
    so the squared residual between two coherent TL curves is multimodal: scanned at 5 m/s
    over 1560–1760 m/s it shows eight local minima, 15 to 35 m/s apart, and its value at 1600 m/s
    (1205) is *lower* than at 1605 (1347) even though 1605 is closer to the truth. Gradient descent on that
    lands wherever it started. The incoherent sum has thrown the phase away, so what is left is the
    smooth mean level, and the misfit falls monotonically to zero from both sides. Matched-field
    inversion smooths or averages before it optimizes for the same reason. The gradients of the
    coherent misfit are correct; the function they describe is the problem.

!!! note "Fix `nmodes`, for the same reason the mesh schedule is fixed"
    The *number* of trapped modes is an integer function of `cb`: over the same span the 150 Hz solve
    goes from 5 modes to 10, one at a time. Each new mode arrives with a finite amount of energy, so
    a misfit summed over whatever number of modes each solve returns is a step function: it jumps by a factor of
    four at 1600→1605 m/s and again at 1675→1680. Summing a *fixed* number of modes is the
    derivative that is well defined, and it is the same integer-valued discontinuity already
    documented for the mesh schedule under [Caveats](@ref).

With those settled, the misfit is a smooth bowl with its minimum at the truth:

```@example ad
@printf("%10s %14s\n", "cb (m/s)", "misfit")
for cb in 1560.0:25.0:1760.0
    @printf("%10.0f %14.4f\n", cb, tl_misfit(cb))
end
```

The optimizer is the same Adam, driven by the same one line of Kraken.jl:

```@example ad
function invert_cb(cb0; iters=200, α=3.0, decay=0.995)
    cb = cb0
    m = 0.0
    v = 0.0
    history = NTuple{3,Float64}[]
    for t in 1:iters
        g = Zygote.gradient(tl_misfit, cb)[1]    # <- the whole of the Kraken.jl side
        m = 0.9 * m + 0.1 * g
        v = 0.9 * v + 0.1 * g^2                  # β₂ = 0.9, not 0.999; see below
        cb -= (α * decay^t) * (m / (1 - 0.9^t)) / (sqrt(v / (1 - 0.9^t)) + 1e-12)
        iszero(t % 50) && push!(history, (t, tl_misfit(cb), cb))
    end
    return cb, history
end

cb_hat, tl_history = invert_cb(1560.0)

@printf("%12s %14s %14s\n", "iteration", "misfit", "cb (m/s)")
@printf("%12d %14.3e %14.2f\n", 0, tl_misfit(1560.0), 1560.0)
for (t, L, cb) in tl_history
    @printf("%12d %14.3e %14.2f\n", t, L, cb)
end
@printf("\ntruth %.1f m/s, recovered %.2f m/s, error %.2f m/s\n", cb_true, cb_hat, cb_hat - cb_true)
```

The fit, at the shallower receiver and 150 Hz. The starting guess is only about 0.2 dB off,
consistently along the whole track. That small, systematic offset is what the inversion has to work
with, and the recovered curve falls on the data:

```@example ad
fig = Figure(; size=(900, 350))
ax = Axis(fig[1, 1]; xlabel="Range (km)", ylabel="Incoherent TL (dB)", yreversed=true,
          title="Receiver at $(tl_depths[1]) m, $(tl_freqs[2]) Hz")
plot_tl_range!(ax, tl_ranges, tl_data[2][1, :]; label="data (cb = $(cb_true) m/s)", linewidth=5, color=(:gray, 0.6))
plot_tl_range!(ax, tl_ranges, tl_of(1560.0, tl_freqs[2], tl_nmodes[2])[1, :]; label="start (cb = 1560 m/s)", linestyle=:dash)
plot_tl_range!(ax, tl_ranges, tl_of(cb_hat, tl_freqs[2], tl_nmodes[2])[1, :];
               label=@sprintf("recovered (cb = %.1f m/s)", cb_hat), color=:black)
axislegend(ax; position=:rt)
fig
```

The 90 m/s error in the starting guess ends as a fraction of 1 m/s. Starting from the other side,
at 1750 m/s, the inversion reaches 1650.0001 m/s with the misfit at ``10^{-11}``.

!!! tip "Adam's β₂ has to be shortened here"
    The textbook ``\beta_2 = 0.999`` fails on this misfit without raising an error. The gradient
    falls by about four orders of magnitude between the starting point and the minimum, and a
    second-moment memory a thousand steps long keeps the large early gradients in `v` long after
    they are gone. `m / √v` collapses to a small fraction of 1, and the step size collapses with it.
    The search still descends monotonically, but the same 200 iterations reach 1622 m/s instead of
    1650. ``\beta_2 = 0.9`` forgets fast enough to follow the gradient's own scale. The same applies
    to any misfit whose gradient spans several decades between the start and the minimum.

!!! warning "Transmission loss does not determine the whole environment"
    Open the water-column profile up as an unknown alongside `cb` (nine parameters instead of one,
    on the same data) and the optimizer drives the misfit from 143 down to ``1.6\times10^{-8}``,
    a fit ten orders of magnitude better than anything above. It lands on `cb = 1593` m/s and a
    profile whose RMS error is *worse* than the prior it started from: 12.0 m/s against 5.3. The
    optimizer did its job. It found a different environment that reproduces the data to round-off.
    Adding a nine-depth vertical array and four frequencies narrows the gap but does not close it.
    This is non-uniqueness in the physics, which more accurate gradients cannot fix; priors,
    regularization and joint inversion of independent data address it. The gradients make each
    attempt cheap. They do not make the problem well posed.

## Where the derivative rules attach

The solver is not traced. Reverse mode enters through explicit rules attached at four seams, so the
internals (in-place cache mutation, an `argmax`, a sign-fixing branch, an iterative refinement)
never had to change:

| Seam | Rule |
|---|---|
| [`solve_for_kr`](@ref) | Implicit function theorem on the Sturm-sequence determinant: ``\partial k_r/\partial\theta = -D_\theta/D_{k_r}`` at the converged root |
| [`inverse_iteration`](@ref) | Eigenvector adjoint, a bordered tridiagonal solve that projects out the ``\psi`` direction. Linear in the mesh size, so a mode-shape gradient costs about what a wavenumber gradient costs |
| [`soundspeed`](@ref), [`density`](@ref) | Linear-interpolant adjoint over the values, the **knots**, and the query depth |
| [`bisection`](@ref) | Non-differentiable: mode counting is integer-valued and piecewise constant, so it correctly contributes zero |
| [`acoustic_field`](@ref), [`transmission_loss`](@ref), [`broadband_field`](@ref) | **No rule.** The field is arithmetic over the wavenumbers and mode shapes, so it is traced, and the cotangent reaches the rules above unchanged. Its depth interpolation is a dozen lines of linear weights rather than an interpolant object, which is what makes tracing it viable |

Because of the implicit-function rule, **the wavenumber derivative does not depend on how hard the
root solver worked.** It differentiates the equation the root satisfies, not
the arithmetic of the search. Forward mode gets the same rule from `ext/KrakenForwardDiffExt.jl`, which
loads with ForwardDiff. After the plain solve, it takes one Newton step in `Dual` arithmetic, which
sets the partials to ``-D_\theta/D_{k_r}`` and leaves the value unchanged. Before that extension
existed, ForwardDiff differentiated the search itself, and at 250 Hz it was over 1% from finite
differences. Now the two modes agree to rounding.

Mooncake reaches the same rules through `ext/KrakenMooncakeExt.jl`, a package extension that loads
only if you load Mooncake. Mooncake does not read `ChainRulesCore` rules directly, so the extension
declares each one a Mooncake primitive.

## Caveats

* **The top-level gradient is conditional on the mesh schedule.** [`kraken_jl`](@ref) refines the
  depth mesh until the wavenumbers converge, and how many levels that takes is an integer decision in
  the parameters. It is held fixed, which is the standard treatment, but it makes the gradient
  discontinuous at parameter values where the level count changes. Finite differences taken across
  such a point measure the jump rather than the slope. If AD and finite differences disagree on a
  thickness or a bottom speed, check whether the mesh moved. `dont_break=true` fixes the schedule by
  construction.
* **A field summed over every trapped mode is a step function of the parameters.** The number of
  trapped modes is an integer that changes as a parameter crosses a cutoff, and the mode arriving at
  that cutoff carries a finite amount of energy, so the field jumps the way it does when the mesh
  schedule changes.
  Pass an explicit `nmodes` when you differentiate a field, low enough that the solve returns at
  least that many everywhere the search goes. The worked example above measures the size of the jump
  it avoids.
* **Only the wavenumbers are Richardson-extrapolated**; mode shapes come from the coarsest mesh. So a
  mode-shape gradient from `kraken_jl` is the level-1 solve's.
* **Do not loosen the solver tolerance for top-level gradients.** The implicit rule is evaluated at
  the level-1 root, which is itself a solver output, so here its tolerance independence does not
  apply: a root converged to `abstol = reltol = 1e-6` moved the rule's gradient by 4.6e-6 and
  ForwardDiff's by 6.7e-6. `kraken_jl`'s default is `1e-10`, at which the two modes agree to `1e-10`.
* **A lossy environment makes `kr` complex, which narrows what you may differentiate.** The rule is
  the standard one: real parameters in, real loss out, complex only in between. `sum(imag, sol.kr)`,
  `sum(real, sol.kr)` and `sum(abs2, sol.kr)` are all fine and every backend agrees on them. A
  *complex-valued* loss is not: Zygote raises `"Output is complex, so the gradient is not defined."`
  rather than picking a conjugation convention. Take `real`, `imag` or `abs2` first. (ForwardDiff
  returns a `Complex` derivative for the same function, so the two backends disagree about whether
  the question is well posed.)
* **Do not differentiate attenuation at ``\alpha = 0``.** ``\mathrm{Im}(k_r)`` is identically zero
  below zero attenuation and linear above it, so that point is a kink with no two-sided derivative,
  and the two modes take *opposite sides* without a warning. ForwardDiff sees the perturbation
  through `is_lossy` and returns the right-hand derivative; Zygote evaluates `is_lossy` on the primal
  and returns the left-hand one, which is zero. Just above zero they agree to ``10^{-15}``.
  Differentiate at a nonzero attenuation, where the derivative is two-sided.
* **Reverse mode over a lossy environment means Zygote.** Mooncake cannot trace the complex
  arithmetic that attenuation introduces and fails with
  `"It is not permissible to bitcast to a differentiable type during AD"`. The failing call is inside
  Mooncake's own `Complex` support, so there is no rule Kraken.jl could add to fix it. Lossless
  solves are unaffected and both backends handle them.

## Example script and derivations

[`examples/reverse_ad.jl`](https://github.com/vardister/Kraken.jl/blob/master/examples/reverse_ad.jl)
in the repository is a runnable script that covers this page in more depth. It includes the
tolerance-independence demonstration, mode-shape functionals, and a case where the layer-thickness
derivatives were validated against central differences of unmodified Fortran `kraken.exe`.

[`docs/notebooks/implicit-differentiation.nb`](https://github.com/vardister/Kraken.jl/blob/master/docs/notebooks/implicit-differentiation.nb)
is a Wolfram notebook that derives each rule symbolically and checks it numerically against this
implementation. It shows where each rule comes from, why the Sturm rescaling cancels, why the
eigenvector solve has to be projected, and what dropping the interpolant term costs.
