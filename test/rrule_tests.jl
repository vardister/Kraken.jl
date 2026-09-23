using Test
using Kraken
using ChainRulesTestUtils
using ChainRulesCore
using DataInterpolations
using FiniteDifferences: FiniteDifferences, central_fdm
using LinearAlgebra
using Random: AbstractRNG

### The rules on their own, through ChainRulesTestUtils --------------------------------------------
#
# `reverse_ad_tests.jl` compares gradients of *compositions*, where a wrong term in one rule can be
# cancelled by a wrong term in another or simply never reached by the cotangent that happens to flow.
# `test_rrule` asks the narrower question: with a random cotangent, does this one pullback reproduce
# finite differences of this one primal, over every argument at once? It checks arguments the solver
# never varies independently — `props.zn_vec` against `cache.a_vec`, say — and it checks the tangent
# *types*, which is how a rule acquires a `Tangent` field that no backend ever accumulates.
#
# **Why this is its own file.** Measured 2026-09-23: the `solve_for_kr` check costs 373 s the first
# time and 44 s every time after, in the same session. That 8.5x gap is not the finite differencing
# — it is Julia compiling `test_rrule`, `to_vec`, `from_vec` and `rand_tangent` for `solve_for_kr`'s
# concrete signature, which names `UnderwaterEnv`'s six type parameters, three `LinearInterpolation`s
# among them. Nothing in the assertion is reducible and the compilation is unavoidable, but it is the
# single largest block in the AD tests — so it lives here, where `KRAKEN_SKIP_TESTS` can leave it out
# while you iterate on something it cannot reach. A full run loses nothing, and CI never skips.
#
# It needs to be told three things about these arguments before it can run at all.

# 1. `rand_tangent` and `to_vec` both recurse through a struct's fields, and the interpolant inside
#    `SampledSSP1D` bottoms out in `DataInterpolations.ExtrapolationType.T` — an `EnumX` value that
#    neither knows what to do with ("Non-struct types are not supported by this fallback"). Declaring
#    the interpolant opaque is the honest answer rather than a workaround: the rules in
#    `src/kraken_ad.jl` already return `f=NoTangent()` for it, because it is fully determined by the
#    `z` and `c` fields sitting beside it and differentiating both would count the profile twice.
ChainRulesTestUtils.rand_tangent(::AbstractRNG, ::DataInterpolations.AbstractInterpolation) = NoTangent()
FiniteDifferences.to_vec(A::DataInterpolations.AbstractInterpolation) = (Float64[], _ -> A)

# 2. `Nz_vec` is a `Vector{Int}`, and the generic struct reconstruction tries to rebuild it out of
#    perturbed `Float64`s (`TypeError: in new, expected Vector{Int64}`). It is the mesh point count —
#    the one part of `AcousticProblemProperties` that is deliberately not differentiable, see the
#    seam note at the top of `src/kraken_core.jl` — so it is carried across unperturbed.
function FiniteDifferences.to_vec(props::AcousticProblemProperties)
    v, back = FiniteDifferences.to_vec((props.freq, props.Δz_vec, props.zn_vec))
    function props_from_vec(x)
        freq, Δz_vec, zn_vec = back(x)
        return AcousticProblemProperties{typeof(freq),eltype(Δz_vec)}(freq, props.Nz_vec, Δz_vec, zn_vec)
    end
    return v, props_from_vec
end

# 3. The step has to be capped. FiniteDifferences' adaptive stepper picks ~1e-3 here, and the entries
#    of `cache.a_vec` are themselves ~2e-3 — so the default probe perturbs the discretized problem by
#    tens of percent, the root leaves the bracket it was handed, and `BracketingNonlinearSolve` says
#    so in a warning. `max_range=1e-6` keeps every probe inside the bracket while leaving five
#    decades of headroom over the 1e-11 at which cancellation would start to matter.
const RRULE_FDM = central_fdm(5, 1; max_range=1e-6)

# The two environments, duplicated from `reverse_ad_tests.jl` rather than shared through a third
# file: they are one line each, and a `test/` helper module for both files to import would be a
# heavier thing to keep correct than the duplication it saves.
pekeris(θ) = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5])...)
onelayer(θ) = UnderwaterEnv(one_layer_env(; c0=θ[1], c1=θ[2], cb=θ[3], ρ0=θ[4], ρ1=θ[5], ρb=θ[6], h0=θ[7], h1=θ[8])...)

θ0 = [1500.0, 1600.0, 1000.0, 1500.0, 100.0]
θ1 = [1500.0, 1550.0, 1600.0, 1000.0, 1500.0, 2000.0, 100.0, 20.0]

# Only the discretized problem is needed here, not a converged root — `reverse_ad_tests.jl`'s `setup`
# also runs `solve_for_kr` and throws the result away, which this file has no use for.
function discretize(envf, θ; freq)
    env = envf(θ)
    props = AcousticProblemProperties(env, freq)
    return env, props, AcousticProblemCache(env, props)
end

@testset "ChainRulesTestUtils" begin
    # `layer_mesh` is the whole rule: `collect(range(z_top, z_bot, Nz))`, three arguments, one of
    # them an `Int` that must come back `NoTangent`. Nothing else in the suite tests it in isolation.
    @testset "layer_mesh" begin
        test_rrule(Kraken.layer_mesh, 0.5, 100.0, 12; check_inferred=false)
    end

    # `atol = 1e-8` and not smaller because of one specific entry: `∂kr/∂ρb` is ~1.3e-6 on
    # `one_layer_env`, three decades below the largest cotangent in the same call, and finite
    # differences resolve it to about 2e-9 absolute. Everything else clears `rtol = 1e-6` on its own.
    @testset "solve_for_kr — $name" for (name, envf, θ) in (("pekeris", pekeris, θ0), ("one_layer", onelayer, θ1))
        # 25 Hz keeps the mesh at a few dozen points: `test_rrule` finite-differences *every*
        # coordinate of every argument, so the cost is quadratic in the mesh, not linear.
        env, props, cache = discretize(envf, θ; freq=25.0)
        span = bisection(env, props, cache)[1, :]
        test_rrule(
            Kraken.solve_for_kr, span, env, props, cache; fdm=RRULE_FDM, rtol=1e-6, atol=1e-8, check_inferred=false
        )
    end

    # `mode_eigenvector` is deliberately absent, and it is worth recording why rather than leaving a
    # gap. `test_rrule` requires a primal that does not mutate its arguments, and this one does:
    # `create_finite_diff_matrix!` shifts `cache.a_vec` in place and unshifts it on the way out. That
    # is invisible to a backend, which sees the cache restored, but not to `to_vec`, whose
    # reconstruction hands the same arrays back to every probe. Measured, with the cotangent fixed so
    # both sides compute the same number: differentiating with respect to `kr` alone, finite
    # differences give 3.9617 and the rule gives 3.96168; add `cache` to the argument list and the
    # same finite difference collapses to 1e-14. Two further blockers sit behind that one — the
    # primal is the last iterate of an iteration and so is not smooth at the scale the adaptive
    # stepper probes (4.3's finding), and `get_g`'s `sqrt(kr² - q²)` makes `kr` inadmissible a short
    # way below the root, which the stepper's magnitude estimation walks straight into
    # (`DomainError`). The direct check on this rule is the hand-stepped "eigenvector adjoint
    # content" testset in `reverse_ad_tests.jl`, which covers every partial family entry by entry.
end
