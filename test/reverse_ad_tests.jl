using Test
using Kraken
using ForwardDiff
using Zygote
using LinearAlgebra

# Reverse-mode AD through the wavenumber solve. `src/kraken_ad.jl` attaches a ChainRules `rrule` to
# `solve_for_kr` that applies the implicit function theorem to the Sturm determinant, so the
# bracketing solver is never traced. These tests check the rule's content directly (against
# ForwardDiff on `det_sturm`, which involves no iteration and so is an exact reference) and then
# end to end (against ForwardDiff on the whole solve).

# --- helpers: rebuild a struct with one field replaced, promoting the rest to its type -----------
# Needed because ForwardDiff has to perturb `cb`, `ρb` and `freq` *without* rebuilding the cache —
# that is exactly the partial derivative the rule computes for the half-space term.

function env_with(env; cb=env.cb, ρb=env.ρb)
    T = promote_type(typeof(cb), typeof(ρb))
    return UnderwaterEnv{typeof(env.c),typeof(env.ρ),T}(
        env.c, env.ρ, T(cb), T(ρb), T.(env.h_vec), T.(env.layer_depth), T(env.depth)
    )
end

function props_with(props, freq)
    return AcousticProblemProperties{typeof(freq),eltype(props.Δz_vec)}(freq, props.Nz_vec, props.Δz_vec, props.zn_vec)
end

function cache_with(a, e, λ)
    T = promote_type(eltype(a), eltype(e), eltype(λ))
    a, e, λ = T.(a), T.(e), T.(λ)
    return AcousticProblemCache(a, e, λ, Tridiagonal(e[2:end], a, e[2:end]))
end

# `v` with entry `k` replaced by `x`, keeping a concrete element type — a plain
# `[j == k ? x : v[j] for j in ...]` infers `Vector{Real}` and `AcousticProblemCache` has no method
# for that.
perturb(v, k, x) = [j == k ? x : oftype(x, v[j]) for j in eachindex(v)]

pekeris(θ) = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5])...)

"""
One converged wavenumber for the Pekeris waveguide described by `θ = [c0, cb, ρ0, ρb, depth]`.

Everything but `solve_for_kr` has to be traceable by the AD backend for this to work at all — the
environment builder, the mesh, and the finite-difference coefficients are all on the differentiable
seam described at the top of `src/kraken_core.jl`.
"""
function kr_at(θ; freq=100.0, mode=1, tol=nothing)
    env = pekeris(θ)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    spans = bisection(env, props, cache)
    isnothing(tol) && return solve_for_kr(spans[mode, :], env, props, cache)
    return solve_for_kr(spans[mode, :], env, props, cache; abstol=tol, reltol=tol)
end

# Componentwise relative error, floored so a gradient entry near zero is not compared with its own
# roundoff (the density derivatives are ~1e-7 while ∂kr/∂c0 is ~1e-4).
relerr(x, y) = maximum(abs.(x .- y) ./ max.(abs.(y), 1e-12))

θ0 = [1500.0, 1600.0, 1000.0, 1500.0, 100.0]

@testset "Reverse-mode AD" begin
    @testset "sturm_sensitivities" begin
        # The rule's actual content, isolated: partial derivatives of the Sturm determinant at a
        # fixed `kr`. No root-finding is involved, so ForwardDiff is an exact reference here.
        for freq in (50.0, 100.0, 250.0)
            env = pekeris(θ0)
            props = AcousticProblemProperties(env, freq)
            cache = AcousticProblemCache(env, props)
            kr = solve_for_kr(bisection(env, props, cache)[1, :], env, props, cache)
            s = Kraken.sturm_sensitivities(kr, env, props, cache)

            # The forward sweep must reproduce the solver's own determinant *exactly*, rescaling
            # included — otherwise the two are not derivatives of the same scalar.
            @test s.D == first(det_sturm(kr, env, props, cache))

            a, e, λ = cache.a_vec, cache.e_vec, cache.λ_scaling
            D(kr_, env_, props_, cache_) = first(det_sturm(kr_, env_, props_, cache_))

            @test s.dkr ≈ ForwardDiff.derivative(x -> D(x, env, props, cache), kr) rtol = 1e-10
            @test s.dcb ≈ ForwardDiff.derivative(x -> D(kr, env_with(env; cb=x), props, cache), env.cb) rtol = 1e-10
            @test s.dρb ≈ ForwardDiff.derivative(x -> D(kr, env_with(env; ρb=x), props, cache), env.ρb) rtol = 1e-10
            @test s.dfreq ≈ ForwardDiff.derivative(x -> D(kr, env, props_with(props, x), cache), props.freq) rtol =
                1e-10

            # The three mesh-length vectors, sampled: first, last (which carries the half-space
            # boundary condition and so takes a different branch), and a few in between.
            N = length(a)
            for k in unique([1, 2, N ÷ 3, N ÷ 2, N - 1, N])
                @test s.da[k] ≈ ForwardDiff.derivative(x -> D(kr, env, props, cache_with(perturb(a, k, x), e, λ)), a[k]) rtol =
                    1e-9
                @test s.de[k] ≈ ForwardDiff.derivative(x -> D(kr, env, props, cache_with(a, perturb(e, k, x), λ)), e[k]) rtol =
                    1e-9 atol = 1e-300
                @test s.dλ[k] ≈ ForwardDiff.derivative(x -> D(kr, env, props, cache_with(a, e, perturb(λ, k, x))), λ[k]) rtol =
                    1e-9
            end
        end
    end

    @testset "wavenumber gradient, Zygote vs ForwardDiff" begin
        # Every mode at 50 and 100 Hz. The 250 Hz case gets its own testset below, because there the
        # two methods genuinely disagree and ForwardDiff is the one that is wrong.
        for freq in (50.0, 100.0)
            env = pekeris(θ0)
            props = AcousticProblemProperties(env, freq)
            nmodes = size(bisection(env, props, AcousticProblemCache(env, props)), 1)
            @test nmodes > 0
            for mode in 1:nmodes
                f = θ -> kr_at(θ; freq=freq, mode=mode)
                @test relerr(Zygote.gradient(f, θ0)[1], ForwardDiff.gradient(f, θ0)) < 1e-9
            end
        end
    end

    @testset "the rule is tolerance-independent, ForwardDiff is not" begin
        # 250 Hz, mode 1. ForwardDiff differentiates the arithmetic of the *last ITP iterate*, so its
        # gradient carries the solver's truncation error; the rule differentiates the exact root and
        # does not. Central differences on the primal settle the disagreement — see the plan's
        # Architecture Decisions.
        f = θ -> kr_at(θ; freq=250.0, mode=1)
        g_zy = Zygote.gradient(f, θ0)[1]
        g_fw = ForwardDiff.gradient(f, θ0)
        g_fd = map(1:length(θ0)) do i
            h = 1e-6 * abs(θ0[i])
            step = [j == i ? h : 0.0 for j in 1:length(θ0)]
            (f(θ0 .+ step) - f(θ0 .- step)) / (2h)
        end

        @test relerr(g_zy, g_fd) < 1e-4
        @test relerr(g_zy, g_fd) < relerr(g_fw, g_fd)   # the rule is the closer of the two

        # Tightening the root tolerance by four orders of magnitude must not move the rule's answer.
        g_tight = Zygote.gradient(θ -> kr_at(θ; freq=250.0, mode=1, tol=1e-10), θ0)[1]
        @test relerr(g_zy, g_tight) < 1e-9
    end
end
