using Test
using Kraken
using ForwardDiff
using Zygote
using LinearAlgebra

# Reverse-mode AD through the solver. `src/kraken_ad.jl` attaches two ChainRules `rrule`s, so that
# neither the bracketing solver nor the inverse iteration is ever traced:
#
#   * `solve_for_kr`     — implicit function theorem on the Sturm determinant
#   * `mode_eigenvector` — eigenvector adjoint via a bordered tridiagonal solve
#
# Each is checked twice: the rule's content in isolation, and then end to end against ForwardDiff on
# the whole solve. The two isolated checks use different references, for a reason worth stating.
# `det_sturm` involves no iteration, so ForwardDiff on it is exact and is the reference for the first
# rule. `mode_eigenvector` *is* an iteration, and ForwardDiff differentiates the arithmetic of its
# last iterate rather than the eigenvector it converges to — so there central differences on the
# primal are the reference, and they are what shows ForwardDiff to be the inaccurate side (see the
# plan's Architecture Decisions). End to end the distinction stops mattering: the mode-shape
# gradients below agree to 1e-11.

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

# A real sediment layer over the half-space, so the mode shapes are tested against something the
# analytic Pekeris solution does not cover: `a_vec` picks up interface averaging, the mesh is built
# in two segments, and the density is piecewise rather than constant.
onelayer(θ) = UnderwaterEnv(one_layer_env(; c0=θ[1], c1=θ[2], cb=θ[3], ρ0=θ[4], ρ1=θ[5], ρb=θ[6], h0=θ[7], h1=θ[8])...)

# Two sediment layers, each with a sound-speed *gradient* inside it. This is the case that exposes
# the known gap at the bottom of this file.
function twolayer(θ)
    return UnderwaterEnv(
        two_layer_slope_env(;
            c0=θ[1],
            c1_1=θ[2],
            c1_2=θ[3],
            c2_1=θ[4],
            c2_2=θ[5],
            cb=θ[6],
            ρ0=θ[7],
            ρ1=θ[8],
            ρ2=θ[9],
            ρb=θ[10],
            h0=θ[11],
            h1=θ[12],
            h2=θ[13],
        )...,
    )
end

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

# Central difference in a single argument, stepped relative to the point itself.
central(f, x; h=1e-6) = (f(x * (1 + h)) - f(x * (1 - h))) / (2 * x * h)

θ0 = [1500.0, 1600.0, 1000.0, 1500.0, 100.0]
θ1 = [1500.0, 1550.0, 1600.0, 1000.0, 1500.0, 2000.0, 100.0, 20.0]
θ2 = [1500.0, 1550.0, 1580.0, 1600.0, 1650.0, 1800.0, 1000.0, 1500.0, 1600.0, 2000.0, 100.0, 20.0, 20.0]

"""
The discretized problem for `envf(θ)` at `freq`, and one converged wavenumber of it.
"""
function setup(envf, θ; freq=100.0, mode=1)
    env = envf(θ)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    kr = solve_for_kr(bisection(env, props, cache)[mode, :], env, props, cache)
    return env, props, cache, kr
end

"""
`∫ψ dz` for one mode — a scalar functional of a *mode shape*, which is what the eigenvector rule has
to get right. Chosen over something like `sum(abs2, ψ)` because it does not scale with the number of
mesh points, so its value is comparable across environments and its finite differences do not jump
when a perturbation moves `Nz`.
"""
function mode_integral(envf, θ; freq=100.0, mode=1)
    env, props, cache, kr = setup(envf, θ; freq=freq, mode=mode)
    _, ψ = inverse_iteration(kr, env, props, cache)
    return Kraken.integral_trapz(ψ, reduce(vcat, props.zn_vec))
end

"""
A loss over *every* mode at once, so the gradient also has to travel the vectorized paths —
`find_kr` and the `Vector` method of `inverse_iteration`, both of which used to fill preallocated
arrays and had to be rewritten for this milestone.
"""
function all_modes_loss(envf, θ; freq=100.0)
    env = envf(θ)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    krs = find_kr(env, props, cache)
    _, modes = inverse_iteration(krs, env, props, cache)
    return sum(abs2, modes) + 1e6 * sum(krs)
end

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

    @testset "eigenvector adjoint content — $name" for (name, envf, θ) in
                                                       (("pekeris", pekeris, θ0), ("one_layer", onelayer, θ1))
        # The rule's actual content, isolated: the pullback of a cotangent on the *unnormalized*
        # eigenvector, against central differences on the primal. Every partial family gets checked
        # — the three cache vectors, and the half-space parameters that reach the matrix only
        # through `get_g`. Central differences settle at ~2e-6 here, which is what sets the
        # tolerance; the rule is the more accurate of the two.
        env, props, cache, kr = setup(envf, θ)
        a, e, λ = copy(cache.a_vec), copy(cache.e_vec), copy(cache.λ_scaling)
        N = length(a)
        w = [cos(1.7k) for k in 1:N]   # an arbitrary but deterministic cotangent

        _, pullback = Zygote.pullback(
            (kr_, env_, props_, cache_) -> mode_eigenvector(kr_, env_, props_, cache_)[2], kr, env, props, cache
        )
        _, denv, dprops, dcache = pullback(w)

        L(env_, props_, cache_) = dot(w, mode_eigenvector(kr, env_, props_, cache_)[2])

        @test denv.cb ≈ central(x -> L(env_with(env; cb=x), props, cache), env.cb) rtol = 1e-5
        @test denv.ρb ≈ central(x -> L(env_with(env; ρb=x), props, cache), env.ρb) rtol = 1e-5
        @test dprops.freq ≈ central(x -> L(env, props_with(props, x), cache), props.freq) rtol = 1e-5

        # `e_vec[1]` multiplies the Sturm sequence's `p0 = 0` and appears in no entry of the
        # matrix, so its cotangent is structurally zero rather than merely small.
        @test dcache.e_vec[1] == 0

        # Looser than the half-space checks above, and the slack belongs to the *reference*: the
        # primal's iteration truncation makes these central differences converge only first order in
        # the step, so at h = 1e-6 they still sit ~1e-5 from their limit on the worst entry. Shrinking
        # the step walks each estimate monotonically onto the rule's value.
        for k in unique([1, 2, N ÷ 3, N ÷ 2, N - 1, N])
            @test dcache.a_vec[k] ≈ central(x -> L(env, props, cache_with(perturb(a, k, x), e, λ)), a[k]) rtol = 1e-4
            @test dcache.e_vec[k] ≈ central(x -> L(env, props, cache_with(a, perturb(e, k, x), λ)), e[k]) rtol = 1e-4 atol =
                1e-300
            @test dcache.λ_scaling[k] ≈ central(x -> L(env, props, cache_with(a, e, perturb(λ, k, x))), λ[k]) rtol =
                1e-4
        end
    end

    @testset "mode shape gradient, Zygote vs ForwardDiff — $name" for (name, envf, θ) in (
        ("pekeris", pekeris, θ0), ("one_layer", onelayer, θ1)
    )
        # End to end, where the two methods do agree: the mode-shape gradient travels the eigenvector
        # rule, the wavenumber rule, the energy normalization (which is *not* in either rule — it is
        # traced), and the environment builder.
        for mode in 1:3
            f = θ_ -> mode_integral(envf, θ_; mode=mode)
            @test relerr(Zygote.gradient(f, θ)[1], ForwardDiff.gradient(f, θ)) < 1e-8
        end

        # …and once through the vectorized paths, all modes at once.
        f_all = θ_ -> all_modes_loss(envf, θ_)
        @test relerr(Zygote.gradient(f_all, θ)[1], ForwardDiff.gradient(f_all, θ)) < 1e-8
    end

    @testset "known gap: a profile that varies inside a layer whose depth is a parameter" begin
        # `two_layer_slope_env` has a sound-speed gradient *within* each sediment layer, so moving a
        # layer boundary moves the interpolant's knots as well as the mesh. `DataInterpolations`'
        # ChainRules rule treats an interpolant's knot vector as non-differentiable, so Zygote
        # silently drops that term — it returns a gradient, just not the whole one. ForwardDiff and
        # central differences both carry it.
        #
        # This is upstream of both rules in this file: it hits the wavenumber gradient too, which is
        # why it does not show up as an eigenvector-rule failure. The fix is an `rrule` for
        # `SampledSSP`/`SampledDensity`, which the plan carries as its own task (4.5).
        #
        # Which side is wrong was settled against Fortran, not against ForwardDiff — finite
        # differences of `kraken.exe`'s own `Re(kᵣ)` across perturbed `.env` files put `∂kr₁/∂h0` at
        # 1.6344e-5, where ForwardDiff says 1.6342e-5 and Zygote says 1.3836e-5; on `h1` and `h2`
        # Zygote has the sign wrong. The numbers are tabulated in the plan's Architecture Decisions.
        # The cause, isolated: the sound speed at a fixed depth inside the first sediment layer, as
        # that layer's upper boundary moves.
        c_at(x) = soundspeed(twolayer([θ2[1:10]; x; θ2[12:13]]).c, 110.0)
        @test ForwardDiff.derivative(c_at, 100.0) ≈ central(c_at, 100.0) rtol = 1e-6
        @test_broken isapprox(Zygote.gradient(c_at, 100.0)[1], central(c_at, 100.0); rtol=1e-6)

        # The consequence, at the level a user sees it: every sound-speed and density derivative is
        # right, and only the layer-thickness ones (11:13) are wrong. Compared against the gradient's
        # own scale rather than entrywise — this environment's wavenumber gradient spans seven orders
        # of magnitude, and the smallest entries are below the precision either method reaches.
        relerr_norm(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))
        for f in (θ_ -> last(setup(twolayer, θ_)), θ_ -> mode_integral(twolayer, θ_))
            g_zy = Zygote.gradient(f, θ2)[1]
            g_fw = ForwardDiff.gradient(f, θ2)
            @test relerr_norm(g_zy[1:10], g_fw[1:10]) < 1e-8
            @test_broken relerr_norm(g_zy[11:13], g_fw[11:13]) < 1e-8
        end
    end
end
