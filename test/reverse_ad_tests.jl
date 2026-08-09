using Test
using Kraken
using ForwardDiff
using Zygote
using Mooncake
using FiniteDiff
using DifferentiationInterface
using ChainRulesTestUtils
using ChainRulesCore
using DataInterpolations
using FiniteDifferences: FiniteDifferences, central_fdm
using LinearAlgebra
using Random: AbstractRNG

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

# The same disagreement measured against the gradient's own scale instead of entrywise. The right
# choice wherever a gradient spans many orders of magnitude — `two_layer_slope_env`'s spans seven,
# and Munk's `∂/∂cb` is 3e-6 against a 1.3 largest entry — because there the smallest entries are
# below the precision *either* method reaches and comparing them entrywise measures nothing but
# their roundoff.
relerr_norm(x, y) = maximum(abs.(x .- y)) / maximum(abs.(y))

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
The tolerances the top-level tests below run `kraken_jl` at, and the reason they are not the defaults.

`kraken_jl`'s default `abstol = reltol = 1e-6` is *looser* than what `find_kr` gets on the refinement
meshes, which are left to NonlinearSolve's own (much tighter) defaults — so at the default the level-1
roots carry ~1e-6 of truncation and everything downstream inherits it. Unlike the isolated
`solve_for_kr` case above, that moves the *rule's* answer too, by 4.6e-6: the implicit-function
derivative is tolerance-independent at a fixed root, but here the root itself has moved. ForwardDiff
shifts by 6.7e-6 over the same change. Tightening to 1e-10 removes both effects and lets the
comparison measure the rules rather than the root solver.
"""
const TOL = (abstol=1e-10, reltol=1e-10)

"""
The sum of every wavenumber `kraken_jl` returns — the whole solve, Richardson extrapolation included.
"""
kr_sum(envf, θ; freq=100.0, kws...) = sum(kraken_jl(envf(θ), freq; TOL..., kws...).kr)

"""
A loss over the whole `NormalModeSolution`, so the gradient travels the extrapolated wavenumbers *and*
the mode shapes. The `1e6` puts the two terms on a comparable scale — `kr ~ 0.4` against mode-shape
energies of order 1.
"""
function solution_loss(envf, θ; freq=100.0, kws...)
    sol = kraken_jl(envf(θ), freq; TOL..., kws...)
    return 1e6 * sum(sol.kr) + sum(abs2, sol.modes)
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

    @testset "kraken_jl end to end — $name" for (name, envf, θ) in
                                                (("pekeris", pekeris, θ0), ("one_layer", onelayer, θ1))
        # The top-level entry point, which is everything the two testsets above cover *plus* the
        # mesh-refinement loop and the Richardson extrapolation. The loop was rewritten for this: it
        # grew `h_meshes` and `krs_all` by `setindex!`/`push!`, which reverse mode cannot follow, and
        # now grows them with `vcat`. `richard_extrap` needed nothing — `LinearSolve` carries its own
        # rules and `sqrt` is `sqrt`.
        for f in (θ_ -> kr_sum(envf, θ_), θ_ -> solution_loss(envf, θ_))
            @test relerr(Zygote.gradient(f, θ)[1], ForwardDiff.gradient(f, θ)) < 1e-8
        end
    end

    @testset "kraken_jl end to end — mesh schedules" begin
        # Each of the loop's three exits reaches the return by a different route, and only the middle
        # one is what the default call does: `n_meshes = 1` returns before the extrapolation exists,
        # the convergence test breaks out partway, and `dont_break` runs every level. All three have
        # to be traceable, because the schedule is chosen by the environment, not by the caller.
        for kws in ((n_meshes=1,), (n_meshes=5,), (dont_break=true, n_meshes=3))
            f = θ_ -> kr_sum(pekeris, θ_; kws...)
            @test relerr(Zygote.gradient(f, θ0)[1], ForwardDiff.gradient(f, θ0)) < 1e-8
        end
    end

    @testset "interpolant rule: the primal value, exactly — $name" for (name, envf, θ) in (
        ("pekeris", pekeris, θ0), ("one_layer", onelayer, θ1), ("two_layer_slope", twolayer, θ2)
    )
        # `linear_interp_partials` recomputes the interpolated value instead of asking the primal for
        # it, so that this comparison exists: the partials are derivatives of the primal only if the
        # rule landed in the same interval. `===` rather than `≈` — the rule reproduces the arithmetic,
        # it does not approximate it.
        #
        # The queries that matter are the mesh points, because `get_z_vec` ends each layer's mesh
        # exactly on that layer's lower boundary, which is exactly where these environments place a
        # duplicated knot. The extra queries cover both extrapolation branches and both endpoints.
        env = envf(θ)
        for freq in (25.0, 100.0, 400.0)
            props = AcousticProblemProperties(env, freq)
            zn = reduce(vcat, props.zn_vec)
            queries = vcat(zn, env.layer_depth, [-5.0, 0.0, env.depth, env.depth + 10.0])
            @test all(z -> soundspeed(env.c, z) === first(Kraken.linear_interp_partials(env.c.f, z)), queries)
            @test all(z -> density(env.ρ, z) === first(Kraken.linear_interp_partials(env.ρ.f, z)), queries)
        end
    end

    @testset "interpolant rule: all five partial families" begin
        # The rule's content in isolation, against ForwardDiff on the same three inputs — the knot
        # depths, the values, and the query point. `two_layer_slope_env` is the one standard
        # environment where all three are simultaneously non-trivial.
        ssp = twolayer(θ2).c
        depths, values = -ssp.z, ssp.c        # the struct stores `z = -depth`; see the rule's comment
        rebuild(d, v) = SampledSSP(d, v)

        for z in (37.0, 110.0, 130.0)         # water column, first sediment layer, second
            @test Zygote.gradient(d -> soundspeed(rebuild(d, values), z), depths)[1] ≈
                ForwardDiff.gradient(d -> soundspeed(rebuild(d, values), z), depths) rtol = 1e-10
            @test Zygote.gradient(v -> soundspeed(rebuild(depths, v), z), values)[1] ≈
                ForwardDiff.gradient(v -> soundspeed(rebuild(depths, v), z), values) rtol = 1e-10
            @test Zygote.gradient(zz -> soundspeed(ssp, zz), z)[1] ≈
                ForwardDiff.derivative(zz -> soundspeed(ssp, zz), z) rtol = 1e-10
        end

        # Constant extrapolation: the nearer endpoint's *value* carries the whole cotangent, and the
        # knots and the query point carry none. Zero here is structural, not small.
        for z in (-5.0, 200.0)
            @test all(iszero, Zygote.gradient(d -> soundspeed(rebuild(d, values), z), depths)[1])
            @test iszero(Zygote.gradient(zz -> soundspeed(ssp, zz), z)[1])
            @test sum(Zygote.gradient(v -> soundspeed(rebuild(depths, v), z), values)[1]) ≈ 1.0
        end

        # The vector method is not the scalar one in a loop — it has its own rule — so it gets its own
        # check that it agrees with summing the scalar rule.
        zs = [17.0, 100.0, 110.0, 137.0]
        @test Zygote.gradient(d -> sum(soundspeed(rebuild(d, values), zs)), depths)[1] ≈
            sum(z -> Zygote.gradient(d -> soundspeed(rebuild(d, values), z), depths)[1], zs) rtol = 1e-12

        # The sign on the constructor. `z = -depth`, so a cotangent on the knots enters the struct
        # negated; getting this backwards is invisible in any environment whose profile is constant
        # within each layer, which is most of them.
        @test Zygote.gradient(d -> sum(SampledSSP(d, values).z), depths)[1] ≈ fill(-1.0, length(depths))
    end

    @testset "a profile that varies inside a layer whose depth is a parameter" begin
        # `two_layer_slope_env` has a sound-speed gradient *within* each sediment layer, so moving a
        # layer boundary moves the interpolant's knots as well as the mesh. This used to be the
        # milestone's one wrong answer: `DataInterpolations`' ChainRules rule holds an interpolant's
        # knots fixed, so Zygote dropped that term silently — `∂kr₁/∂h0` came out 15% low and `∂h1`,
        # `∂h2` had the wrong sign, with nothing raised. The rules in `src/kraken_ad.jl` replace it.
        #
        # Which side was wrong was settled against Fortran, not against ForwardDiff: finite
        # differences of `kraken.exe`'s own `Re(kᵣ)` across perturbed `.env` files put `∂kr₁/∂h0` at
        # 1.6344e-5, against ForwardDiff's 1.6342e-5 and the old Zygote's 1.3836e-5. The full table is
        # in the plan's Architecture Decisions; turning that one-off measurement into a standing test
        # is plan task 4.7, so what is pinned *here* is agreement with ForwardDiff, which the Fortran
        # run vouched for.
        #
        # The cause, isolated: the sound speed at a fixed depth inside the first sediment layer, as
        # that layer's upper boundary moves. It read 0.0 under Zygote and -1.5 under everything else.
        c_at(x) = soundspeed(twolayer([θ2[1:10]; x; θ2[12:13]]).c, 110.0)
        @test ForwardDiff.derivative(c_at, 100.0) ≈ central(c_at, 100.0) rtol = 1e-6
        @test Zygote.gradient(c_at, 100.0)[1] ≈ central(c_at, 100.0) rtol = 1e-6

        # And end to end, on every parameter — the three thicknesses included. Compared against the
        # gradient's own scale rather than entrywise: this environment's wavenumber gradient spans
        # seven orders of magnitude, and the smallest entries are below the precision either method
        # reaches.
        for f in (θ_ -> last(setup(twolayer, θ_)), θ_ -> mode_integral(twolayer, θ_), θ_ -> kr_sum(twolayer, θ_))
            g_zy = Zygote.gradient(f, θ2)[1]
            g_fw = ForwardDiff.gradient(f, θ2)
            @test relerr_norm(g_zy, g_fw) < 1e-8
            @test relerr_norm(g_zy[11:13], g_fw[11:13]) < 1e-8
        end
    end
end

### The multi-backend table -----------------------------------------------------------------------
#
# Everything above pins Zygote against ForwardDiff by hand. What follows makes that comparison a
# table — backends down one axis, environments and targets across the others — so that adding a
# backend is one entry in `BACKENDS` rather than a copy of a testset. `DifferentiationInterface`
# supplies the uniform `gradient(f, backend, θ)` that makes that possible; each of the four has its
# own calling convention otherwise (`Zygote.gradient` returns a tuple, Mooncake wants a prepared
# cache, `FiniteDiff` takes its step configuration in the type), and the table would be four
# near-copies of one loop.
#
# What each axis is for:
#
#   * **Backends.** Zygote and Mooncake reach the same rules by different routes, and the route is
#     the thing being tested — a rule can be correct and still be reached wrongly. Zygote consumes
#     `ChainRulesCore.rrule`s directly. Mooncake does not: it derives its own rules from the IR
#     unless a signature is declared a primitive, so `ext/KrakenMooncakeExt.jl` declares each of
#     ours and translates between the two tangent representations. Without that extension Mooncake
#     does not merely run slowly, it stops with a `bitcast` error inside `range`.
#   * **References.** ForwardDiff differentiates the same code by an unrelated mechanism, so it
#     catches an error in a rule's content. Finite differences differentiate no code at all, so they
#     catch what the two AD systems could get wrong together — a rule that is internally consistent
#     and still not the derivative of the primal. That is exactly the shape of the interpolant gap
#     4.5 closed, and finite differences (of Fortran, there) are what settled it.
#
# Task 4.7 adds the third kind of reference: Fortran's own numbers, finite-differenced across
# perturbed `.env` files, which is the only oracle here that shares no code with Kraken.jl at all.

"""
Munk profile as a function of the five parameters that define it: `θ = [c_ref, ε, z_axis, cb, ρb]`.

`munk_env()` takes no arguments — the canonical profile is one fixed curve — so a differentiable
parameterization has to be introduced here rather than imported. At `θ3` this reproduces
`munk_env()` exactly, which the first testset below asserts.

The knot depths are constants, and they are built by comprehension rather than with `range`: a
`range` of `Float64`s is a `StepRangeLen` over `Base.TwicePrecision`, and Zygote has no adjoint for
that constructor even when — as here — its arguments carry no derivative. `src/kraken_ad.jl` meets
the same wall in `layer_mesh`, where the depths *are* parameters, and answers it with a rule.
"""
function munk(θ; depth=5000.0, nsamples=51)
    c_ref, ϵ, z_axis, cb, ρb = θ
    ρ0 = 1000.0
    zs = [depth * (k - 1) / (nsamples - 1) for k in 1:nsamples]
    profile(z) = (ẑ=2 * (z - z_axis) / z_axis; c_ref * (1 + ϵ * (ẑ - 1 + exp(-ẑ))))
    ssp = reduce(vcat, [[z profile(z) 0.0 ρ0 0.0 0.0] for z in zs])
    return UnderwaterEnv(ssp, [0.0 0.0 depth], [0.0 343.0 0.0 0.00121 0.0 0.0; depth cb 0.0 ρb 0.0 0.0])
end

θ3 = [1500.0, 0.00737, 1300.0, 1600.0, 1500.0]

# Zygote reads the rules in `src/kraken_ad.jl` directly; Mooncake reads them through
# `ext/KrakenMooncakeExt.jl`. Adding a third is one line here.
const BACKENDS = [("Zygote", AutoZygote()), ("Mooncake", AutoMooncake())]

# `munk` runs at 10 Hz rather than 100: 5000 m of water at 100 Hz is a 6600-point mesh with hundreds
# of modes, and the finite-difference reference alone would then need ten of those solves. At 10 Hz
# it is 625 points and 20 modes, which still exercises every path in the table.
#
# The last column lists the parameters a *finite difference* can measure, and it is a shorter list
# than θ for a reason that is a documented property of the solver rather than a limitation of the
# test. `Nz_vec = max(10, ceil(h / Δz))` with `Δz = c_b / (20 f)`, so the mesh point count is a step
# function of `cb` and of every layer thickness — and of nothing else. Task 4.4 settled what the
# gradient means there: it is the derivative at a *fixed* mesh schedule, which is the standard
# treatment and the only one that is well defined, but it makes the primal genuinely discontinuous
# in those parameters. A central difference straddles the jumps: on `∂(∫ψ dz)/∂depth` for Pekeris it
# reads 5.29 at a 1e-6 relative step, 1.82 at 1e-5 and 1.47 at 1e-4, wandering toward — and never
# reaching — the 1.4356 that ForwardDiff, Zygote and Mooncake all agree on. Differencing the
# remaining parameters is a real check; differencing those is measuring the mesh.
const ENVIRONMENTS = [
    ("pekeris", pekeris, θ0, 100.0, [1, 3, 4]),          # θ0 = [c0, cb, ρ0, ρb, depth]
    ("one_layer", onelayer, θ1, 100.0, [1, 2, 4, 5, 6]),  # θ1 = [c0, c1, cb, ρ0, ρ1, ρb, h0, h1]
    ("munk", munk, θ3, 10.0, [1, 2, 3, 5]),               # θ3 = [c_ref, ε, z_axis, cb, ρb]
]

# The three things a user differentiates: the wavenumbers, a mode shape, and a scalar loss over both
# at once. All three go through `kraken_jl` or its two rule-bearing pieces, so between them they
# cover every seam in `src/kraken_core.jl`'s list.
const TARGETS = [
    ("wavenumbers", (envf, θ, freq) -> kr_sum(envf, θ; freq=freq)),
    ("mode shape", (envf, θ, freq) -> mode_integral(envf, θ; freq=freq)),
    ("scalar loss", (envf, θ, freq) -> solution_loss(envf, θ; freq=freq)),
]

@testset "Reverse-mode AD across backends" begin
    @testset "the parameterized Munk profile reproduces munk_env()" begin
        # The table below is only a test of the Munk environment if `munk(θ3)` is that environment.
        reference = UnderwaterEnv(munk_env()...)
        built = munk(θ3)
        @test built.c.c ≈ reference.c.c
        @test built.c.z == reference.c.z
        @test (built.cb, built.ρb, built.depth) == (reference.cb, reference.ρb, reference.depth)
    end

    @testset "$env_name / $target_name" for (env_name, envf, θ, freq, differenceable) in ENVIRONMENTS,
        (target_name, target) in TARGETS

        f = θ_ -> target(envf, θ_, freq)

        # ForwardDiff at 1e-7 rather than the 1e-8 the hand-written testsets above use, because at
        # 20 modes ForwardDiff is the side that moves: it differentiates the arithmetic of the root
        # solver's last iterate, so its answer carries the iteration's truncation (4.2's finding,
        # sharpened in 4.3). It shows only on `munk`, and only on the wavenumber sum — 8.9e-9 there
        # against 1e-12 everywhere else in this table.
        g_forward = DifferentiationInterface.gradient(f, AutoForwardDiff(), θ)

        # Central differences, and only on the parameters that do not move the mesh — see the
        # comment on `ENVIRONMENTS`. FiniteDiff's default is *forward* differences, which are not
        # accurate enough to be a reference here at any parameter (3.7% out on the Pekeris
        # wavenumber sum).
        g_finite = DifferentiationInterface.gradient(f, AutoFiniteDiff(; fdtype=Val(:central)), θ)

        gradients = map(BACKENDS) do (backend_name, backend)
            g = DifferentiationInterface.gradient(f, backend, θ)
            @testset "$backend_name" begin
                @test length(g) == length(θ)
                @test all(isfinite, g)
                @test relerr_norm(g, g_forward) < 1e-7
                @test relerr_norm(g[differenceable], g_finite[differenceable]) < 1e-5
            end
            return g
        end

        # The backends evaluate the *same* rules, so they should agree far more closely than either
        # agrees with a reference — anything else means a backend is reaching the rules differently,
        # which is exactly the failure mode `ext/KrakenMooncakeExt.jl` could introduce. Observed:
        # 3e-14 and below.
        @testset "backends agree with each other" begin
            for k in 2:length(gradients)
                @test relerr_norm(gradients[k], gradients[1]) < 1e-11
            end
        end
    end
end

### The rules on their own, through ChainRulesTestUtils --------------------------------------------
#
# Every comparison so far is a gradient of a *composition*, where a wrong term in one rule can be
# cancelled by a wrong term in another or simply never reached by the cotangent that happens to flow.
# `test_rrule` asks the narrower question: with a random cotangent, does this one pullback reproduce
# finite differences of this one primal, over every argument at once? It checks arguments the solver
# never varies independently — `props.zn_vec` against `cache.a_vec`, say — and it checks the tangent
# *types*, which is how a rule acquires a `Tangent` field that no backend ever accumulates.
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

@testset "ChainRulesTestUtils" begin
    # `layer_mesh` is the whole rule: `collect(range(z_top, z_bot, Nz))`, three arguments, one of
    # them an `Int` that must come back `NoTangent`. Nothing else in this file tests it in isolation.
    @testset "layer_mesh" begin
        test_rrule(Kraken.layer_mesh, 0.5, 100.0, 12; check_inferred=false)
    end

    # `atol = 1e-8` and not smaller because of one specific entry: `∂kr/∂ρb` is ~1.3e-6 on
    # `one_layer_env`, three decades below the largest cotangent in the same call, and finite
    # differences resolve it to about 2e-9 absolute. Everything else clears `rtol = 1e-6` on its own.
    @testset "solve_for_kr — $name" for (name, envf, θ) in (("pekeris", pekeris, θ0), ("one_layer", onelayer, θ1))
        # 25 Hz keeps the mesh at a few dozen points: `test_rrule` finite-differences *every*
        # coordinate of every argument, so the cost is quadratic in the mesh, not linear.
        env, props, cache, _ = setup(envf, θ; freq=25.0)
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
    # content" testset above, which covers every partial family entry by entry.
end
