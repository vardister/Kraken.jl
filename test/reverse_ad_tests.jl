using Test
using Kraken
using ForwardDiff
using Zygote
using Mooncake
using FiniteDiff
using DifferentiationInterface
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

function env_with(env; cb=env.cb, ρb=env.ρb, αb=env.αb)
    T = promote_type(typeof(cb), typeof(ρb), typeof(αb))
    return UnderwaterEnv{typeof(env.c),typeof(env.ρ),T,typeof(env.α),typeof(env.top_bc),typeof(env.bottom_bc)}(
        env.c,
        env.ρ,
        T(cb),
        T(ρb),
        T.(env.h_vec),
        T.(env.layer_depth),
        T(env.depth),
        env.α,
        T(αb),
        env.atten_units,
        env.top_bc,
        env.bottom_bc,
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

    @testset "the rule is tolerance-independent, and ForwardDiff now uses it too" begin
        # 250 Hz, mode 1. Until plan task 6.11 ForwardDiff differentiated the arithmetic of the *last
        # ITP iterate*, and here it was 1.3% off `∂kr₁/∂cb`: this testset asserted that the rule was
        # the closer of the two to central differences. 6.11 gave ForwardDiff the same implicit-function
        # derivative (`ext/KrakenForwardDiffExt.jl`), so the two now agree to rounding.
        f = θ -> kr_at(θ; freq=250.0, mode=1)
        g_zy = Zygote.gradient(f, θ0)[1]
        g_fw = ForwardDiff.gradient(f, θ0)
        g_fd = map(1:length(θ0)) do i
            h = 1e-6 * abs(θ0[i])
            step = [j == i ? h : 0.0 for j in 1:length(θ0)]
            (f(θ0 .+ step) - f(θ0 .- step)) / (2h)
        end

        @test relerr(g_zy, g_fd) < 1e-4
        @test relerr(g_fw, g_zy) < 1e-9

        # Tightening the root tolerance by four orders of magnitude must move neither answer.
        tight = θ -> kr_at(θ; freq=250.0, mode=1, tol=1e-10)
        @test relerr(g_zy, Zygote.gradient(tight, θ0)[1]) < 1e-9
        @test relerr(g_fw, ForwardDiff.gradient(tight, θ0)) < 1e-9
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

    # Every environment above spells its interfaces `z0`, `z0 + eps(z0)`. The documented alternative is
    # two rows at exactly `z0`, and before plan task 6.9 the primal and this rule disagreed there: the
    # primal read the water column's last mesh point from the medium *below*, while the rule picked the
    # interval ending at the knot. So Zygote differentiated a function the solver did not compute, and
    # ∂kᵣ/∂c_water on this shelf was 5% wrong on one mesh and 15% wrong through `kraken_jl`. 6.9 fixed the
    # primal by separating equal depths when the profile is built. These tests check that the rule and
    # the primal now agree for the same-depth spelling too.
    shelf(cw) = UnderwaterEnv(
        [
            0.0 cw 0.0 1030.0 0.0 0.0
            74.4 cw 0.0 1030.0 0.0 0.0
            74.4 1445.0 0.0 1600.0 0.0 0.0
            85.4 1555.0 0.0 1600.0 0.0 0.0
        ],
        [0.0 0.0 74.4; 0.0 0.0 85.4],
        [0.0 343.0 0.0 1.21 0.0 0.0; 85.4 1830.0 0.0 2000.0 0.0 0.0],
    )

    @testset "interpolant rule: the primal value, exactly — interface at one depth" begin
        env = shelf(1468.3)
        for freq in (25.0, 100.0)
            props = AcousticProblemProperties(env, freq)
            queries = vcat(reduce(vcat, props.zn_vec), env.layer_depth, [0.0, env.depth])
            @test all(z -> soundspeed(env.c, z) === first(Kraken.linear_interp_partials(env.c.f, z)), queries)
            @test all(z -> density(env.ρ, z) === first(Kraken.linear_interp_partials(env.ρ.f, z)), queries)
        end
        # The same holds for a profile built directly, without going through `UnderwaterEnv`.
        ssp = SampledSSP([0.0, 50.0, 50.0, 100.0], [1500.0, 1490.0, 1600.0, 1650.0])
        @test soundspeed(ssp, 50.0) === first(Kraken.linear_interp_partials(ssp.f, 50.0)) === 1490.0
    end

    @testset "a gradient across an interface at one depth" begin
        # ∂kᵣ₁/∂c_water, the parameter the old disagreement hit hardest. Central differences use no AD
        # at all. On one mesh the bound is set by the root tolerance: `kraken_jl` solves the coarse
        # mesh to `abstol = 1e-6`, which leaves the rule 4.5e-7 from the differences (the 4.4 note).
        for kws in ((n_meshes=1,), (;))
            f = cw -> real(kraken_jl(shelf(cw), 50.0; kws...).kr[1])
            fd = (f(1468.3 + 1e-3) - f(1468.3 - 1e-3)) / 2e-3
            @test abs(Zygote.gradient(f, 1468.3)[1] / fd - 1) < 1e-6
            @test abs(ForwardDiff.derivative(f, 1468.3) / fd - 1) < 1e-6
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

### Milestone 5: differentiating a lossy solve -----------------------------------------------------
#
# `modal_attenuation` needed no rule of its own — it is traced arithmetic sitting downstream of the
# two functions that do carry rules. These tests are what makes that claim checkable rather than
# hopeful, and they pin the one thing attenuation genuinely changes: `kr` is now complex, so what
# counts as a differentiable loss is narrower than it was.

# `θ = [c0, cb, ρ0, ρb, depth, α0, αb]` — the Pekeris parameters plus both attenuations. Built the
# same non-mutating way as every other environment here; `pekeris_env`'s matrix literals promote, so
# a `Dual` or a tracked value in `α0`/`αb` reaches the solver intact.
pekeris_lossy(θ) = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5], α0=θ[6], αb=θ[7])...)

# 0.1 dB/λ in the water and 0.3 in the seabed: lossy enough that `Im(kᵣ)` is far above roundoff,
# weak enough that the near-cutoff mode is not in the saturating regime documented in
# `test/README.md` (which would make the finite-difference reference the inaccurate side).
θ_lossy = [1500.0, 1600.0, 1000.0, 1500.0, 100.0, 0.1, 0.3]

lossy_kr(θ; freq=100.0) = kraken_jl(pekeris_lossy(θ), freq; TOL...).kr

@testset "Reverse-mode AD through attenuation" begin
    @testset "the environment is actually lossy and actually complex" begin
        # If this fails, everything below is silently testing the lossless solver.
        env = pekeris_lossy(θ_lossy)
        @test is_lossy(env)
        @test eltype(lossy_kr(θ_lossy)) === ComplexF64
        @test all(imag.(lossy_kr(θ_lossy)) .< 0)
        # ...and the lossless corner of the same parameterization is not complex.
        @test eltype(lossy_kr([1500.0, 1600.0, 1000.0, 1500.0, 100.0, 0.0, 0.0])) === Float64
    end

    @testset "d(imag kᵣ)/dα — $name" for (name, target) in (
        ("summed", θ -> sum(imag, lossy_kr(θ))),
        ("first mode", θ -> imag(lossy_kr(θ)[1])),
        ("last mode", θ -> imag(lossy_kr(θ)[end])),
    )
        # The plan's acceptance criterion. Both attenuation entries must carry a real gradient —
        # a silently-zero ∂/∂α is precisely what a missing interpolant rule would produce.
        g_forward = ForwardDiff.gradient(target, θ_lossy)
        g_zygote = Zygote.gradient(target, θ_lossy)[1]

        @test relerr_norm(g_zygote, g_forward) < 1e-8
        @test g_forward[6] != 0                      # ∂/∂α0, the water column
        @test g_forward[7] != 0                      # ∂/∂αb, the half-space
        @test all(isfinite, g_zygote)

        # More loss means more attenuation, and `Im(kᵣ)` is negative — so both derivatives are
        # negative. Sign errors here are the failure this whole milestone is exposed to.
        @test g_forward[6] < 0
        @test g_forward[7] < 0
    end

    @testset "against central differences, which use no AD at all" begin
        # ForwardDiff and Zygote agreeing proves the rules are consistent, not that they are right —
        # both read the same primal. `α0` and `αb` move no mesh (unlike `cb` and `depth`; see the
        # note on ENVIRONMENTS), so a central difference in them is a clean independent reference.
        target = θ -> sum(imag, lossy_kr(θ))
        g_forward = ForwardDiff.gradient(target, θ_lossy)
        for k in (6, 7)
            fd = central(x -> target([θ_lossy[1:(k - 1)]; x; θ_lossy[(k + 1):end]]), θ_lossy[k])
            @test isapprox(fd, g_forward[k]; rtol=1e-5)
        end
    end

    @testset "real-valued losses over a complex kᵣ — $name" for (name, target) in (
        ("real part", θ -> sum(real, lossy_kr(θ))),
        ("abs2", θ -> sum(abs2, lossy_kr(θ))),
        ("modulus", θ -> sum(abs, lossy_kr(θ))),
    )
        # The documented convention: real parameters in, real loss out, complex only in between.
        # Under that restriction the gradient is unambiguous and every backend agrees.
        g_forward = ForwardDiff.gradient(target, θ_lossy)
        @test relerr_norm(Zygote.gradient(target, θ_lossy)[1], g_forward) < 1e-8
        @test all(isfinite, g_forward)
    end

    @testset "a complex-valued loss is refused, not guessed" begin
        # There is no single real gradient of a complex output, and Zygote says so rather than
        # silently picking a convention. Pinned because the *absence* of a convention is the
        # decision — see the note at the end of src/kraken_ad.jl.
        err = try
            Zygote.gradient(θ -> lossy_kr(θ)[1], θ_lossy)
            nothing
        catch e
            e
        end
        @test err !== nothing
        @test occursin("complex", lowercase(sprint(showerror, err)))

        # Taking a real part first is the fix, and it works.
        @test Zygote.gradient(θ -> real(lossy_kr(θ)[1]), θ_lossy)[1] isa Vector{Float64}
    end

    @testset "at α = 0 the two modes take opposite sides of the kink" begin
        # `imag(kr)` is identically zero for α ≤ 0 and linear above it, so α = 0 is a genuine kink
        # with no two-sided derivative. Both backends silently return a *different* one-sided one,
        # and that disagreement is the thing worth pinning — it is exactly what gets mistaken for a
        # broken rule.
        #
        # ForwardDiff: `iszero` on a `Dual` sees the seed, so `is_lossy` takes the lossy branch and
        #              the whole attenuation path is evaluated -> right-hand derivative.
        # Zygote:      `is_lossy` runs on the primal Float64s, returns false, and the attenuation
        #              path never reaches the tape -> left-hand derivative, which is zero.
        target = θ -> sum(imag, lossy_kr(θ))
        θ_zero = [1500.0, 1600.0, 1000.0, 1500.0, 100.0, 0.0, 0.0]

        g_forward_0 = ForwardDiff.gradient(target, θ_zero)
        g_zygote_0 = Zygote.gradient(target, θ_zero)[1]

        @test target(θ_zero) == 0.0                    # the primal really is the lossless one
        @test g_forward_0[6] < -1e-3                   # right-hand derivative, and substantial
        @test g_forward_0[7] < -1e-6
        @test all(iszero, g_zygote_0)                  # left-hand derivative
        # Everything that is not an attenuation is zero on both sides — a lossless solve has no
        # imaginary part to be sensitive to anything.
        @test all(iszero, g_forward_0[1:5])

        # Step off the kink by a hair and the disagreement vanishes entirely: same derivative to
        # fifteen digits. That is what makes the split above a property of the point, not of either
        # backend's accuracy.
        θ_eps = [1500.0, 1600.0, 1000.0, 1500.0, 100.0, 1e-6, 1e-6]
        g_forward_eps = ForwardDiff.gradient(target, θ_eps)
        g_zygote_eps = Zygote.gradient(target, θ_eps)[1]
        @test relerr_norm(g_zygote_eps, g_forward_eps) < 1e-12
        # ...and forward mode's answer *at* the kink is that same right-hand derivative.
        @test isapprox(g_forward_0[6], g_forward_eps[6]; rtol=1e-6)
        @test isapprox(g_forward_0[7], g_forward_eps[7]; rtol=1e-6)
    end

    @testset "Zygote traces the lossy solve" begin
        # The lossless path is covered across every backend by the table above. Here only Zygote is
        # asserted, because Mooncake cannot do it — see the next testset.
        target = θ -> sum(imag, lossy_kr(θ))
        g_forward = DifferentiationInterface.gradient(target, AutoForwardDiff(), θ_lossy)
        g = DifferentiationInterface.gradient(target, AutoZygote(), θ_lossy)
        @test length(g) == length(θ_lossy)
        @test all(isfinite, g)
        @test relerr_norm(g, g_forward) < 1e-7
    end

    # Opt-in, and the only testset in this file that is. It costs ~45 s — two cold Mooncake
    # derivations, one of which is thrown away when it fails — to establish that a limitation
    # *still* holds. That is a tripwire for the day Mooncake gains complex support: worth having on
    # every CI run, worth nothing on the twentieth local run of the hour.
    # `KRAKEN_RUN_KNOWN_LIMITATIONS=true` runs it, and `.github/workflows/CI.yml` sets it.
    #
    # Note what goes with it when it is skipped: the lossless-corner check at the bottom, which is
    # the only place Mooncake meets *this* parameterization. Mooncake over lossless environments is
    # covered thoroughly by the backend table above, so that is a duplicate rather than a gap.
    #
    # The guard is outside the `@testset` rather than an early return inside it: a testset body is a
    # `let` block, not a function, so `return` there would escape the enclosing scope.
    if get(ENV, "KRAKEN_RUN_KNOWN_LIMITATIONS", "false") != "true"
        @info "SKIPPING the Mooncake known-limitation tripwire. Set KRAKEN_RUN_KNOWN_LIMITATIONS=true to run it."
    else
        @testset "Mooncake cannot trace the complex path (known limitation)" begin
            # `add_attenuation` introduces complex arithmetic — `sqrt` of a `Complex` — and Mooncake
            # refuses it:
            #
            #   ArgumentError: It is not permissible to bitcast to a differentiable type during AD ...
            #
            # The offending call is inside Mooncake's own handling of `Complex`, not in
            # src/kraken_ad.jl, so no rule here would fix it. Pinned as a *specific* failure rather than
            # left out, so that Mooncake gaining complex support shows up as this test starting to pass
            # rather than as nobody noticing.
            target = θ -> sum(imag, lossy_kr(θ))
            err = try
                DifferentiationInterface.gradient(target, AutoMooncake(), θ_lossy)
                nothing
            catch e
                e
            end
            @test err !== nothing
            @test occursin("bitcast", lowercase(sprint(showerror, err)))

            # Mooncake is fine on the same environment as long as the loss is real-valued throughout —
            # i.e. the lossless corner. This is what confines the limitation to the complex path.
            θ_lossless = [1500.0, 1600.0, 1000.0, 1500.0, 100.0, 0.0, 0.0]
            real_target = θ -> sum(real, lossy_kr(θ))
            g_mooncake = DifferentiationInterface.gradient(real_target, AutoMooncake(), θ_lossless)
            g_forward = DifferentiationInterface.gradient(real_target, AutoForwardDiff(), θ_lossless)
            @test relerr_norm(g_mooncake, g_forward) < 1e-7
        end
    end

    @testset "a lossy sediment layer, not just a lossy half-space" begin
        # `one_layer_env`'s α1 goes through the *volume* integral rather than the half-space term, so
        # this covers the other branch of `modal_attenuation` — and it is the physically common case.
        onelayer_lossy(θ) = UnderwaterEnv(
            one_layer_env(; c0=θ[1], c1=θ[2], cb=θ[3], ρ0=θ[4], ρ1=θ[5], ρb=θ[6], h0=θ[7], h1=θ[8], α1=θ[9])...
        )
        θ1_lossy = [1500.0, 1550.0, 1600.0, 1000.0, 1500.0, 2000.0, 100.0, 20.0, 0.4]
        target = θ -> sum(imag, kraken_jl(onelayer_lossy(θ), 100.0; TOL...).kr)

        g_forward = ForwardDiff.gradient(target, θ1_lossy)
        @test relerr_norm(Zygote.gradient(target, θ1_lossy)[1], g_forward) < 1e-8
        @test g_forward[9] < 0
        # The sediment attenuation is the dominant sensitivity of the modal loss, by a wide margin.
        @test abs(g_forward[9]) > 10 * maximum(abs, g_forward[1:8])
    end
end

### The rules on their own, through ChainRulesTestUtils --------------------------------------------
#
# Moved to `rrule_tests.jl`. `test_rrule` checks a pullback against finite differences over every
# argument at once, which is a different question from the compositions above — and its
# `solve_for_kr` case costs 373 s cold against 44 s warm, essentially all of it Julia compiling
# `test_rrule` for the solver's concrete signature. Its own file makes that block skippable through
# `KRAKEN_SKIP_TESTS` while iterating here. Nothing was dropped.

### Milestone 6: every boundary condition and every interpolation mode ----------------------------
#
# Each boundary condition rewrites the first or last rows of the tridiagonal system both rules
# differentiate through, and each interpolation mode changes the map from profile parameters to those
# rows. Everything above runs on the default pair (pressure-release over a half-space, c-linear), so
# without this section a non-default configuration could lose its gradient silently. One already had:
# a rigid top over a rigid bottom returned an all-NaN gradient until the eigenvector pullback learned
# to survive an exactly-zero pivot (6.6).

"""
One medium over `θ = [c₁, c₂, c₃, c₄, z₂, cb, ρb, depth]`: sound-speed knots at `0`, `z₂`,
`0.7·depth` and `depth`, so the profile varies inside the water column and one knot depth is a
parameter. A single medium is what a `:cubic_spline` profile needs (6.3), so all three interpolation
modes run on the same environment.
"""
function m6env(θ; top_bc=PressureRelease(), bottom_bc=AcousticHalfspace(), ssp_interp=:c_linear)
    c1, c2, c3, c4, z2, cb, ρb, depth = θ
    ssp = [
        0.0 c1 0.0 1000.0 0.0 0.0
        z2 c2 0.0 1000.0 0.0 0.0
        0.7*depth c3 0.0 1000.0 0.0 0.0
        depth c4 0.0 1000.0 0.0 0.0
    ]
    sspHS = [0.0 343.0 0.0 0.00121 0.0 0.0; depth cb 0.0 ρb 0.0 0.0]
    return UnderwaterEnv(ssp, [0.0 0.0 depth], sspHS; top_bc, bottom_bc, ssp_interp)
end

θ6 = [1500.0, 1490.0, 1510.0, 1520.0, 40.0, 1700.0, 1800.0, 100.0]

# What a central difference can measure: the four sound speeds, the interior knot depth and the
# half-space density. `cb` and `depth` set the mesh point count (see the comment on `ENVIRONMENTS`).
const M6_DIFFERENCEABLE = [1, 2, 3, 4, 5, 7]

# 50 Hz: 3 to 7 modes depending on the boundaries, none of them near a cutoff.
m6loss(θ; kws...) = solution_loss(θ_ -> m6env(θ_; kws...), θ; freq=50.0)

function central_gradient(f, θ, idx; h=1e-6)
    return map(idx) do i
        step = [j == i ? h * abs(θ[i]) : 0.0 for j in eachindex(θ)]
        return (f(θ .+ step) - f(θ .- step)) / (2 * step[i])
    end
end

# One configuration per new option also goes through Mooncake, which reaches the same rules through
# `ext/KrakenMooncakeExt.jl`. Every configuration would recompile Mooncake's rules a dozen times.
const M6_MOONCAKE = [
    (RigidBoundary(), AcousticHalfspace(), :c_linear),
    (PressureRelease(), RigidBoundary(), :c_linear),
    (PressureRelease(), PressureRelease(), :c_linear),
    (PressureRelease(), AcousticHalfspace(), :n2_linear),
]

@testset "Reverse-mode AD across boundary conditions and SSP interpolation" begin
    @testset "$top over $bottom, $interp" for top in (PressureRelease(), RigidBoundary()),
        bottom in (AcousticHalfspace(), RigidBoundary(), PressureRelease()),
        interp in (:c_linear, :n2_linear)

        f = θ_ -> m6loss(θ_; top_bc=top, bottom_bc=bottom, ssp_interp=interp)
        @test !isempty(kraken_jl(m6env(θ6; top_bc=top, bottom_bc=bottom, ssp_interp=interp), 50.0).kr)

        g_zygote = Zygote.gradient(f, θ6)[1]
        g_forward = ForwardDiff.gradient(f, θ6)
        @test all(isfinite, g_zygote)

        # Central differences share no code with either AD system, so they are the arbiter. Measured
        # 1.6e-8 to 3.2e-7 across all twelve configurations.
        @test relerr_norm(g_zygote[M6_DIFFERENCEABLE], central_gradient(f, θ6, M6_DIFFERENCEABLE)) < 1e-6

        # Only a half-space bottom has a `g` term; over any other `cb` and `ρb` reach nothing.
        if !(bottom isa AcousticHalfspace)
            @test g_zygote[6] == 0 && g_zygote[7] == 0
        end

        if top isa PressureRelease
            # Measured ≤ 5e-12.
            @test relerr_norm(g_zygote, g_forward) < 1e-8
        else
            # Under a rigid top ForwardDiff stops being a reference. It differentiates the root
            # solver's last iterate rather than the root (4.2), and under a rigid top that error is
            # erratic rather than ~1e-12: 3e-13 over a half-space at 50 Hz but 6e-4 at 100 Hz, 3e-6
            # over a rigid bottom, 1.4e-3 over a vacuum. Per mode on a refinement mesh, the rule matches central
            # differences to every printed digit while ForwardDiff is 0.3% out on some modes, and a
            # tighter solver tolerance moves that error without removing it. The rule is right, as
            # the central-difference check above says; this bound only catches a gross break.
            @test relerr_norm(g_zygote, g_forward) < 1e-2
        end

        if (top, bottom, interp) in M6_MOONCAKE
            g_mooncake = DifferentiationInterface.gradient(f, AutoMooncake(), θ6)
            @test relerr_norm(g_mooncake, g_zygote) < 1e-11
        end
    end

    @testset "a cubic spline: reverse mode refuses, ForwardDiff goes through — $top over $bottom" for top in (
            PressureRelease(), RigidBoundary()
        ),
        bottom in (AcousticHalfspace(), RigidBoundary(), PressureRelease())

        # 6.3 left reverse mode without a rule for a spline's coefficients, so it throws rather than
        # returning the linear interpolant's derivative. ForwardDiff is the way it names, so that is
        # what gets checked, against central differences. The rigid-top bound is looser for the
        # reason given above; measured ≤ 6.3e-7 under a pressure-release top and ≤ 5.7e-5 under a
        # rigid one.
        f = θ_ -> m6loss(θ_; top_bc=top, bottom_bc=bottom, ssp_interp=:cubic_spline)
        @test_throws ArgumentError Zygote.gradient(f, θ6)
        g_forward = ForwardDiff.gradient(f, θ6)
        @test all(isfinite, g_forward)
        @test relerr_norm(g_forward[M6_DIFFERENCEABLE], central_gradient(f, θ6, M6_DIFFERENCEABLE)) <
            (top isa PressureRelease ? 1e-5 : 1e-3)
    end
end

### Milestone 7: differentiating the pressure field ------------------------------------------------
#
# The field introduces no solver of its own. `acoustic_field` is arithmetic over `sol.kr` and
# `sol.modes`, and both of those already carry rules, so 7.5's question is not *which rule does the
# field need* but *does anything in `src/kraken_field.jl` break the tape*. Three things could have:
# the `[f(i, j) for i in …, j in …]` comprehensions, the generator `mode_sum` consumes, and the
# `searchsortedlast` that picks an interpolation bracket out of knots that are themselves parameters.
#
# None of them do, and that is worth recording because the plan predicted otherwise. 7.5 named
# mode-shape interpolation as the likely obstruction, on the grounds that `DataInterpolations`
# constructors are not reverse-mode friendly. Task 7.1 had already written that interpolation as a
# dozen lines of linear weights instead of a `DataInterpolations` object — chosen to match
# `field.exe`'s `Weight`, not for AD — and that is the whole reason there was nothing to fix. The
# `searchsortedlast` is admissible for the same reason `bisection` is: it returns an integer, so it
# carries no cotangent, and the bracket it picks is locally constant in the parameters.
#
# Broadband synthesis was a different story and did need fixing: `broadband_field` assembled its
# output by writing planes into an `Array{ComplexF64,3}(undef, …)`, which reverse mode refuses
# outright and which strips a `Dual` before forward mode ever sees it. Both it and the three-argument
# `synthesize_pulse` now stack by comprehension, and the last testset here is what says so.

"""
`θ = [c₁ … c₁₀, cb]` — a ten-point sound-speed profile over a half-space whose sound speed is the
eleventh parameter. This is 7.5's acceptance criterion verbatim, and the pairing is deliberate: the
profile reaches the field only through the wavenumbers and mode shapes, while `cb` also enters
[`Kraken.halfspace_decay`](@ref) and the half-space term of the normalization directly, so a field
gradient that was silently dropping one of the two routes would still look plausible on the other.

The knot depths are a comprehension rather than a `range`, for the reason `munk` above gives.
"""
function field_env(θ; depth=100.0)
    n = length(θ) - 1
    cvec = θ[1:n]
    z = [depth * (k - 1) / (n - 1) for k in 1:n]
    ssp = hcat(z, cvec, zero(cvec), fill(1000.0, n), zero(cvec), zero(cvec))
    sspHS = [0.0 343.0 0.0 0.00121 0.0 0.0; depth θ[end] 0.0 1500.0 0.0 0.0]
    return UnderwaterEnv(ssp, [0.0 0.0 depth], sspHS)
end

"""
`field_env`'s parameterization with a water-column attenuation spliced in ahead of `cb`:
`θ = [c₁ … c₁₀, α₀, cb]`, over a seabed that is lossy at a fixed 0.3 dB/λ.

A separate builder rather than a keyword on `field_env`, because the parameter *vector* is a
different length — and a named function rather than a closure in the testset, because the version
written inline formatted into a `(a=…; b=…)` block whose parse is worth not having to think about.
"""
function lossy_field_env(θ; depth=100.0)
    n = length(θ) - 2
    cvec = θ[1:n]
    z = [depth * (k - 1) / (n - 1) for k in 1:n]
    ssp = hcat(z, cvec, zero(cvec), fill(1000.0, n), fill(θ[end - 1], n), zero(cvec))
    sspHS = [0.0 343.0 0.0 0.00121 0.0 0.0; depth θ[end] 0.0 1500.0 0.3 0.0]
    return UnderwaterEnv(ssp, [0.0 0.0 depth], sspHS)
end

function lossy_field_loss(θ; freq=100.0)
    sol = kraken_jl(lossy_field_env(θ), freq; TOL...)
    return sum(abs2, transmission_loss(sol, FIELD_RANGES, FIELD_ZS, 50.0; nmodes=3))
end

# A downward-refracting profile with structure in it, not a constant: over an isovelocity column the
# mode shapes are sinusoids whose sensitivity to each knot is nearly the same number, and a gradient
# that mixed the knots up would still pass.
θ7 = [[1500.0 + 20.0 * (k - 1) / 9 - 6.0 * sinpi(2 * (k - 1) / 9) for k in 1:10]; 1600.0]

# The field geometry every target below shares. Five ranges rather than one so the gradient has to
# travel the range axis, and a receiver both in the water column and below the bottom so the
# half-space tail of `mode_amplitude` is on the tape too.
const FIELD_RANGES = collect(2_000.0:1_000.0:6_000.0)
const FIELD_ZS = 25.0
const FIELD_ZR = [50.0, 120.0]

# `nmodes = 3` and not "all of them", because the mode *count* is a step function of θ — the same
# integer-valued discontinuity documented for the mesh schedule on `ENVIRONMENTS` above. Summing a
# fixed number of modes is the derivative that is well defined; summing however many the solve
# returned is not. At θ7 and 100 Hz the solve finds 5, so 3 is comfortably inside.
field_solution(θ; freq=100.0) = kraken_jl(field_env(θ), freq; TOL...)

const FIELD_TARGETS = [
    ("mode amplitudes", θ -> sum(abs2, mode_amplitudes(field_solution(θ), FIELD_ZR; nmodes=3))),
    ("pressure", θ -> sum(abs2, acoustic_field(field_solution(θ), FIELD_RANGES, FIELD_ZS, FIELD_ZR; nmodes=3))),
    (
        "transmission loss",
        θ -> sum(abs2, transmission_loss(field_solution(θ), FIELD_RANGES, FIELD_ZS, FIELD_ZR; nmodes=3)),
    ),
    (
        "incoherent TL",
        θ -> sum(
            abs2, transmission_loss(field_solution(θ), FIELD_RANGES, FIELD_ZS, FIELD_ZR; nmodes=3, mode=:incoherent)
        ),
    ),
]

@testset "Reverse-mode AD through the pressure field" begin
    @testset "the environment supports the field the targets ask for" begin
        # If this fails every comparison below is measuring a degenerate case rather than a field.
        sol = field_solution(θ7)
        @test length(sol.kr) >= 3
        @test sol.env.bottom_bc isa AcousticHalfspace     # …so FIELD_ZR's 120 m is admissible
        @test all(isfinite, transmission_loss(sol, FIELD_RANGES, FIELD_ZS, FIELD_ZR; nmodes=3))
    end

    @testset "$name — Zygote vs ForwardDiff" for (name, target) in FIELD_TARGETS
        g_forward = ForwardDiff.gradient(target, θ7)
        g_zygote = Zygote.gradient(target, θ7)[1]

        @test length(g_zygote) == length(θ7)
        @test all(isfinite, g_zygote)
        # 7.5's acceptance criterion is 1e-6; measured 1e-12 and below on all four targets.
        @test relerr_norm(g_zygote, g_forward) < 1e-6

        # Both routes into the field have to be live. A field gradient that reached the mode shapes
        # but not the half-space would leave the last entry at exactly zero, and vice versa.
        @test all(!iszero, g_zygote[1:10])
        @test !iszero(g_zygote[11])
    end

    @testset "against central differences, which use no AD at all" begin
        # Only the profile entries: `cb` sets `Δz = cb / (20f)` and so the mesh point count, which
        # makes the primal a step function of it — see the comment on `ENVIRONMENTS`. Differencing
        # the ten sound speeds is a real independent check; differencing `cb` would measure the mesh.
        target = FIELD_TARGETS[3][2]     # transmission loss, the criterion's own quantity
        g_forward = ForwardDiff.gradient(target, θ7)
        @test relerr_norm(central_gradient(target, θ7, 1:10), g_forward[1:10]) < 1e-5
    end

    @testset "Mooncake reaches the same field gradient" begin
        # Zygote and Mooncake read the same rules by different routes (see the backend table above),
        # and the field is new ground for that route: it is the first place the cotangent flows
        # through comprehensions and a generator sum before it reaches a rule.
        target = FIELD_TARGETS[3][2]
        g_zygote = Zygote.gradient(target, θ7)[1]
        g_mooncake = DifferentiationInterface.gradient(target, AutoMooncake(), θ7)
        @test relerr_norm(g_mooncake, g_zygote) < 1e-11
    end

    @testset "a receiver below the bottom differentiates through the half-space tail" begin
        # `mode_amplitude`'s `Φ[end, m] * exp(-γ(z - D))` branch is the only place `cb` reaches a
        # *receiver*, and it is reached only for a depth past the last mesh point. Comparing a
        # deep-receiver target against a water-column one isolates it: the deep one must be strictly
        # more sensitive to `cb`, or the exponential tail is not on the tape.
        shallow = θ -> sum(abs2, transmission_loss(field_solution(θ), FIELD_RANGES, FIELD_ZS, 50.0; nmodes=3))
        deep = θ -> sum(abs2, transmission_loss(field_solution(θ), FIELD_RANGES, FIELD_ZS, 150.0; nmodes=3))

        g_shallow = Zygote.gradient(shallow, θ7)[1]
        g_deep = Zygote.gradient(deep, θ7)[1]
        @test relerr_norm(g_deep, ForwardDiff.gradient(deep, θ7)) < 1e-6
        @test abs(g_deep[11]) > abs(g_shallow[11])
    end

    @testset "a lossy field, where the pressure is genuinely complex" begin
        # Attenuation makes `kr` complex, so `exp(-i kr r)` acquires a decay and the field's own
        # arithmetic goes complex on a path the lossless targets never exercise. The convention is
        # the one M5 documents — real parameters in, real loss out — and `transmission_loss` is real
        # by construction, so it is admissible without further comment.
        θ_lossy7 = [θ7[1:10]; 0.1; θ7[11]]
        target = lossy_field_loss

        @test is_lossy(lossy_field_env(θ_lossy7))
        g_forward = ForwardDiff.gradient(target, θ_lossy7)
        g_zygote = Zygote.gradient(target, θ_lossy7)[1]
        @test all(isfinite, g_zygote)
        @test relerr_norm(g_zygote, g_forward) < 1e-6
        # More water-column loss means more TL: the attenuation derivative is positive and real.
        @test g_forward[11] > 0
    end

    @testset "broadband synthesis is differentiable too" begin
        # `broadband_field` used to fill an `Array{ComplexF64,3}(undef, …)` plane by plane, and the
        # three-argument `synthesize_pulse` a `Array{Float64,3}` trace by trace. Both are mutation,
        # which Zygote refuses outright, and both pinned an element type that a `Dual` cannot fit
        # into — so the header's claim that the hand-written transform was "differentiable with no
        # extra rules for task 7.5" was true of the transform and false of the function in front of
        # it. Task 7.5 restated both as comprehensions; this is the test that keeps them that way.
        bfreqs = collect(80.0:20.0:140.0)
        times = collect(0.0:0.05:0.4)
        pulse(θ) = sum(
            abs2,
            synthesize_pulse(
                broadband_field(field_env(θ), bfreqs, 5_000.0, FIELD_ZS, 50.0; nmodes=3, abstol=1e-10, reltol=1e-10),
                bfreqs,
                times,
            ),
        )

        @test pulse(θ7) > 0
        g_forward = ForwardDiff.gradient(pulse, θ7)
        g_zygote = Zygote.gradient(pulse, θ7)[1]
        @test all(isfinite, g_zygote)
        @test relerr_norm(g_zygote, g_forward) < 1e-6
    end
end
