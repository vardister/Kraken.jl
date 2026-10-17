# Unit tests for the solver internals behind `kraken_jl`'s speed: the safeguarded Newton root solver
# and its starting guesses, the per-mesh energy weights, the tridiagonal factorization the inverse
# iteration uses, `sturm_and_slope`, and `primal_value`. The end-to-end tests elsewhere (Fortran
# cross-validation, AD agreement) would catch a regression in any of these, but only as a mismatch
# several layers away. These name the piece that broke.

using Test
using Kraken
using LinearAlgebra
using Random
using ForwardDiff
using ChainRulesTestUtils

@testset "Solver internals" begin
    pekeris = UnderwaterEnv(pekeris_env()...)
    onelayer = UnderwaterEnv(one_layer_env()...)
    twolayer = UnderwaterEnv(two_layer_slope_env()...)
    munk = UnderwaterEnv(munk_env()...)

    # Every top and bottom boundary the solver accepts, on the same column, so the boundary-dependent
    # branches (`has_top_node`, `close_medium`, the rigid-rigid starting vector) are each reached.
    ssp, layers, sspHS = pekeris_env()
    boundary_envs = [
        UnderwaterEnv(ssp, layers, sspHS; top_bc=top, bottom_bc=bottom) for top in (PressureRelease(), RigidBoundary())
        for bottom in (PressureRelease(), RigidBoundary(), AcousticHalfspace())
    ]

    setup(env, f; factor=1) = (p=AcousticProblemProperties(env, f; factor=factor); (p, AcousticProblemCache(env, p)))
    relerr(a, b) = maximum(abs.(a .- b) ./ abs.(b))

    # `Brent()` reaches Kraken.jl through NonlinearSolve, which the test environment does not
    # declare, so it is named through the package.
    @testset "SafeguardedNewton finds Brent's roots" begin
        for (env, f) in ((pekeris, 100.0), (onelayer, 100.0), (twolayer, 100.0), (munk, 50.0))
            p, c = setup(env, f)
            newton = find_kr(env, p, c)
            brent = find_kr(env, p, c; method=Kraken.Brent(), abstol=1e-12, reltol=1e-12)
            @test length(newton) == length(brent)
            @test relerr(newton, brent) < 1e-12
        end
        for env in boundary_envs
            p, c = setup(env, 100.0)
            @test relerr(find_kr(env, p, c), find_kr(env, p, c; method=Kraken.Brent(), abstol=1e-12, reltol=1e-12)) <
                1e-12
        end
    end

    @testset "starting guesses change no root" begin
        for (env, f) in ((pekeris, 100.0), (twolayer, 100.0), (munk, 50.0))
            p, c = setup(env, f)
            p2, c2 = setup(env, f; factor=2)
            coarse = find_kr(env, p, c)
            @test relerr(find_kr(env, p2, c2; guesses=coarse), find_kr(env, p2, c2)) < 1e-12
            # Fewer guesses than modes: the rest are solved from their brackets as usual.
            @test relerr(find_kr(env, p2, c2; guesses=coarse[1:2]), find_kr(env, p2, c2)) < 1e-12

            # A guess far outside the bracket is clamped into it and still finds that bracket's mode.
            spans = bisection(env, p2, c2)
            for m in (1, size(spans, 1))
                plain = solve_for_kr(spans[m, :], env, p2, c2)
                @test abs(solve_for_kr(spans[m, :], env, p2, c2; guess=10 * spans[m, 2]) - plain) / plain < 1e-12
                @test abs(solve_for_kr(spans[m, :], env, p2, c2; guess=0.0) - plain) / plain < 1e-12
            end
        end
    end

    @testset "mode_energy_weights reproduce the medium-by-medium trapezoid" begin
        # The integral as `mode_energy` computed it before the weights: each medium's closed mesh,
        # closed mode and closed density, integrated separately.
        function trapezoid_energy(ψ, env, p)
            zz = Kraken.medium_mesh(env, p)
            vv = Kraken.medium_mode(ψ, env, p)
            ρρ = Kraken.medium_property(Kraken.density, env.ρ, env, p)
            return sum(i -> Kraken.integral_trapz(abs2.(vv[i]) ./ ρρ[i], zz[i]), eachindex(zz))
        end
        for env in [boundary_envs; onelayer; twolayer]
            p = AcousticProblemProperties(env, 100.0)
            ψ = sin.(range(0.0, 7.0; length=sum(p.Nz_vec))) .+ 0.1
            @test Kraken.mode_energy(ψ, env, p) ≈ trapezoid_energy(ψ, env, p) rtol = 1e-13
            @test length(Kraken.mode_energy_weights(env, p)) == sum(p.Nz_vec)
        end
    end

    @testset "mode_energy: rule and allocation" begin
        rng = Xoshiro(20261016)
        ψ = randn(rng, 40)
        w = rand(rng, 40)
        test_rrule(Kraken.mode_energy, ψ, w)
        @test ForwardDiff.gradient(x -> Kraken.mode_energy(x, w), ψ) ≈ 2 .* w .* ψ
        allocated(f, args...) = (f(args...); @allocated f(args...))
        @test allocated(Kraken.mode_energy, ψ, w) == 0
    end

    @testset "shifted_lu! and shifted_ldiv! solve like lu" begin
        rng = Xoshiro(7)
        n = 60
        ws = Kraken.eigen_workspace(Float64, n)
        for A in (
            Tridiagonal(randn(rng, n - 1), randn(rng, n) .+ 4, randn(rng, n - 1)),   # pivots rarely
            Tridiagonal(ones(n - 1), fill(1e-8, n), ones(n - 1)),                    # pivots on every row
        )
            A0 = copy(A)
            b = randn(rng, n)
            x = Kraken.shifted_ldiv!(copy(b), Kraken.shifted_lu!(ws, A))
            @test x ≈ lu(A) \ b rtol = 1e-12
            @test A == A0
        end
        @test_throws SingularException Kraken.shifted_lu!(ws, Tridiagonal(zeros(n - 1), zeros(n), zeros(n - 1)))
    end

    @testset "mode_eigenvector allocates little beyond its result" begin
        p, c = setup(munk, 50.0)
        kr = first(find_kr(munk, p, c))
        Kraken.mode_eigenvector(kr, munk, p, c)
        # The returned vector is N Float64s. Before the task-local workspace this was 253 KiB on this
        # 3125-point mesh, about ten vectors; it measured 28.2 KiB after.
        @test (@allocated Kraken.mode_eigenvector(kr, munk, p, c)) < 2 * sizeof(Float64) * sum(p.Nz_vec)
    end

    @testset "sturm_and_slope is det_sturm and its derivative" begin
        for (env, f) in ((pekeris, 100.0), (twolayer, 100.0), (munk, 50.0), (boundary_envs[4], 100.0))
            p, c = setup(env, f)
            spans = bisection(env, p, c)
            for kr in range(spans[end, 1], spans[1, 2]; length=25)
                D, dD = Kraken.sturm_and_slope(kr, env, p, c)
                @test D == first(det_sturm(kr, env, p, c))
                @test dD ≈ ForwardDiff.derivative(k -> first(det_sturm(k, env, p, c)), kr) rtol = 1e-10
            end
        end
    end

    @testset "primal_value strips Duals" begin
        d = ForwardDiff.Dual(1.5, 2.0)
        @test Kraken.primal_value(d) === 1.5
        @test Kraken.primal_value(ForwardDiff.Dual(d, d)) === 1.5
        @test Kraken.primal_value(1.5) === 1.5
        # The reason it exists: ForwardDiff orders two Duals with equal values by their partials, so
        # a branch on `<` can go one way for a Dual and the other for its plain value.
        @test ForwardDiff.Dual(1.0, -1.0) < ForwardDiff.Dual(1.0, 1.0)
        @test !(1.0 < 1.0)
    end
end
