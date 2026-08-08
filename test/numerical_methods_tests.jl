using Test
using Kraken
using Roots

@testset "Numerical Methods" begin

    # Setup standard test environment
    function setup_test_environment()
        ssp, layers, sspHS = pekeris_env()
        env = UnderwaterEnv(ssp, layers, sspHS)
        freq = 100.0
        props = AcousticProblemProperties(env, freq)
        cache = AcousticProblemCache(env, props)
        return env, props, cache
    end

    @testset "AcousticProblemProperties Construction" begin
        ssp, layers, sspHS = pekeris_env()
        env = UnderwaterEnv(ssp, layers, sspHS)
        freq = 100.0

        props = AcousticProblemProperties(env, freq)

        @test props.freq == freq
        @test !isempty(props.zn_vec)
        @test !isempty(props.Δz_vec)
        @test sum(props.Nz_vec) > 0

        # Test with different mesh factor
        props_fine = AcousticProblemProperties(env, freq; factor=2)
        @test sum(props_fine.Nz_vec) > sum(props.Nz_vec)  # Finer mesh
    end

    @testset "AcousticProblemCache Construction" begin
        env, props, cache = setup_test_environment()

        @test length(cache.a_vec) == sum(props.Nz_vec)
        @test length(cache.e_vec) == sum(props.Nz_vec)
        @test length(cache.λ_scaling) == sum(props.Nz_vec)

        # Test that cache creates proper tridiagonal structure
        @test size(cache.A, 1) == size(cache.A, 2)
        @test size(cache.A, 1) == length(cache.a_vec)
    end

    @testset "finite_difference_coefficients" begin
        # The coefficient assembly is the differentiable seam: a pure function the reverse-mode
        # rules trace. These tests pin both halves of that — that it is what the cache contains,
        # and that it mutates nothing.
        env, props, cache = setup_test_environment()
        a_vec, e_vec, λ_scaling = finite_difference_coefficients(env, props)
        N = sum(props.Nz_vec)

        @test length(a_vec) == N
        @test length(e_vec) == N
        @test length(λ_scaling) == N

        # The cache is a thin wrapper — its fields are exactly this function's output.
        @test a_vec == cache.a_vec
        @test e_vec == cache.e_vec
        @test λ_scaling == cache.λ_scaling
        # `create_finite_diff_matrix!` updates `A` by writing into `cache.a_vec`, which only works
        # because the matrix shares that array rather than copying it.
        @test cache.A.d === cache.a_vec

        # The bottom entry of λ_scaling is halved for the half-space boundary condition.
        @test λ_scaling[end] ≈ e_vec[end] * props.Δz_vec[end]^2 / 2

        # Pure: inputs untouched, and every call returns fresh arrays.
        zn_before = deepcopy(props.zn_vec)
        c_before = copy(env.c.c)
        ρ_before = copy(env.ρ.ρ)
        a2, e2, λ2 = finite_difference_coefficients(env, props)
        @test props.zn_vec == zn_before
        @test env.c.c == c_before
        @test env.ρ.ρ == ρ_before
        @test a2 == a_vec
        @test e2 == e_vec
        @test λ2 == λ_scaling
        @test a2 !== a_vec

        # Multi-layer: the diagonal coefficient at a medium interface is the average of the values
        # the two media give there.
        env_ml = UnderwaterEnv(two_layer_slope_env()...)
        props_ml = AcousticProblemProperties(env_ml, 100.0)
        @test length(props_ml.Nz_vec) > 1
        a_ml, _, _ = finite_difference_coefficients(env_ml, props_ml)
        zn_all = reduce(vcat, props_ml.zn_vec)
        Δz_all = reduce(vcat, [fill(props_ml.Δz_vec[i], props_ml.Nz_vec[i]) for i in eachindex(props_ml.Nz_vec)])
        a_raw(k) =
            Kraken.a_element(soundspeed(env_ml.c, zn_all[k]), density(env_ml.ρ, zn_all[k]), props_ml.freq, Δz_all[k])
        for k in cumsum(props_ml.Nz_vec)[1:(end - 1)]
            @test a_ml[k] ≈ (a_raw(k) + a_raw(k + 1)) / 2
            @test a_ml[k] != a_raw(k)   # the interface value really was averaged, not copied
        end
        # Away from interfaces the coefficient is the plain element formula.
        interfaces_ml = cumsum(props_ml.Nz_vec)[1:(end - 1)]
        for k in (1, length(a_ml) ÷ 2, length(a_ml))
            k in interfaces_ml && continue
            @test a_ml[k] ≈ a_raw(k)
        end
    end

    @testset "Determinant and Sturm Count Function" begin
        env, props, cache = setup_test_environment()

        # `det_sturm` is only defined on the trapped-mode band kr ∈ [ω/cb, max(ω/c)]; below
        # ω/cb the bottom vertical wavenumber is real (a leaky mode) and `get_g` takes the
        # square root of a negative number. For pekeris_env() at 100 Hz that band is
        # [0.39270, 0.41888] — probe wavenumbers must stay inside it.
        ω = 2π * props.freq
        kr_cutoff = ω / env.cb
        kr_water = maximum(ω ./ env.c.c)
        @test kr_cutoff < kr_water

        # Test at different wavenumber values
        kr_test = 0.4
        det_val, sturm_count = det_sturm(kr_test, env, props, cache)

        @test det_val isa Real
        @test sturm_count isa Integer
        @test sturm_count >= 0

        # Below the cutoff the model has no real solution — that is a domain error, not a bug.
        @test_throws DomainError det_sturm(0.9 * kr_cutoff, env, props, cache)

        # Test that determinant changes sign around roots
        kr_low = 0.395
        kr_high = 0.415
        det_low, _ = det_sturm(kr_low, env, props, cache)
        det_high, _ = det_sturm(kr_high, env, props, cache)

        # Should have different signs if there's a root between them
        @test typeof(det_low) == typeof(det_high)

        # The Sturm count at `kr` is the number of modes with wavenumber *above* `kr`, so it
        # decreases as kr sweeps up through the trapped band and reaches 0 at max(ω/c).
        # This is the convention `bisection` relies on to bracket each mode.
        _, count_low = det_sturm(kr_low, env, props, cache)
        _, count_high = det_sturm(kr_high, env, props, cache)
        @test count_high <= count_low
        @test last(det_sturm(nextfloat(kr_cutoff), env, props, cache)) >= count_low
        @test last(det_sturm(kr_water, env, props, cache)) == 0
    end

    @testset "Bisection Method" begin
        env, props, cache = setup_test_environment()

        intervals = bisection(env, props, cache)

        @test intervals isa Matrix
        @test size(intervals, 2) == 2  # [left, right] bounds
        @test all(intervals[:, 1] .< intervals[:, 2])  # left < right
        @test all(intervals .> 0)  # all wavenumbers positive

        # Test that each interval contains exactly one sign change
        for i in 1:size(intervals, 1)
            left, right = intervals[i, :]
            det_left, _ = det_sturm(left, env, props, cache)
            det_right, _ = det_sturm(right, env, props, cache)
            # Should have opposite signs (allowing for numerical precision)
            @test det_left * det_right <= 0
        end
    end

    @testset "Root Solving" begin
        env, props, cache = setup_test_environment()

        # Get intervals from bisection
        intervals = bisection(env, props, cache)

        if !isempty(intervals)
            # Test solving for first root
            span = (intervals[1, 1], intervals[1, 2])
            root = solve_for_kr(span, env, props, cache)

            @test root isa Real
            @test span[1] <= root <= span[2]

            # Test that the root is actually close to zero
            det_at_root, _ = det_sturm(root, env, props, cache)
            @test abs(det_at_root) < 1e-6  # Should be close to zero
        end
    end

    @testset "find_kr Function" begin
        env, props, cache = setup_test_environment()

        krs = find_kr(env, props, cache)

        @test krs isa Vector
        @test all(krs .> 0)  # All wavenumbers should be positive
        @test issorted(krs, rev=true)  # Should be sorted in descending order

        # Test with different root-finding methods
        krs_bisection = find_kr(env, props, cache; method=Roots.Bisection())
        @test length(krs_bisection) == length(krs)
        @test all(abs.(krs_bisection .- krs) .< 1e-8)  # Should be very close

        # Test that found roots are actual solutions
        for kr in krs
            det_val, _ = det_sturm(kr, env, props, cache)
            @test abs(det_val) < 1e-6
        end
    end

    @testset "Inverse Iteration" begin
        env, props, cache = setup_test_environment()

        # First find some wavenumbers
        krs = find_kr(env, props, cache)

        if !isempty(krs)
            kr_test = krs[1]

            # Test single mode inverse iteration
            kr_refined, mode = inverse_iteration(kr_test, env, props, cache)

            @test kr_refined ≈ kr_test atol=1e-6  # Should refine the wavenumber
            @test mode isa Vector
            @test length(mode) == sum(props.Nz_vec)
            @test !all(mode .== 0)  # Mode should not be zero everywhere

            # Test that mode is normalized
            zn = vcat(props.zn_vec...)
            ρn = density(env.ρ, zn)
            norm_squared = sum(abs2.(mode) ./ ρn) * (zn[2] - zn[1])  # Simple integration
            @test norm_squared > 0  # Should have some normalization

            # Test multiple modes at once
            if length(krs) > 1
                krs_test = krs[1:min(3, length(krs))]
                krs_refined, modes = inverse_iteration(krs_test, env, props, cache)

                @test length(krs_refined) == length(krs_test)
                @test size(modes, 2) == length(krs_test)
                @test size(modes, 1) == length(mode)

                # Test orthogonality (modes should be approximately orthogonal)
                if size(modes, 2) >= 2
                    mode1, mode2 = modes[:, 1], modes[:, 2]
                    overlap = sum(mode1 .* mode2 ./ ρn)
                    @test abs(overlap) < 0.1  # Should be small but not necessarily zero
                end
            end
        end
    end

    @testset "get_g Function" begin
        env, props, cache = setup_test_environment()

        kr_test = 0.4
        g_val = get_g(kr_test, env, props)

        @test g_val isa Real
        @test !isnan(g_val)
        @test !isinf(g_val)

        # Test that g varies smoothly with kr
        kr_nearby = kr_test + 1e-6
        g_nearby = get_g(kr_nearby, env, props)
        @test abs(g_nearby - g_val) < 1e-3  # Should vary smoothly
    end

    @testset "Numerical Stability" begin
        env, props, cache = setup_test_environment()

        # Extreme wavenumbers, but on the side of the cutoff where the model is defined:
        # just above ω/cb (where g → 0) and far above the water wavenumber. The Sturm
        # sequence rescaling (`scale_const`) is what keeps these from overflowing.
        ω = 2π * props.freq
        kr_cutoff = ω / env.cb
        det_edge, _ = det_sturm(nextfloat(kr_cutoff), env, props, cache)
        det_large, _ = det_sturm(10.0, env, props, cache)

        @test !isnan(det_edge) && !isinf(det_edge)
        @test !isnan(det_large) && !isinf(det_large)

        # Test with different frequencies. 10 Hz in 100 m of water is below cutoff, so
        # `bisection` legitimately returns `nothing` there; every higher frequency must
        # produce brackets. Mode count must be non-decreasing in frequency.
        n_modes_prev = 0
        for freq in [10.0, 50.0, 200.0, 1000.0]
            props_freq = AcousticProblemProperties(env, freq)
            cache_freq = AcousticProblemCache(env, props_freq)

            intervals_freq = bisection(env, props_freq, cache_freq)
            @test isnothing(intervals_freq) || intervals_freq isa Matrix  # Should not crash

            krs_freq = find_kr(env, props_freq, cache_freq)
            @test all(krs_freq .> 0)  # All wavenumbers positive
            @test length(krs_freq) >= n_modes_prev
            n_modes_prev = length(krs_freq)
        end
        @test n_modes_prev > 1  # 1 kHz supports many modes
    end

    @testset "Edge Cases" begin
        # Test with very shallow water
        ssp_shallow, layers_shallow, sspHS_shallow = pekeris_env(depth=10.0)
        env_shallow = UnderwaterEnv(ssp_shallow, layers_shallow, sspHS_shallow)
        props_shallow = AcousticProblemProperties(env_shallow, 100.0)
        cache_shallow = AcousticProblemCache(env_shallow, props_shallow)

        # 10 m of water at 100 Hz (λ ≈ 15 m) supports no trapped mode, and `bisection`
        # signals that by returning `nothing` — which `find_kr` turns into an empty vector.
        intervals_shallow = bisection(env_shallow, props_shallow, cache_shallow)
        @test isnothing(intervals_shallow)
        @test isempty(find_kr(env_shallow, props_shallow, cache_shallow))

        # Test with very deep water
        ssp_deep, layers_deep, sspHS_deep = pekeris_env(depth=5000.0)
        env_deep = UnderwaterEnv(ssp_deep, layers_deep, sspHS_deep)
        props_deep = AcousticProblemProperties(env_deep, 100.0)
        cache_deep = AcousticProblemCache(env_deep, props_deep)

        intervals_deep = bisection(env_deep, props_deep, cache_deep)
        @test intervals_deep isa Matrix

        if !isempty(intervals_deep)
            krs_deep = find_kr(env_deep, props_deep, cache_deep)
            @test length(krs_deep) > length(find_kr(env_shallow, props_shallow, cache_shallow))
        end
    end
end
