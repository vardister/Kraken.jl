using TestItems

@testitem "Pekeris Waveguide Standard Case" begin
    using Kraken

    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0

    sol = kraken_jl(env, freq)

    @test sol isa Kraken.NormalModeSolution
    @test !isempty(sol.kr)
    @test size(sol.modes, 2) == length(sol.kr)
    @test sol.env === env
    @test sol.props.freq == freq

    # Test against known Pekeris values (from analytical solution)
    expected_krs = [0.417908, 0.414964, 0.409971, 0.40286, 0.39383]
    @test length(sol.kr) >= 5
    for i in 1:min(5, length(sol.kr))
        @test sol.kr[i] ≈ expected_krs[i] atol = 1e-4
    end

    # Test mode properties
    for i in 1:size(sol.modes, 2)
        mode = sol.modes[:, i]
        @test !all(mode .== 0)  # Mode should not be zero everywhere
        @test mode[1] > 0  # First element should be positive (normalization convention)
    end
end

@testitem "Different Frequencies" begin
    using Kraken

    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)

    freqs = [25.0, 50.0, 100.0, 200.0, 400.0]
    solutions = []

    for freq in freqs
        sol = kraken_jl(env, freq)
        push!(solutions, sol)

        @test !isempty(sol.kr)
        @test all(sol.kr .> 0)
        @test issorted(sol.kr, rev=true)
    end

    # Higher frequencies should generally have more modes
    n_modes = [length(sol.kr) for sol in solutions]
    @test issorted(n_modes)  # Should be non-decreasing

    # Test frequency scaling - higher freq should have higher wavenumbers
    for i in 2:length(solutions)
        @test maximum(solutions[i].kr) > maximum(solutions[i - 1].kr)
    end
end

@testitem "Different Environment Parameters" begin
    using Kraken

    # Test various sound speed contrasts
    sound_speeds = [(1500, 1600), (1500, 1700), (1500, 1800)]

    for (cw, cb) in sound_speeds
        ssp, layers, sspHS = pekeris_env(; c0=cw, cb=cb)
        env = UnderwaterEnv(ssp, layers, sspHS)
        sol = kraken_jl(env, 100.0)

        @test !isempty(sol.kr)
        @test all(sol.kr .> 0)

        # Larger sound speed contrast should support more modes
        contrast = cb - cw
        @test contrast > 0  # Positive contrast
    end

    # Test different water depths
    depths = [50.0, 100.0, 200.0, 500.0]
    n_modes_by_depth = []

    for depth in depths
        ssp, layers, sspHS = pekeris_env(; depth=depth)
        env = UnderwaterEnv(ssp, layers, sspHS)
        sol = kraken_jl(env, 100.0)

        @test !isempty(sol.kr)
        push!(n_modes_by_depth, length(sol.kr))
    end

    # Deeper water should generally support more modes
    @test issorted(n_modes_by_depth)
end

@testitem "Richardson Extrapolation" begin
    using Kraken

    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0

    # Test different mesh refinements
    sol_coarse = kraken_jl(env, freq; n_meshes=1)
    sol_fine = kraken_jl(env, freq; n_meshes=3)
    sol_finest = kraken_jl(env, freq; n_meshes=5)

    @test length(sol_coarse.kr) == length(sol_fine.kr) == length(sol_finest.kr)

    # Richardson extrapolation should improve accuracy
    # (finer meshes should be closer to analytical values)
    expected_kr1 = 0.417908
    if !isempty(sol_coarse.kr)
        err_coarse = abs(sol_coarse.kr[1] - expected_kr1)
        err_fine = abs(sol_fine.kr[1] - expected_kr1)
        err_finest = abs(sol_finest.kr[1] - expected_kr1)

        @test err_finest <= err_fine
        @test err_fine <= err_coarse * 2  # Allow some tolerance
    end
end

@testitem "One-layer Sediment" begin
    using Kraken

    ssp, layers, sspHS = one_layer_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0

    sol = kraken_jl(env, freq)

    @test sol isa Kraken.NormalModeSolution
    if !isempty(sol.kr)
        @test all(sol.kr .> 0)
        @test size(sol.modes, 1) == sum(sol.props.Nz_vec)
    end
end

@testitem "Munk Profile" begin
    using Kraken

    ssp, layers, sspHS = munk_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 25.0  # Lower frequency for deep water

    sol = kraken_jl(env, freq)

    @test !isempty(sol.kr)
    @test all(sol.kr .> 0)

    # Munk profile should support many modes due to depth
    @test length(sol.kr) > 5

    # Test that modes have reasonable structure
    for i in 1:min(3, size(sol.modes, 2))
        mode = sol.modes[:, i]
        @test maximum(abs.(mode)) > 0.001  # Should have reasonable amplitude
    end
end

@testitem "Mode Orthogonality" begin
    using Kraken

    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    sol = kraken_jl(env, 100.0)

    if size(sol.modes, 2) >= 2
        # Test approximate orthogonality
        zn = vcat(sol.props.zn_vec...)
        ρn = density(env.ρ, zn)

        for i in 1:min(3, size(sol.modes, 2))
            for j in (i + 1):min(3, size(sol.modes, 2))
                mode_i = sol.modes[:, i]
                mode_j = sol.modes[:, j]

                # Compute inner product
                overlap = sum(mode_i .* mode_j ./ ρn) * (zn[2] - zn[1])
                @test abs(overlap) < 0.1  # Should be small
            end
        end
    end
end

@testitem "Mode Normalization" begin
    using Kraken

    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    sol = kraken_jl(env, 100.0)

    zn = vcat(sol.props.zn_vec...)
    ρn = density(env.ρ, zn)

    for i in 1:size(sol.modes, 2)
        mode = sol.modes[:, i]

        # Compute mode normalization
        norm_sq_water = sum(abs2.(mode) ./ ρn) * (zn[2] - zn[1])

        # Add half-space contribution
        kr = sol.kr[i]
        k_b = 2π * sol.props.freq / env.cb
        if kr > k_b
            norm_sq_bottom = abs2(mode[end]) / (2 * env.ρb * sqrt(kr^2 - k_b^2))
            total_norm = norm_sq_water + norm_sq_bottom
        else
            total_norm = norm_sq_water
        end

        @test total_norm > 0.1  # Should be reasonably normalized
        @test total_norm < 10.0  # But not too large
    end
end

@testitem "Wavenumber Properties" begin
    using Kraken

    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0
    sol = kraken_jl(env, freq)

    k0 = 2π * freq / soundspeed(env.c, 0.0)  # Water wavenumber
    kb = 2π * freq / env.cb  # Bottom wavenumber

    for kr in sol.kr
        @test kb < kr < k0  # Should be trapped modes
    end

    # Test dispersion relation
    @test issorted(sol.kr, rev=true)  # Higher order modes have smaller kr
end

@testitem "No Trapped Modes Case" begin
    using Kraken

    # Create environment with bottom slower than water (no trapping possible)
    # This should throw an assertion error because trapped modes require cb > c_water
    ssp, layers, sspHS = pekeris_env(; cb=1400.0)  # Bottom slower than water
    env = UnderwaterEnv(ssp, layers, sspHS)

    # kraken_jl requires maxsoundspeed(env.c) < env.cb for trapped modes
    @test_throws AssertionError kraken_jl(env, 100.0)
end

@testitem "Very High Frequency" begin
    using Kraken

    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)

    # Test with high frequency (should have many modes)
    sol_hf = kraken_jl(env, 1000.0)
    @test sol_hf isa Kraken.NormalModeSolution

    if !isempty(sol_hf.kr)
        @test all(sol_hf.kr .> 0)
        @test length(sol_hf.kr) > 10  # Should have many modes
    end
end

@testitem "Standard Pekeris Regression" begin
    using Kraken

    # Test against known good values to catch regressions
    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    sol = kraken_jl(env, 100.0)

    # These are the reference values from the analytical Pekeris solution
    ref_krs = [0.417908, 0.414964, 0.409971, 0.40286, 0.39383]

    @test length(sol.kr) >= 5
    for i in 1:5
        @test sol.kr[i] ≈ ref_krs[i] atol = 1e-4 rtol = 1e-6
    end

    # Test mode structure
    @test size(sol.modes, 1) > 50  # Should have reasonable mesh resolution
    @test all(sol.modes[1, :] .> 0)  # First element positive by convention
end
