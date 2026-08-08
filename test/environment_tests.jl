using Test
using TestItems
using Kraken

@testitem "SampledSSP1D Construction and Evaluation" begin
    using Kraken

    depths = [0.0, 50.0, 100.0]
    speeds = [1500.0, 1480.0, 1490.0]

    ssp = SampledSSP(depths, speeds)

    # Test interpolation at known points
    @test soundspeed(ssp, 0.0) ≈ 1500.0
    @test soundspeed(ssp, 50.0) ≈ 1480.0
    @test soundspeed(ssp, 100.0) ≈ 1490.0

    # Test interpolation at intermediate points
    @test soundspeed(ssp, 25.0) ≈ 1490.0  # Linear interpolation

    # Test maximum sound speed
    @test maxsoundspeed(ssp) == 1500.0
end

@testitem "SampledDensity1D Construction and Evaluation" begin
    using Kraken

    depths = [0.0, 50.0, 100.0]
    densities = [1000.0, 1020.0, 1030.0]

    ρ = SampledDensity(depths, densities)

    # Test interpolation at known points
    @test density(ρ, 0.0) ≈ 1000.0
    @test density(ρ, 50.0) ≈ 1020.0
    @test density(ρ, 100.0) ≈ 1030.0

    # Test interpolation at intermediate points
    @test density(ρ, 25.0) ≈ 1010.0  # Linear interpolation
end

@testitem "UnderwaterEnv Construction" begin
    using Kraken

    # Test with Pekeris environment
    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)

    @test env.depth == 100.0
    @test env.cb ≈ 1600.0
    @test env.ρb ≈ 1500.0
    @test soundspeed(env.c, 0.0) ≈ 1500.0
    @test density(env.ρ, 0.0) ≈ 1000.0

    # Test different parameters
    ssp2, layers2, sspHS2 = pekeris_env(c0=1520.0, cb=1650.0, depth=200.0)
    env2 = UnderwaterEnv(ssp2, layers2, sspHS2)

    @test env2.depth == 200.0
    @test env2.cb ≈ 1650.0
    @test soundspeed(env2.c, 0.0) ≈ 1520.0
end

@testitem "Standard Environment Functions" begin
    using Kraken
    using Test

    @testset "Pekeris Environment" begin
        ssp, layers, sspHS = pekeris_env()

        # Check matrix dimensions and basic structure
        @test size(ssp, 2) == 6  # depth, c, α, ρ, etc.
        @test size(layers, 2) == 3
        @test size(sspHS, 2) == 6
        @test sspHS[2, 2] ≈ 1600.0  # Bottom sound speed
        @test sspHS[2, 4] ≈ 1500.0  # Bottom density

        # Test custom parameters
        ssp_custom, _, sspHS_custom = pekeris_env(c0=1520.0, cb=1700.0, ρ0=1050.0, ρb=2000.0)
        @test ssp_custom[1, 2] ≈ 1520.0
        @test sspHS_custom[2, 2] ≈ 1700.0
        @test sspHS_custom[2, 4] ≈ 2000.0
    end

    @testset "One Layer Environment" begin
        ssp, layers, sspHS = one_layer_env()

        @test size(ssp, 1) == 4  # Water surface, water bottom, sediment top, sediment bottom
        @test size(layers, 1) == 2  # Two layers
        @test ssp[1, 2] ≈ 1500.0  # Water sound speed
        @test ssp[3, 2] ≈ 1550.0  # Sediment sound speed
    end

    @testset "Munk Profile" begin
        ssp, layers, sspHS = munk_env()

        @test size(ssp, 1) > 10  # Should have many depth points
        @test minimum(ssp[:, 2]) < 1490.0  # Should have sound speed minimum
        @test maximum(ssp[:, 2]) > 1520.0  # Should have sound speed maximum

        # Test that Munk profile has the characteristic minimum
        depths = -ssp[:, 1]
        speeds = ssp[:, 2]
        min_idx = argmin(speeds)
        min_depth = depths[min_idx]
        @test min_depth > 1000.0 && min_depth < 1500.0  # Munk minimum around 1300m
    end
end

@testitem "Error Handling" begin
    using Kraken

    # Test with invalid inputs
    @test_throws BoundsError SampledSSP(Float64[], Float64[])  # Empty arrays

    # Test environment with mismatched dimensions
    ssp_bad = [0.0 1500.0 0.0 1000.0 0.0]  # Wrong number of columns
    layers_good, sspHS_good = pekeris_env()[2:3]
    @test_throws Exception UnderwaterEnv(ssp_bad, layers_good, sspHS_good)
end
