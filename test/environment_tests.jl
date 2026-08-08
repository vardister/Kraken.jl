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

        # Depths in column 1 are already positive (z increases downwards from the
        # surface); do NOT negate them.
        depths = ssp[:, 1]
        speeds = ssp[:, 2]
        min_idx = argmin(speeds)

        # The Munk profile is *defined* by its sound-channel axis at z = 1300 m, where
        # ẑ = 2(z - 1300)/1300 vanishes and the perturbation term ε(ẑ - 1 + exp(-ẑ))
        # is therefore exactly zero. So the minimum is exactly 1500.0 m/s at exactly
        # 1300 m -- not "around" either value. If this fails, the profile coefficients
        # have moved; do not repair it by loosening the bounds.
        @test speeds[min_idx] ≈ 1500.0 atol = 1e-10
        @test depths[min_idx] ≈ 1300.0

        # Sound speed increases monotonically away from the axis in both directions.
        @test issorted(speeds[1:min_idx]; rev=true)
        @test issorted(speeds[min_idx:end])
        @test maximum(speeds) > 1520.0

        # The two assertions above are insensitive to ε, because the ε term vanishes
        # identically at the axis. Pinning the endpoints is what guards ε itself:
        # ε = 0.00737 -> c(0) = 1548.52, whereas ε = 0.0080 -> c(0) = 1541.86.
        @test speeds[1] ≈ 1548.5210151736783 atol = 1e-6
        @test speeds[end] ≈ 1551.9107368195284 atol = 1e-6
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

# ---------------------------------------------------------------------------------------------
# Regression tests for the five latent bugs fixed in plan task 2.5. Each one fails on the code as
# it stood before that commit; the failure mode is noted so a future "simplification" that
# reintroduces the bug is recognisable from the test output alone.
# ---------------------------------------------------------------------------------------------

@testitem "B1: showing a density profile does not throw" begin
    using Kraken

    ρ = SampledDensity([0.0, 50.0, 100.0], [1000.0, 1020.0, 1030.0])

    # Before the fix, the show method printed a nonexistent `.type` field, so *any* display of a
    # density profile — including `env.ρ` at the REPL — threw an ErrorException.
    s = sprint(show, ρ)
    @test occursin("SampledDensity1D", s)
    @test occursin("3 points", s)
    @test !occursin("type", s)

    # And it must survive the paths that call show implicitly.
    @test sprint(print, ρ) == s
    @test !isempty(sprint(show, UnderwaterEnv(pekeris_env()...).ρ))
end

@testitem "B2: pressure_f evaluates the modal functions" begin
    using Kraken

    env = PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0)
    freq = 100.0
    krs = Kraken.find_kr(env, freq)
    @test !isempty(krs)

    # Before the fix this called a 5-argument `get_modal_function`, which does not exist — every
    # call raised MethodError. Only the 3-argument method (returning a closure) is defined.
    p = pressure_f(env, krs, freq, 1000.0, 25.0, 50.0)
    @test p isa Complex
    @test isfinite(real(p)) && isfinite(imag(p))
    @test abs(p) > 0

    # Documented short-circuits still hold.
    @test pressure_f(env, krs, 0.0, 1000.0, 25.0, 50.0) == 0.0 + 0.0im
    @test pressure_f(env, eltype(krs)[], freq, 1000.0, 25.0, 50.0) == 0.0 + 0.0im

    # Pressure falls off with range (the 1/sqrt(8πr) geometric factor dominates).
    @test abs(pressure_f(env, krs, freq, 10_000.0, 25.0, 50.0)) < abs(pressure_f(env, krs, freq, 1_000.0, 25.0, 50.0))
end

@testitem "B3: det_sturm has no inert stop_at_k option" begin
    using Kraken

    env = UnderwaterEnv(pekeris_env()...)
    props = AcousticProblemProperties(env, 100.0)
    cache = AcousticProblemCache(env, props)

    # The keyword used to exist and do nothing at all: its body was the expression statement
    # `p2, mode_count` rather than a `return`, so the "early exit" ran the loop to completion and
    # returned the full-length result. Rather than leave a silently-inert option, the parameter is
    # gone — passing it is now an error the caller can see.
    @test_throws MethodError det_sturm(0.41, env, props, cache; stop_at_k=3)
    @test_throws MethodError det_sturm(0.41, env, props, cache; return_det=true)

    # The one supported keyword, `scale`, is load-bearing rather than cosmetic. Pekeris at 100 Hz
    # has 5 modes and 125 mesh points; without rescaling the Sturm sequence underflows to exactly
    # 0.0 partway down and the mode count comes out wrong. Pin both behaviours so nobody "cleans
    # up" scale_const on the assumption that it only guards pathological inputs.
    krs = find_kr(env, props, cache)
    d_scaled, n_scaled = det_sturm(0.41, env, props, cache; scale=true)
    d_raw, n_raw = det_sturm(0.41, env, props, cache; scale=false)

    @test n_scaled == count(>(0.41), krs)   # 2 modes lie above kr = 0.41
    @test d_scaled != 0
    @test d_raw == 0.0                      # underflowed
    @test n_raw < n_scaled                  # and lost a mode as a result
end

@testitem "B4: bisection indexes kLeft/kRight in bounds" begin
    using Kraken

    # `kLeft[Δn]` and `kRight[Δn + 1]` are indexed with no guard. The branch they sit in is only
    # reached when Δn >= mm >= 1, and Δn <= n_max, so both are in bounds — but that invariant is
    # implicit in the loop structure and easy to break. Sweep environments and frequencies wide
    # enough to exercise mode counts from 0 to ~100 and assert nothing throws.
    for (name, envtuple) in (
        ("pekeris", pekeris_env()),
        ("pekeris shallow", pekeris_env(depth=10.0)),
        ("pekeris deep", pekeris_env(depth=2000.0)),
        ("pekeris fast bottom", pekeris_env(cb=2500.0)),
        ("one layer", one_layer_env()),
        ("one layer slope", one_layer_slope_env()),
        ("two layer slope", two_layer_slope_env()),
        ("munk", munk_env()),
    )
        env = UnderwaterEnv(envtuple...)
        for freq in (5.0, 25.0, 100.0, 400.0)
            props = AcousticProblemProperties(env, freq)
            cache = AcousticProblemCache(env, props)
            intervals = bisection(env, props, cache)
            @test isnothing(intervals) || (intervals isa Matrix && size(intervals, 2) == 2)
            if intervals isa Matrix
                # Every bracket must be ordered and inside the trapped band.
                ω = 2π * freq
                @test all(intervals[:, 1] .<= intervals[:, 2])
                @test all(intervals .>= ω / env.cb)
            end
        end
    end
end

@testitem "B5: both UnderwaterEnv constructors agree on depth" begin
    using Kraken

    for envtuple in (pekeris_env(), one_layer_env(), one_layer_slope_env(), two_layer_slope_env(), munk_env())
        ssp, layers, sspHS = envtuple
        direct = UnderwaterEnv(ssp, layers, sspHS)
        viafortran = UnderwaterEnv(UnderwaterEnvFORTRAN(ssp, layers, sspHS))

        # The two constructors used to read `depth` from different places — layers[end, 3] here and
        # ssp[end, 1] there — so they could disagree. Both now use layers[end, 3].
        @test direct.depth == viafortran.depth
        @test direct.depth == layers[end, 3]
        @test direct.layer_depth == viafortran.layer_depth
        @test direct.h_vec == viafortran.h_vec
    end

    # Make the disagreement observable: an ssp table sampled past the last layer boundary. `layers`
    # is authoritative, so `depth` follows it and not the deeper ssp row.
    ssp, layers, sspHS = pekeris_env(depth=100.0)
    ssp_long = vcat(ssp, [150.0 1500.0 0.0 1000.0 0.0 0.0])
    @test UnderwaterEnv(ssp_long, layers, sspHS).depth == 100.0
    @test UnderwaterEnv(UnderwaterEnvFORTRAN(ssp_long, layers, sspHS)).depth == 100.0
end
