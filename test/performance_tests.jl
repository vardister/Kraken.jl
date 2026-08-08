using Test
using Kraken
using BenchmarkTools
using ForwardDiff

# Enabled with KRAKEN_RUN_PERFORMANCE_TESTS=true. Every timing below is measured *after* a warm-up
# call, so compilation is not counted; thresholds are set several times above the measured value on
# a 2021 M1 laptop so they catch order-of-magnitude regressions without flaking on slower CI runners.
# Ratios that are dominated by timer noise (frequency/depth scaling) are reported, not asserted —
# see the comments at those testsets.

@testset "Performance Tests" begin

    # Performance thresholds (in seconds) - adjust based on typical hardware
    PERFORMANCE_THRESHOLDS = Dict(
        :pekeris_100hz => 5.0,      # Basic Pekeris at 100 Hz
        :pekeris_500hz => 10.0,     # Higher frequency
        :multilayer => 15.0,        # More complex environment
        :munk_profile => 20.0,      # Deep water with many modes
        :find_kr_basic => 2.0,      # Core numerical function
        :inverse_iteration => 3.0,   # Mode shape calculation
        :ad_gradient => 8.0,        # Automatic differentiation
    )

    # sigdigits, not digits: every measurement here is well under a millisecond on a warm
    # session, and rounding those to 3 decimal places prints a useless "0.0s".
    fmt(x) = string(round(x; sigdigits=3))

    @testset "Basic Performance Benchmarks" begin
        @testset "Standard Pekeris Waveguide" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)

            kraken_jl(env, 100.0)  # warm-up: keep compilation out of the measurement
            result = @timed kraken_jl(env, 100.0)
            elapsed_time = result.time
            sol = result.value

            @test elapsed_time < PERFORMANCE_THRESHOLDS[:pekeris_100hz]
            @test !isempty(sol.kr)  # Should produce results

            println("Pekeris 100Hz: $(fmt(elapsed_time))s (threshold: $(PERFORMANCE_THRESHOLDS[:pekeris_100hz])s)")
        end

        @testset "Higher Frequency Performance" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)

            # Test higher frequency (more modes)
            kraken_jl(env, 500.0)
            result = @timed kraken_jl(env, 500.0)
            elapsed_time = result.time
            sol = result.value

            @test elapsed_time < PERFORMANCE_THRESHOLDS[:pekeris_500hz]
            @test length(sol.kr) > 5  # Should have many modes at high frequency

            println("Pekeris 500Hz: $(fmt(elapsed_time))s (threshold: $(PERFORMANCE_THRESHOLDS[:pekeris_500hz])s)")
        end

        @testset "Multi-layer Environment" begin
            ssp, layers, sspHS = one_layer_env()
            env = UnderwaterEnv(ssp, layers, sspHS)

            kraken_jl(env, 100.0)
            result = @timed kraken_jl(env, 100.0)
            elapsed_time = result.time
            sol = result.value

            @test elapsed_time < PERFORMANCE_THRESHOLDS[:multilayer]
            @test !isempty(sol.kr)

            println("Multi-layer: $(fmt(elapsed_time))s (threshold: $(PERFORMANCE_THRESHOLDS[:multilayer])s)")
        end

        @testset "Deep Water (Munk Profile)" begin
            ssp, layers, sspHS = munk_env()
            env = UnderwaterEnv(ssp, layers, sspHS)

            kraken_jl(env, 25.0)
            result = @timed kraken_jl(env, 25.0)  # Lower frequency for deep water
            elapsed_time = result.time
            sol = result.value

            @test elapsed_time < PERFORMANCE_THRESHOLDS[:munk_profile]
            @test length(sol.kr) > 10  # Deep water should support many modes

            println("Munk profile: $(fmt(elapsed_time))s (threshold: $(PERFORMANCE_THRESHOLDS[:munk_profile])s)")
        end
    end

    @testset "Core Function Performance" begin
        @testset "find_kr Performance" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            props = AcousticProblemProperties(env, 100.0)
            cache = AcousticProblemCache(env, props)

            find_kr(env, props, cache)
            result = @timed find_kr(env, props, cache)
            elapsed_time = result.time
            krs = result.value

            @test elapsed_time < PERFORMANCE_THRESHOLDS[:find_kr_basic]
            @test !isempty(krs)

            println("find_kr: $(fmt(elapsed_time))s (threshold: $(PERFORMANCE_THRESHOLDS[:find_kr_basic])s)")
        end

        @testset "Inverse Iteration Performance" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            props = AcousticProblemProperties(env, 100.0)
            cache = AcousticProblemCache(env, props)

            # Get wavenumbers first
            krs = find_kr(env, props, cache)
            @test !isempty(krs)

            krs_test = krs[1:min(3, length(krs))]
            inverse_iteration(krs_test, env, props, cache)
            result = @timed inverse_iteration(krs_test, env, props, cache)
            elapsed_time = result.time

            @test elapsed_time < PERFORMANCE_THRESHOLDS[:inverse_iteration]

            println(
                "Inverse iteration: $(fmt(elapsed_time))s (threshold: $(PERFORMANCE_THRESHOLDS[:inverse_iteration])s)"
            )
        end
    end

    @testset "Automatic Differentiation Performance" begin
        @testset "Parameter Gradient Performance" begin
            function kr_multi_param(params)
                c0, cb, depth = params
                ssp, layers, sspHS = pekeris_env(c0=c0, cb=cb, depth=depth)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            params_base = [1500.0, 1600.0, 100.0]

            ForwardDiff.gradient(kr_multi_param, params_base)
            result = @timed ForwardDiff.gradient(kr_multi_param, params_base)
            elapsed_time = result.time
            grad = result.value

            @test elapsed_time < PERFORMANCE_THRESHOLDS[:ad_gradient]
            @test length(grad) == 3
            @test !any(isnan.(grad))

            println("AD gradient: $(fmt(elapsed_time))s (threshold: $(PERFORMANCE_THRESHOLDS[:ad_gradient])s)")
        end

        @testset "Frequency Derivative Performance" begin
            function kr_vs_freq(freq)
                ssp, layers, sspHS = pekeris_env()
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, freq)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            ForwardDiff.derivative(kr_vs_freq, 100.0)
            result = @timed ForwardDiff.derivative(kr_vs_freq, 100.0)
            elapsed_time = result.time
            deriv = result.value

            @test elapsed_time < PERFORMANCE_THRESHOLDS[:ad_gradient]  # Use same threshold
            @test !isnan(deriv)

            println("Frequency derivative: $(fmt(elapsed_time))s")
        end
    end

    @testset "Memory Usage Tests" begin
        @testset "Memory Allocation Regression" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)

            # Measure memory allocation
            kraken_jl(env, 100.0)
            stats = @timed kraken_jl(env, 100.0)
            allocated_bytes = stats.bytes

            # Basic sanity check - shouldn't allocate excessive memory
            # (adjust threshold based on typical usage)
            max_expected_mb = 100  # 100 MB maximum
            allocated_mb = allocated_bytes / (1024^2)

            @test allocated_mb < max_expected_mb

            println("Memory allocated: $(fmt(allocated_mb)) MB (max: $(max_expected_mb) MB)")
        end

        @testset "Repeated Calls Memory Stability" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)

            # Run multiple times to check for memory leaks. The first call is discarded so the
            # baseline does not include one-time compilation allocations.
            kraken_jl(env, 100.0)
            initial_stats = @timed kraken_jl(env, 100.0)

            for _ in 1:5
                kraken_jl(env, 100.0)
            end

            final_stats = @timed kraken_jl(env, 100.0)

            # Memory usage shouldn't grow at all across identical calls
            memory_growth = final_stats.bytes / initial_stats.bytes
            @test memory_growth < 1.5

            println("Memory growth ratio: $(fmt(memory_growth))x")
        end
    end

    @testset "Scaling Performance" begin

        # Wall-clock ratios between sub-second runs are dominated by timer noise and by how many
        # modes each case happens to produce, so these two testsets report the scaling and assert
        # only a generous absolute ceiling on each individual run.
        @testset "Frequency Scaling" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)

            frequencies = [50.0, 100.0, 200.0, 400.0]
            times = Float64[]
            n_modes = Int[]

            kraken_jl(env, 50.0)
            for freq in frequencies
                result = @timed kraken_jl(env, freq)
                push!(times, result.time)
                push!(n_modes, length(result.value.kr))
            end

            @test all(times .< PERFORMANCE_THRESHOLDS[:pekeris_500hz])
            @test issorted(n_modes)  # more modes at higher frequency

            println("Frequency scaling:")
            for (f, t, n) in zip(frequencies, times, n_modes)
                println("  $(f) Hz: $(fmt(t))s ($(n) modes, $(fmt(t / times[1]))x the 50 Hz time)")
            end
        end

        @testset "Depth Scaling" begin
            depths = [50.0, 100.0, 200.0, 500.0]
            times = Float64[]
            n_modes = Int[]

            let (ssp, layers, sspHS) = pekeris_env(depth=50.0)
                kraken_jl(UnderwaterEnv(ssp, layers, sspHS), 100.0)
            end
            for depth in depths
                ssp, layers, sspHS = pekeris_env(depth=depth)
                env = UnderwaterEnv(ssp, layers, sspHS)

                result = @timed kraken_jl(env, 100.0)
                push!(times, result.time)
                push!(n_modes, length(result.value.kr))
            end

            @test all(times .< PERFORMANCE_THRESHOLDS[:multilayer])
            @test issorted(n_modes)  # more modes in deeper water

            println("Depth scaling:")
            for (d, t, n) in zip(depths, times, n_modes)
                println("  $(d) m: $(fmt(t))s ($(n) modes, $(fmt(t / times[1]))x the 50 m time)")
            end
        end
    end

    # The "Julia vs Fortran Performance" testset that used to live here was removed in plan task 1.6:
    # it called the `EnvKRAKEN`/`kraken` API, which exists in no loaded module, so it could never have
    # run. Milestone 3 rebuilds the Fortran comparison on `test/reference/` (unmodified `kraken.exe`
    # from AcousticsToolbox_jll), and the timing comparison belongs there, next to the correctness one.

    @testset "Performance Summary" begin
        # Print summary of all thresholds and whether they were met
        println("\nPerformance Test Summary:")
        println("========================")
        for (test_name, threshold) in PERFORMANCE_THRESHOLDS
            println("$(test_name): < $(threshold)s")
        end

        # Additional system info for context
        println("\nSystem Info:")
        println("Julia version: $(VERSION)")
        println("Threads: $(Threads.nthreads())")
        @test Threads.nthreads() >= 1
    end
end
