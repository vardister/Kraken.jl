using Test
using Kraken
using BenchmarkTools

@testset "Performance Tests" begin
    
    # Performance thresholds (in seconds) - adjust based on typical hardware
    const PERFORMANCE_THRESHOLDS = Dict(
        :pekeris_100hz => 5.0,      # Basic Pekeris at 100 Hz
        :pekeris_500hz => 10.0,     # Higher frequency
        :multilayer => 15.0,        # More complex environment
        :munk_profile => 20.0,      # Deep water with many modes
        :find_kr_basic => 2.0,      # Core numerical function
        :inverse_iteration => 3.0,   # Mode shape calculation
        :ad_gradient => 8.0,        # Automatic differentiation
    )
    
    @testset "Basic Performance Benchmarks" begin
        
        @testset "Standard Pekeris Waveguide" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            
            # Benchmark basic case
            result = @timed kraken_jl(env, 100.0)
            elapsed_time = result.time
            sol = result.value
            
            @test elapsed_time < PERFORMANCE_THRESHOLDS[:pekeris_100hz]
            @test !isempty(sol.kr)  # Should produce results
            
            println("Pekeris 100Hz: $(elapsed_time:.3f)s (threshold: $(PERFORMANCE_THRESHOLDS[:pekeris_100hz])s)")
        end
        
        @testset "Higher Frequency Performance" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            
            # Test higher frequency (more modes)
            result = @timed kraken_jl(env, 500.0)
            elapsed_time = result.time
            sol = result.value
            
            @test elapsed_time < PERFORMANCE_THRESHOLDS[:pekeris_500hz]
            @test length(sol.kr) > 5  # Should have many modes at high frequency
            
            println("Pekeris 500Hz: $(elapsed_time:.3f)s (threshold: $(PERFORMANCE_THRESHOLDS[:pekeris_500hz])s)")
        end
        
        @testset "Multi-layer Environment" begin
            ssp, layers, sspHS = one_layer_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            
            result = @timed kraken_jl(env, 100.0)
            elapsed_time = result.time
            sol = result.value
            
            @test elapsed_time < PERFORMANCE_THRESHOLDS[:multilayer]
            @test !isempty(sol.kr)
            
            println("Multi-layer: $(elapsed_time:.3f)s (threshold: $(PERFORMANCE_THRESHOLDS[:multilayer])s)")
        end
        
        @testset "Deep Water (Munk Profile)" begin
            ssp, layers, sspHS = munk_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            
            result = @timed kraken_jl(env, 25.0)  # Lower frequency for deep water
            elapsed_time = result.time
            sol = result.value
            
            @test elapsed_time < PERFORMANCE_THRESHOLDS[:munk_profile]
            @test length(sol.kr) > 10  # Deep water should support many modes
            
            println("Munk profile: $(elapsed_time:.3f)s (threshold: $(PERFORMANCE_THRESHOLDS[:munk_profile])s)")
        end
    end
    
    @testset "Core Function Performance" begin
        
        @testset "find_kr Performance" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            props = AcousticProblemProperties(env, 100.0)
            cache = AcousticProblemCache(env, props)
            
            result = @timed find_kr(env, props, cache)
            elapsed_time = result.time
            krs = result.value
            
            @test elapsed_time < PERFORMANCE_THRESHOLDS[:find_kr_basic]
            @test !isempty(krs)
            
            println("find_kr: $(elapsed_time:.3f)s (threshold: $(PERFORMANCE_THRESHOLDS[:find_kr_basic])s)")
        end
        
        @testset "Inverse Iteration Performance" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            props = AcousticProblemProperties(env, 100.0)
            cache = AcousticProblemCache(env, props)
            
            # Get wavenumbers first
            krs = find_kr(env, props, cache)
            
            if !isempty(krs)
                result = @timed inverse_iteration(krs[1:min(3, length(krs))], env, props, cache)
                elapsed_time = result.time
                
                @test elapsed_time < PERFORMANCE_THRESHOLDS[:inverse_iteration]
                
                println("Inverse iteration: $(elapsed_time:.3f)s (threshold: $(PERFORMANCE_THRESHOLDS[:inverse_iteration])s)")
            end
        end
    end
    
    @testset "Automatic Differentiation Performance" begin
        
        @testset "Parameter Gradient Performance" begin
            using ForwardDiff
            
            function kr_multi_param(params)
                c0, cb, depth = params
                ssp, layers, sspHS = pekeris_env(c0=c0, cb=cb, depth=depth)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end
            
            params_base = [1500.0, 1600.0, 100.0]
            
            result = @timed ForwardDiff.gradient(kr_multi_param, params_base)
            elapsed_time = result.time
            grad = result.value
            
            @test elapsed_time < PERFORMANCE_THRESHOLDS[:ad_gradient]
            @test length(grad) == 3
            @test !any(isnan.(grad))
            
            println("AD gradient: $(elapsed_time:.3f)s (threshold: $(PERFORMANCE_THRESHOLDS[:ad_gradient])s)")
        end
        
        @testset "Frequency Derivative Performance" begin
            using ForwardDiff
            
            function kr_vs_freq(freq)
                ssp, layers, sspHS = pekeris_env()
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, freq)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end
            
            result = @timed ForwardDiff.derivative(kr_vs_freq, 100.0)
            elapsed_time = result.time
            deriv = result.value
            
            @test elapsed_time < PERFORMANCE_THRESHOLDS[:ad_gradient]  # Use same threshold
            @test !isnan(deriv)
            
            println("Frequency derivative: $(elapsed_time:.3f)s")
        end
    end
    
    @testset "Memory Usage Tests" begin
        
        @testset "Memory Allocation Regression" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            
            # Measure memory allocation
            stats = @timed kraken_jl(env, 100.0)
            allocated_bytes = stats.bytes
            
            # Basic sanity check - shouldn't allocate excessive memory
            # (adjust threshold based on typical usage)
            max_expected_mb = 100  # 100 MB maximum
            allocated_mb = allocated_bytes / (1024^2)
            
            @test allocated_mb < max_expected_mb
            
            println("Memory allocated: $(allocated_mb:.2f) MB (max: $(max_expected_mb) MB)")
        end
        
        @testset "Repeated Calls Memory Stability" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            
            # Run multiple times to check for memory leaks
            initial_stats = @timed kraken_jl(env, 100.0)
            
            for i in 1:5
                stats = @timed kraken_jl(env, 100.0)
            end
            
            final_stats = @timed kraken_jl(env, 100.0)
            
            # Memory usage shouldn't grow significantly
            memory_growth = final_stats.bytes / initial_stats.bytes
            @test memory_growth < 2.0  # Less than 2x growth
            
            println("Memory growth ratio: $(memory_growth:.2f)x")
        end
    end
    
    @testset "Scaling Performance" begin
        
        @testset "Frequency Scaling" begin
            ssp, layers, sspHS = pekeris_env()
            env = UnderwaterEnv(ssp, layers, sspHS)
            
            frequencies = [50.0, 100.0, 200.0, 400.0]
            times = Float64[]
            n_modes = Int[]
            
            for freq in frequencies
                result = @timed kraken_jl(env, freq)
                push!(times, result.time)
                push!(n_modes, length(result.value.kr))
            end
            
            # Performance should scale reasonably with frequency
            # (not exponentially)
            for i in 2:length(times)
                scaling_factor = times[i] / times[1]
                freq_ratio = frequencies[i] / frequencies[1]
                
                # Time shouldn't scale worse than quadratically with frequency
                @test scaling_factor < freq_ratio^2
            end
            
            println("Frequency scaling:")
            for (f, t, n) in zip(frequencies, times, n_modes)
                println("  $(f) Hz: $(t:.3f)s ($(n) modes)")
            end
        end
        
        @testset "Depth Scaling" begin
            depths = [50.0, 100.0, 200.0, 500.0]
            times = Float64[]
            
            for depth in depths
                ssp, layers, sspHS = pekeris_env(depth=depth)
                env = UnderwaterEnv(ssp, layers, sspHS)
                
                result = @timed kraken_jl(env, 100.0)
                push!(times, result.time)
            end
            
            # Deeper water takes longer but shouldn't be exponential
            for i in 2:length(times)
                scaling_factor = times[i] / times[1]
                depth_ratio = depths[i] / depths[1]
                
                # Time shouldn't scale worse than quadratically with depth
                @test scaling_factor < depth_ratio^2
            end
            
            println("Depth scaling:")
            for (d, t) in zip(depths, times)
                println("  $(d) m: $(t:.3f)s")
            end
        end
    end
    
    @testset "Comparison Benchmarks" begin
        
        @testset "Julia vs Fortran Performance" begin
            ssp, layers, sspHS = pekeris_env()
            env_julia = UnderwaterEnv(ssp, layers, sspHS)
            
            zrc = [0.0, 100.0]
            zsr = 50.0
            env_fortran = EnvKRAKEN(ssp, layers, sspHS, zrc, zsr)
            
            freq = 100.0
            n_modes = 5
            
            # Benchmark both implementations
            julia_result = @timed kraken_jl(env_julia, freq)
            fortran_result = @timed kraken(env_fortran, freq; n_modes=n_modes)
            
            julia_time = julia_result.time
            fortran_time = fortran_result.time
            
            # Both should complete in reasonable time
            @test julia_time < 10.0
            @test fortran_time < 10.0
            
            # Results should be consistent
            julia_kr = julia_result.value.kr[1]
            fortran_kr = fortran_result.value["kr_real"][1]
            @test julia_kr ≈ fortran_kr atol=1e-6
            
            println("Julia vs Fortran:")
            println("  Julia: $(julia_time:.3f)s")
            println("  Fortran: $(fortran_time:.3f)s")
            println("  Ratio: $(julia_time/fortran_time:.2f)x")
        end
    end
    
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
    end
end
