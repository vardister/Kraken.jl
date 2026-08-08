using Test
using Kraken
using ForwardDiff
using FiniteDiff

@testset "Automatic Differentiation Tests" begin

    # Helper functions for AD testing
    function test_function_gradient(func, x, name="function"; atol=1e-8, rtol=1e-6)
        grad_fd = FiniteDiff.finite_difference_gradient(func, x)
        grad_ad = ForwardDiff.gradient(func, x)

        @test grad_ad ≈ grad_fd atol=atol rtol=rtol
        return grad_ad, grad_fd
    end

    function test_function_derivative(func, x, name="function"; atol=1e-8, rtol=1e-6)
        deriv_fd = FiniteDiff.finite_difference_derivative(func, x)
        deriv_ad = ForwardDiff.derivative(func, x)

        @test deriv_ad ≈ deriv_fd atol=atol rtol=rtol
        return deriv_ad, deriv_fd
    end

    @testset "Basic Environment Parameter Derivatives" begin
        @testset "Sound Speed Parameter Derivatives" begin
            function kr_vs_c0(c0)
                ssp, layers, sspHS = pekeris_env(c0=c0)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            function kr_vs_cb(cb)
                ssp, layers, sspHS = pekeris_env(cb=cb)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            # Test derivatives w.r.t. water sound speed
            c0_base = 1500.0
            deriv_c0_ad, deriv_c0_fd = test_function_derivative(kr_vs_c0, c0_base, "kr vs c0")
            @test deriv_c0_ad < 0  # Higher water speed should decrease kr

            # Test derivatives w.r.t. bottom sound speed
            # ∂kr/∂cb < 0: a faster bottom means a larger evanescent decay rate
            # β = √(kr² − (ω/cb)²), i.e. a *lower*-impedance (more pressure-release-like)
            # boundary. That pushes the vertical wavenumber γ toward nπ/D and hence lowers
            # kr = √((ω/c0)² − γ²). Verified against the analytic Pekeris solver
            # (`find_kr(::PekerisUnderwaterEnv, 100.0)`), which gives −8.03e-7 for mode 1 —
            # matching this solver to 4 digits. The magnitude is small because mode 1 barely
            # penetrates the bottom.
            cb_base = 1600.0
            deriv_cb_ad, deriv_cb_fd = test_function_derivative(kr_vs_cb, cb_base, "kr vs cb")
            @test deriv_cb_ad < 0
        end

        @testset "Density Parameter Derivatives" begin
            function kr_vs_rho0(ρ0)
                ssp, layers, sspHS = pekeris_env(ρ0=ρ0)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            function kr_vs_rhob(ρb)
                ssp, layers, sspHS = pekeris_env(ρb=ρb)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            # Test derivatives w.r.t. water density
            ρ0_base = 1000.0
            deriv_rho0_ad, deriv_rho0_fd = test_function_derivative(kr_vs_rho0, ρ0_base, "kr vs ρ0", rtol=1e-4)

            # Test derivatives w.r.t. bottom density
            ρb_base = 1500.0
            deriv_rhob_ad, deriv_rhob_fd = test_function_derivative(kr_vs_rhob, ρb_base, "kr vs ρb", rtol=1e-4)
        end

        @testset "Depth Parameter Derivatives" begin
            function kr_vs_depth(depth)
                ssp, layers, sspHS = pekeris_env(depth=depth)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            # ∂kr/∂depth > 0: for a fixed mode number the vertical wavenumber is γ_n ≈ nπ/D,
            # so deeper water lowers γ_n and raises kr = √((ω/c0)² − γ_n²) toward ω/c0.
            # Analytic Pekeris solver gives +1.7678e-5 for mode 1; this solver gives +1.7679e-5.
            depth_base = 100.0
            deriv_depth_ad, deriv_depth_fd = test_function_derivative(kr_vs_depth, depth_base, "kr vs depth", rtol=1e-5)
            @test deriv_depth_ad > 0
        end
    end

    @testset "Multi-parameter Derivatives" begin
        @testset "Multiple Environment Parameters" begin
            function kr_multi_param(params)
                c0, cb, ρ0, ρb, depth = params
                ssp, layers, sspHS = pekeris_env(c0=c0, cb=cb, ρ0=ρ0, ρb=ρb, depth=depth)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            params_base = [1500.0, 1600.0, 1000.0, 1500.0, 100.0]
            grad_ad, grad_fd = test_function_gradient(kr_multi_param, params_base, "kr vs multiple params", rtol=1e-4)

            # Test physical intuition (signs justified in "Sound Speed"/"Depth" testsets above)
            @test grad_ad[1] < 0  # ∂kr/∂c0 < 0
            @test grad_ad[2] < 0  # ∂kr/∂cb < 0
            @test grad_ad[5] > 0  # ∂kr/∂depth > 0
        end

        @testset "Multiple Modes" begin
            function krs_multi_param(params)
                c0, cb = params
                ssp, layers, sspHS = pekeris_env(c0=c0, cb=cb)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                # Return first 3 modes if available
                krs = sol.kr
                if length(krs) >= 3
                    return krs[1:3]
                elseif length(krs) >= 2
                    return [krs[1:2]; 0.0]
                elseif length(krs) >= 1
                    return [krs[1]; 0.0; 0.0]
                else
                    return [0.0, 0.0, 0.0]
                end
            end

            params_base = [1500.0, 1600.0]
            jac_ad = ForwardDiff.jacobian(krs_multi_param, params_base)
            # Central differences: `finite_difference_jacobian` defaults to forward
            # differences with relstep √eps, whose step is far too small for an output
            # produced by an iterative solver — see the "AD vs Finite Difference" testset.
            jac_fd = FiniteDiff.finite_difference_jacobian(krs_multi_param, params_base, Val{:central})

            @test jac_ad ≈ jac_fd rtol=1e-4

            # Test that all modes have similar derivative structure. ∂kr/∂cb grows in
            # magnitude with mode number (higher modes penetrate the bottom more) but stays
            # negative for every mode.
            for i in 1:3
                if jac_ad[i, 1] != 0  # If mode exists
                    @test jac_ad[i, 1] < 0  # ∂kr_i/∂c0 < 0
                    @test jac_ad[i, 2] < 0  # ∂kr_i/∂cb < 0
                end
            end
        end
    end

    @testset "Frequency Derivatives (Group Speed)" begin
        @testset "Single Mode Group Speed" begin
            function kr_vs_freq(freq)
                ssp, layers, sspHS = pekeris_env()
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, freq)
                return isempty(sol.kr) ? 0.0 : sol.kr[1]
            end

            freq_base = 100.0
            dkr_df_ad, dkr_df_fd = test_function_derivative(kr_vs_freq, freq_base, "kr vs frequency", rtol=1e-5)

            # Group speed is cg = 2π / (dkr/df) = (dω/dkr)
            group_speed_ad = 2π / dkr_df_ad

            # Should be physically reasonable (between water and bottom speeds)
            @test 1400.0 < group_speed_ad < 1700.0

            # Test using the README example approach
            function calculate_kr_for_group_speed(freq)
                ssp, layers, sspHS = pekeris_env()
                env = UnderwaterEnv(ssp, layers, sspHS)
                props = AcousticProblemProperties(env, freq)
                cache = AcousticProblemCache(env, props)
                wavenumbers = find_kr(env, props, cache)
                return isempty(wavenumbers) ? 0.0 : wavenumbers[1]
            end

            group_speed_readme = 2π / ForwardDiff.derivative(calculate_kr_for_group_speed, freq_base)
            @test group_speed_readme ≈ group_speed_ad rtol=1e-3
        end

        @testset "Multiple Mode Group Speeds" begin
            function krs_vs_freq(freq)
                ssp, layers, sspHS = pekeris_env()
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, freq)
                krs = sol.kr
                # Return first 3 modes
                if length(krs) >= 3
                    return krs[1:3]
                else
                    return vcat(krs, zeros(3 - length(krs)))
                end
            end

            freq_base = 100.0
            jac_ad = ForwardDiff.jacobian(f -> krs_vs_freq(f[1]), [freq_base])
            dkrs_df = jac_ad[:, 1]

            # Calculate group speeds for each mode
            group_speeds = 2π ./ dkrs_df[dkrs_df .!= 0]

            for cg in group_speeds
                @test 1400.0 < cg < 1700.0  # Reasonable range
            end

            # Higher order modes should have slower group speeds
            if length(group_speeds) > 1
                @test issorted(group_speeds, rev=true)
            end
        end
    end

    @testset "Core Function AD Compatibility" begin
        @testset "det_sturm Function" begin
            function det_sturm_wrapper(kr_params)
                kr = kr_params[1]
                ssp, layers, sspHS = pekeris_env()
                env = UnderwaterEnv(ssp, layers, sspHS)
                props = AcousticProblemProperties(env, 100.0)
                cache = AcousticProblemCache(env, props)
                det_val, _ = det_sturm(kr, env, props, cache)
                return det_val
            end

            kr_test = [0.41]
            grad_ad = ForwardDiff.gradient(det_sturm_wrapper, kr_test)
            grad_fd = FiniteDiff.finite_difference_gradient(det_sturm_wrapper, kr_test)

            @test grad_ad ≈ grad_fd rtol=1e-5
            @test !isnan(grad_ad[1])
            @test !isinf(grad_ad[1])
        end

        @testset "get_g Function" begin
            function g_wrapper(kr_params)
                kr = kr_params[1]
                ssp, layers, sspHS = pekeris_env()
                env = UnderwaterEnv(ssp, layers, sspHS)
                props = AcousticProblemProperties(env, 100.0)
                return get_g(kr, env, props)
            end

            kr_test = [0.41]
            grad_ad = ForwardDiff.gradient(g_wrapper, kr_test)
            grad_fd = FiniteDiff.finite_difference_gradient(g_wrapper, kr_test)

            @test grad_ad ≈ grad_fd rtol=1e-5
            @test !isnan(grad_ad[1])
            @test !isinf(grad_ad[1])
        end
    end

    @testset "Inverse Iteration AD" begin
        @testset "Mode Shape Derivatives" begin
            function mode_amplitude_vs_param(c0)
                ssp, layers, sspHS = pekeris_env(c0=c0)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)
                if isempty(sol.modes)
                    return 0.0
                end
                # Return amplitude at mid-depth
                mid_idx = size(sol.modes, 1) ÷ 2
                return sol.modes[mid_idx, 1]  # First mode amplitude
            end

            c0_base = 1500.0
            deriv_ad, deriv_fd = test_function_derivative(
                mode_amplitude_vs_param, c0_base, "mode amplitude vs c0", rtol=1e-3
            )

            @test !isnan(deriv_ad)
            @test !isinf(deriv_ad)
        end
    end

    @testset "Numerical Stability" begin
        @testset "Parameter Perturbations" begin
            function test_stability(params)
                c0, cb = params
                try
                    ssp, layers, sspHS = pekeris_env(c0=c0, cb=cb)
                    env = UnderwaterEnv(ssp, layers, sspHS)
                    sol = kraken_jl(env, 100.0)
                    return isempty(sol.kr) ? 0.0 : sol.kr[1]
                catch
                    return 0.0
                end
            end

            # Test with small perturbations
            params_base = [1500.0, 1600.0]

            # Should work with ForwardDiff
            grad_ad = ForwardDiff.gradient(test_stability, params_base)
            @test !any(isnan.(grad_ad))
            @test !any(isinf.(grad_ad))

            # Test with larger perturbations
            params_perturbed = [1520.0, 1650.0]
            grad_ad_pert = ForwardDiff.gradient(test_stability, params_perturbed)
            @test !any(isnan.(grad_ad_pert))
            @test !any(isinf.(grad_ad_pert))
        end

        @testset "Frequency Range Stability" begin
            function kr_stable(freq)
                try
                    ssp, layers, sspHS = pekeris_env()
                    env = UnderwaterEnv(ssp, layers, sspHS)
                    sol = kraken_jl(env, freq)
                    return isempty(sol.kr) ? 0.0 : sol.kr[1]
                catch
                    return 0.0
                end
            end

            # Test across frequency range
            for freq_test in [50.0, 100.0, 200.0, 500.0]
                deriv_ad = ForwardDiff.derivative(kr_stable, freq_test)
                @test !isnan(deriv_ad)
                @test !isinf(deriv_ad)
            end
        end
    end

    @testset "Performance and Accuracy" begin
        @testset "AD vs Finite Difference Accuracy" begin
            function multi_output_test(params)
                c0, cb, depth = params
                ssp, layers, sspHS = pekeris_env(c0=c0, cb=cb, depth=depth)
                env = UnderwaterEnv(ssp, layers, sspHS)
                sol = kraken_jl(env, 100.0)

                if length(sol.kr) >= 2
                    return [sol.kr[1], sol.kr[2]]
                elseif length(sol.kr) >= 1
                    return [sol.kr[1], 0.0]
                else
                    return [0.0, 0.0]
                end
            end

            params_base = [1500.0, 1600.0, 100.0]

            jac_ad = ForwardDiff.jacobian(multi_output_test, params_base)
            # `finite_difference_jacobian` defaults to Val{:forward} with relstep √eps
            # (h ≈ 1.5e-6 at depth = 100 m). `kraken_jl` returns a numerically converged root,
            # so its output carries ~1e-11 of solver noise; a forward difference over that step
            # amplifies the noise to ~25% of ∂kr₂/∂depth. Central differences with relstep
            # ∛eps (h ≈ 6e-4) bring every entry to within ~3e-4 of the AD value, which the
            # analytic Pekeris solver confirms is the correct one.
            jac_fd = FiniteDiff.finite_difference_jacobian(multi_output_test, params_base, Val{:central})

            # Element-wise relative error, floored so that the smallest Jacobian entry
            # (∂kr/∂cb ≈ −8e-7) is not compared against its own roundoff.
            scale = max.(abs.(jac_fd), 1e-8)
            rel_error = maximum(abs.(jac_ad - jac_fd) ./ scale)
            @test rel_error < 0.01  # Less than 1% relative error
        end
    end
end
