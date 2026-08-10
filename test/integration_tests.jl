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

@testitem "M5.2: a lossless environment is untouched by the attenuation path" begin
    using Kraken

    env = UnderwaterEnv(pekeris_env()...)
    sol = kraken_jl(env, 100.0)

    # No attenuation anywhere means real wavenumbers, of the same concrete element type as before
    # attenuation existed -- not `ComplexF64` with a zero imaginary part.
    @test !is_lossy(env)
    @test eltype(sol.kr) === Float64
    @test all(isreal, sol.kr)

    # An environment that spells its zero attenuation out explicitly takes the same branch and must
    # produce the *identical* vector, bit for bit -- this is what the `is_lossy` guard buys.
    ssp, layers, sspHS = pekeris_env()
    ssp[:, 5] .= 0.0
    sspHS[2, 5] = 0.0
    @test kraken_jl(UnderwaterEnv(ssp, layers, sspHS), 100.0).kr == sol.kr

    # Declaring different units cannot change a lossless answer either.
    for u in (:nepers_per_m, :dB_per_m, :dB_per_kmHz, :Q, :loss_parameter)
        @test kraken_jl(UnderwaterEnv(ssp, layers, sspHS; atten_units=u), 100.0).kr == sol.kr
    end

    # The single-mesh path goes through `add_attenuation` too.
    @test eltype(kraken_jl(env, 100.0; n_meshes=1).kr) === Float64
end

@testitem "M5.2: half-space attenuation gives decaying modes that scale with αb" begin
    using Kraken

    freq = 100.0
    function with_αb(αb)
        ssp, layers, sspHS = pekeris_env()
        sspHS[2, 5] = αb
        return UnderwaterEnv(ssp, layers, sspHS)
    end

    base = kraken_jl(with_αb(0.0), freq)
    sols = [kraken_jl(with_αb(a), freq) for a in (0.25, 0.5, 1.0)]

    for s in sols
        @test eltype(s.kr) === ComplexF64
        @test length(s.kr) == length(base.kr)
        # Negative imaginary part: with KRAKEN's sign convention that is a mode losing energy with
        # range. A positive one would be a mode that grows, which is the classic sign-error symptom.
        @test all(imag.(s.kr) .< 0)
        # The real part is a second-order effect of a first-order imaginary perturbation, so it
        # barely moves.
        @test real.(s.kr) ≈ base.kr rtol = 1e-6
    end

    # Perturbation theory is linear in the attenuation -- but only where the perturbation is small
    # compared with what it perturbs, and for the half-space term that comparison is against γ²
    # rather than against kᵣ². The half-space contribution goes as Im√(γ² + 2iω a_b/c_b), and γ → 0
    # at the bottom cutoff, so the *least* trapped mode is the one where linearity fails first. For
    # this waveguide at 100 Hz, 2ω a_b/c_b is 28-70% of γ² for modes 1-4 and 6.4× γ² for mode 5.
    #
    # That is a real property of the method (and of `kraken.exe`, which evaluates the same
    # expression), not an artifact: it is exactly the "first order, and it degrades" caveat that
    # makes `krakenc.exe` the right tool for strongly attenuating bottoms.
    # Concretely, three things hold, and together they say "linear where the theory applies,
    # saturating where it does not":
    #
    #  1. every mode is *sub*linear -- the square root grows more slowly than its tangent, so
    #     doubling αb always gives strictly less than twice the loss;
    #  2. the best-trapped modes are linear to well under a percent;
    #  3. the departure grows monotonically with mode number, i.e. with proximity to cutoff.
    ratio2 = abs.(imag.(sols[2].kr)) ./ (2 .* abs.(imag.(sols[1].kr)))
    ratio4 = abs.(imag.(sols[3].kr)) ./ (4 .* abs.(imag.(sols[1].kr)))

    @test all(ratio2 .< 1)
    @test all(ratio4 .< 1)
    @test all(ratio4 .< ratio2)               # more attenuation, further from linear
    @test imag.(sols[2].kr[1:3]) ≈ 2 .* imag.(sols[1].kr[1:3]) rtol = 1e-2
    @test issorted(ratio2; rev=true)          # worse the closer to cutoff
    @test issorted(ratio4; rev=true)
    @test ratio2[1] > 0.99                    # best-trapped mode: linear to 1%
    @test ratio2[end] < 0.9                   # near-cutoff mode: visibly saturated

    # Monotone mode by mode, not just on average.
    @test all(imag.(sols[2].kr) .< imag.(sols[1].kr))
    @test all(imag.(sols[3].kr) .< imag.(sols[2].kr))

    # Higher modes strike the bottom at steeper angles and so spend more of their energy in the
    # lossy half-space: the loss increases with mode number.
    @test issorted(imag.(sols[3].kr); rev=true)
end

@testitem "M5.2: the perturbation integral matches its closed form" begin
    using Kraken

    # With loss in the water column only, and c and ρ uniform there, the perturbation integral
    # collapses to something computable by hand:
    #
    #     Im(kᵣ) = -ω a (1 - E_hs) / (c kᵣ),    E_hs = ψ(D)² / (2 ρ_b γ)
    #
    # where `a` is the attenuation in nepers/m and `E_hs` the fraction of the mode's energy sitting
    # in the half-space (the normalization makes the two fractions sum to one). Nothing below reuses
    # `modal_attenuation`, so this is a check on the formula and not on its own arithmetic.
    αw, c0, freq = 0.3, 1500.0, 100.0
    ω = 2π * freq
    ssp, layers, sspHS = pekeris_env(; c0=c0)
    ssp[:, 5] .= αw
    env = UnderwaterEnv(ssp, layers, sspHS)

    sol = kraken_jl(env, freq)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    kr_coarse, ψ = inverse_iteration(find_kr(env, props, cache), env, props, cache)

    a = attenuation_nepers_per_m(αw, c0, freq, :dB_per_wavelength)
    predicted = map(eachindex(kr_coarse)) do m
        kr = kr_coarse[m]
        γ = sqrt(kr^2 - (ω / env.cb)^2)
        E_hs = ψ[end, m]^2 / (2 * env.ρb * γ)
        return -ω * a * (1 - E_hs) / (c0 * kr)
    end

    @test imag.(sol.kr) ≈ predicted rtol = 1e-3

    # A mode travelling almost horizontally in a uniform lossy medium loses very nearly the medium's
    # own attenuation per metre -- the classic sanity check on the units.
    @test imag(sol.kr[1]) ≈ -a rtol = 1e-2

    # Loss is loss however the units spell it: the same physical attenuation declared in nepers/m
    # must give the same answer as the dB/wavelength version above.
    ssp_n, layers_n, sspHS_n = pekeris_env(; c0=c0)
    ssp_n[:, 5] .= a
    same = kraken_jl(UnderwaterEnv(ssp_n, layers_n, sspHS_n; atten_units=:nepers_per_m), freq)
    @test same.kr ≈ sol.kr rtol = 1e-12
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
