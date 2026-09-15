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

@testitem "SampledAttenuation1D Construction and Evaluation" begin
    using Kraken

    depths = [0.0, 50.0, 100.0]
    alphas = [0.0, 0.2, 0.6]

    α = SampledAttenuation(depths, alphas)

    @test attenuation(α, 0.0) ≈ 0.0
    @test attenuation(α, 50.0) ≈ 0.2
    @test attenuation(α, 100.0) ≈ 0.6
    @test attenuation(α, 25.0) ≈ 0.1     # linear, like the other two profiles
    @test attenuation(α, 75.0) ≈ 0.4

    # Constant extrapolation past the ends, matching SampledSSP/SampledDensity.
    @test attenuation(α, 200.0) ≈ 0.6
    @test !occursin("type", sprint(show, α))    # `show` reaches its own fields, not a missing one
end

@testitem "M5.1: attenuation units convert to nepers/m as CRCI does" begin
    using Kraken

    # Every expected value below is the arithmetic written out in `CRCI` (Acoustics Toolbox
    # `misc/AttenMod.f90`), transcribed independently here rather than calling the implementation --
    # a test that reuses the function's own constants only checks that it is self-consistent.
    c = 1500.0
    freq = 250.0
    ω = 2π * freq
    α = 0.5

    @test attenuation_nepers_per_m(α, c, freq, :nepers_per_m) == α
    @test attenuation_nepers_per_m(α, c, freq, :dB_per_m) == α / 8.6858896
    @test attenuation_nepers_per_m(α, c, freq, :dB_per_kmHz) == α * freq / 8685.8896
    @test attenuation_nepers_per_m(α, c, freq, :dB_per_wavelength) == α * freq / (8.6858896 * c)
    @test attenuation_nepers_per_m(α, c, freq, :Q) == ω / (2 * c * α)
    @test attenuation_nepers_per_m(α, c, freq, :loss_parameter) == α * ω / c

    # Sanity in physical terms: 1 dB/wavelength at 1500 m/s and 250 Hz is a 6 m wavelength, so the
    # loss is 1 dB per 6 m = 1/(6 * 8.6859) nepers/m.
    @test attenuation_nepers_per_m(1.0, c, freq, :dB_per_wavelength) ≈ 1 / ((c / freq) * 8.6858896)

    # The frequency-independent conventions really are frequency independent, and the other four
    # really are not.
    for u in (:nepers_per_m, :dB_per_m)
        @test attenuation_nepers_per_m(α, c, 10.0, u) == attenuation_nepers_per_m(α, c, 1000.0, u)
    end
    for u in (:dB_per_kmHz, :dB_per_wavelength, :loss_parameter)
        @test attenuation_nepers_per_m(α, c, 1000.0, u) ≈ 100 * attenuation_nepers_per_m(α, c, 10.0, u)
    end
    # Q is a *quality* factor: bigger Q means less loss, and loss still scales with frequency.
    @test attenuation_nepers_per_m(1000.0, c, freq, :Q) < attenuation_nepers_per_m(10.0, c, freq, :Q)

    # KRAKEN's own degenerate guards: these mean "lossless", not "divide by zero".
    @test attenuation_nepers_per_m(0.0, c, freq, :Q) == 0.0
    @test attenuation_nepers_per_m(α, 0.0, freq, :dB_per_wavelength) == 0.0
    @test attenuation_nepers_per_m(α, 0.0, freq, :loss_parameter) == 0.0
    @test all(isfinite, [attenuation_nepers_per_m(0.0, c, freq, u) for u in values(ATTENUATION_UNIT_CHARS)])

    # A lossless environment is lossless under every convention.
    for u in values(ATTENUATION_UNIT_CHARS)
        @test attenuation_nepers_per_m(0.0, c, freq, u) == 0.0
    end

    @test_throws ArgumentError attenuation_nepers_per_m(α, c, freq, :furlongs_per_fortnight)

    # The character table is the one `ReadTopOpt` accepts, and covers all six conventions.
    @test ATTENUATION_UNIT_CHARS['N'] === :nepers_per_m
    @test ATTENUATION_UNIT_CHARS['M'] === :dB_per_m          # dB per METRE -- not dB/km
    @test ATTENUATION_UNIT_CHARS['F'] === :dB_per_kmHz
    @test ATTENUATION_UNIT_CHARS['W'] === :dB_per_wavelength
    @test ATTENUATION_UNIT_CHARS['Q'] === :Q
    @test ATTENUATION_UNIT_CHARS['L'] === :loss_parameter
    @test length(ATTENUATION_UNIT_CHARS) == 6
    @test DEFAULT_ATTENUATION_UNITS === :dB_per_wavelength

    # dB/(m kHz) and dB/(km Hz) are the same quantity, which is why the one symbol serves both names.
    @test attenuation_nepers_per_m(1.0, c, 1000.0, :dB_per_kmHz) ≈ 1.0 / 8.6858896
end

@testitem "M5.1: environments carry their attenuation and know when they are lossy" begin
    using Kraken

    env = UnderwaterEnv(pekeris_env()...)
    @test !is_lossy(env)
    @test env.αb == 0.0
    @test all(iszero, env.α.α)
    @test env.atten_units === DEFAULT_ATTENUATION_UNITS

    # `pekeris_env` writes αb into sspHS[2, 5] and αp into ssp[:, 5]; both must reach the env.
    ssp, layers, sspHS = pekeris_env()
    sspHS[2, 5] = 0.5
    lossy = UnderwaterEnv(ssp, layers, sspHS)
    @test is_lossy(lossy)
    @test lossy.αb == 0.5
    @test occursin("lossy", sprint(show, lossy))

    ssp2, layers2, sspHS2 = pekeris_env()
    ssp2[:, 5] .= 0.01
    watery = UnderwaterEnv(ssp2, layers2, sspHS2)
    @test is_lossy(watery)
    @test attenuation(watery.α, 50.0) ≈ 0.01

    # The units are a construction-time choice and are carried, not guessed.
    nepers = UnderwaterEnv(ssp, layers, sspHS; atten_units=:nepers_per_m)
    @test nepers.atten_units === :nepers_per_m
    @test occursin("nepers_per_m", sprint(show, nepers))

    # The UnderwaterEnvFORTRAN path has to agree with the matrix path -- it is the same file layout.
    envf = UnderwaterEnv(UnderwaterEnvFORTRAN(ssp, layers, sspHS); atten_units=:nepers_per_m)
    @test envf.αb == lossy.αb
    @test envf.atten_units === :nepers_per_m
    @test envf.α.α == lossy.α.α
end

@testitem "M6.1: environments carry their boundary conditions" begin
    using Kraken

    ssp, layers, sspHS = pekeris_env()

    # The defaults are the one configuration the solver had before boundary conditions were types.
    env = UnderwaterEnv(ssp, layers, sspHS)
    @test env.top_bc === PressureRelease()
    @test env.bottom_bc === AcousticHalfspace()
    @test isconcretetype(typeof(env))
    shown = sprint(show, env)
    @test occursin("top: PressureRelease()", shown)
    @test occursin("bottom: AcousticHalfspace()", shown)

    # Spelling the defaults out is the same environment, down to the last bit of every wavenumber.
    explicit = UnderwaterEnv(ssp, layers, sspHS; top_bc=PressureRelease(), bottom_bc=AcousticHalfspace())
    @test typeof(explicit) === typeof(env)
    @test kraken_jl(explicit, 100.0).kr == kraken_jl(env, 100.0).kr

    # Both constructors carry a non-default choice through ...
    rigid = UnderwaterEnv(ssp, layers, sspHS; top_bc=RigidBoundary(), bottom_bc=PressureRelease())
    @test rigid.top_bc === RigidBoundary()
    @test rigid.bottom_bc === PressureRelease()
    @test occursin("top: RigidBoundary()", sprint(show, rigid))
    envf = UnderwaterEnv(
        UnderwaterEnvFORTRAN(ssp, layers, sspHS); top_bc=RigidBoundary(), bottom_bc=PressureRelease()
    )
    @test envf.top_bc === rigid.top_bc
    @test envf.bottom_bc === rigid.bottom_bc

    # ... and a top acoustic half-space, which is out of scope (see task 6.1's outcome in the plan),
    # constructs but is refused at solve time rather than quietly solved as something else.
    for bottom in (AcousticHalfspace(), RigidBoundary(), PressureRelease())
        halfspace_top = UnderwaterEnv(ssp, layers, sspHS; top_bc=AcousticHalfspace(), bottom_bc=bottom)
        @test_throws ArgumentError AcousticProblemProperties(halfspace_top, 100.0)
        @test_throws ArgumentError kraken_jl(halfspace_top, 100.0)
    end
end

@testitem "M6.2: rigid and vacuum boundaries reproduce the analytic isovelocity waveguide" begin
    using Kraken

    # An isovelocity column between two perfect boundaries has closed-form modes: ψ is a sine or cosine
    # of kz·z with kz fixed by the pair, kr = √(k² − kz²), and `∫ψ²/ρ dz = 1` makes the amplitude
    # √(2ρ/D). `pekeris_env`'s half-space row is simply not used by a rigid or vacuum bottom.
    c0, ρ0, D, freq = 1500.0, 1000.0, 100.0, 100.0
    k = 2π * freq / c0
    amplitude = sqrt(2ρ0 / D)
    ssp, layers, sspHS = pekeris_env(; c0=c0, ρ0=ρ0, depth=D)

    # (top, bottom) => (kz of mode m, normalized mode at depth z). Rigid over rigid starts with the plane
    # wave ψ = √(ρ/D), kz = 0: its kr is k *exactly*, on the upper bound of `bisection`'s search, which
    # is the case `kr_search_max` exists for.
    cases = [
        (PressureRelease(), RigidBoundary()) => (m -> (m - 0.5) * π / D, (kz, z) -> amplitude * sin(kz * z)),
        (PressureRelease(), PressureRelease()) => (m -> m * π / D, (kz, z) -> amplitude * sin(kz * z)),
        (RigidBoundary(), PressureRelease()) => (m -> (m - 0.5) * π / D, (kz, z) -> amplitude * cos(kz * z)),
        (RigidBoundary(), RigidBoundary()) =>
            (m -> (m - 1) * π / D, (kz, z) -> (iszero(kz) ? sqrt(ρ0 / D) : amplitude) * cos(kz * z)),
    ]
    for ((top, bottom), (kz_of, ψ_of)) in cases
        env = UnderwaterEnv(ssp, layers, sspHS; top_bc=top, bottom_bc=bottom)
        sol = kraken_jl(env, freq)

        # Every mode below cutoff is found, and nothing else — a perfect bottom traps them all.
        @test length(sol.kr) == count(m -> kz_of(m) < k, 1:100)
        kr_exact = sqrt.(k^2 .- kz_of.(eachindex(sol.kr)) .^ 2)
        @test maximum(abs.(sol.kr .- kr_exact) ./ kr_exact) < 1e-6

        # The mesh reaches exactly the boundaries that are unknowns, and no further.
        z = reduce(vcat, sol.props.zn_vec)
        @test (first(z) == 0) == (top isa RigidBoundary)
        @test (last(z) ≈ D) == !(bottom isa PressureRelease)

        # Mode shapes, including their normalization. On an isovelocity column the sampled sine or
        # cosine is an *exact* eigenvector of the finite-difference operator, so even the coarse-mesh
        # shapes agree to ~1e-9 (measured); 1e-6 leaves room without hiding a wrong end condition, which
        # would show up at O(1). The solver fixes the sign so `ψ[1] ≥ 0`, as both analytic shapes are.
        for m in 1:3
            @test maximum(abs.(sol.modes[:, m] .- ψ_of.(kz_of(m), z))) < 1e-6 * amplitude
        end
    end

    # A uniformly lossy column between perfect boundaries has no half-space term, so the perturbation of
    # every mode is Im(ω²/c̃²) times its normalization integral, which is one — provided the ends close
    # the attenuation integral exactly as they close the energy integral.
    α0 = 0.1  # dB per wavelength, the default units
    lossy_ssp, lossy_layers, lossy_sspHS = pekeris_env(; c0=c0, ρ0=ρ0, depth=D, α0=α0)
    ω = 2π * freq
    α_nepers = attenuation_nepers_per_m(α0, c0, freq, DEFAULT_ATTENUATION_UNITS)
    expected = imag(ω^2 / Kraken.complex_soundspeed(c0, α_nepers, ω)^2)
    for (top, bottom) in first.(cases)
        sol = kraken_jl(UnderwaterEnv(lossy_ssp, lossy_layers, lossy_sspHS; top_bc=top, bottom_bc=bottom), freq)
        @test all(m -> isapprox(imag(sol.kr[m]^2), expected; rtol=1e-8), eachindex(sol.kr))
    end
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
