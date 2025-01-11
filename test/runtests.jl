using KRAKEN
using Test


@testset "kraken_pekeris.jl" begin
    env = PekerisUnderwaterEnv(1500.0, 1600.0, 1000.0, 1500.0, 100.0)
    p = [env.c1, env.cb, env.ρ1, env.ρb, env.depth]
    krs = find_kr(env, 100.0, p)
    krs = round.(krs, digits=6)
    # test the first 5 values of krs for the known pekeris waveguide
    @test all((krs .- [0.417908, 0.414964, 0.409971, 0.40286, 0.39383]) .< 1e4)
end


@testset "kraken_core.jl" begin
    ssp, layers, sspHS = pekeris_test_dict_KRAKEN()
    env = UnderwaterEnv(ssp, layers, sspHS)
    props = AcousticProblemProperties(env, 100.0)
    cache = AcousticProblemCache(env, props)
    krs = find_kr(env, props, cache)
    krs = round.(krs, digits=6)
    # test the first 5 values of krs for the known pekeris waveguide
    @test all((krs .- [0.417908, 0.414964, 0.409971, 0.40286, 0.39383]) .< 1e4)
end