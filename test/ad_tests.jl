using Test
using Kraken

using ForwardDiff
using Enzyme
using SciMLSensitivity

# using Infiltrator
# using Cthulhu

# test AD on the bisection method
function test_ad_bisection(c0=1500.0 <: Real, cb=1600.0 <: Real)
    ssp, layers, sspHS = pekeris_test_dict_KRAKEN(c0=c0, cb=cb)
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    # @infiltrate
    return bisection(env, props, cache)
end

function test_ad_bisection_inp(x)
    ssp, layers, sspHS = pekeris_test_dict_KRAKEN(c0=x[1], cb=x[2])
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    # @infiltrate
    return bisection(env, props, cache)
end

function test_ad_bisection_single(x)
    ssp, layers, sspHS = pekeris_test_dict_KRAKEN(c0=x)
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    # @infiltrate
    return bisection(env, props, cache)
end

@testset "AD with Bisection function" begin
    # Currently testing only Pekeris waveguide ONLY
    # quick check that the function works
    ret = test_ad_bisection(1500, 1600)
    # diff with ForwardDiff
    # Using a single parameter input (derivative)
    @test begin
        try
            # Diff using Enzyme
            # Forward AD
            println("Forward jacobian with Enzyme")
            jacobian(Forward, test_ad_bisection_inp, [1500.0, 1600.0])
            println("Forward Derivative with Enzyme")
            gradient(Forward, test_ad_bisection_single, 1500.0)
            # Using an input vector (jacobian)
            println("Jacobian with ForwardDiff")
            ForwardDiff.jacobian(test_ad_bisection_inp, [1500.0, 1600.0])
            println("Derivative with ForwardDiff")
            ForwardDiff.derivative(test_ad_bisection_single, 1500.0)
            true
        catch
            false
        end
    end


    @test begin
        try
            # Reverse
            println("Reverse jacobiant with Enzyme")
            jacobian(Reverse, test_ad_bisection_inp, [1500.0, 1600.0])
            true
        catch
            false
        end
    end
end


### Test AD on `find_kr` function (the hardest one)
function test_ad_find_kr_inp(x)
    ssp, layers, sspHS = pekeris_test_dict_KRAKEN(c0=x[1], cb=x[2])
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    return find_kr(env, props, cache)
end

function test_ad_find_kr_inp_single(x)
    ssp, layers, sspHS = pekeris_test_dict_KRAKEN(c0=x)
    env = UnderwaterEnv(ssp, layers, sspHS)
    freq = 100.0
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    return find_kr(env, props, cache)
end


ret = test_ad_find_kr_inp([1500, 1600])
ForwardDiff.jacobian(test_ad_find_kr_inp, [1500.0, 1600.0])
jacobian(Forward, test_ad_find_kr_inp, [1500.0, 1600.0])
gradient(Forward, test_ad_find_kr_inp_single, 1500.0)
# jacobian(Reverse, test_ad_find_kr_inp, [1500.0, 1600.0])




# @testset "AD with Bisection function" begin
#     # Currently testing only Pekeris waveguide ONLY
#     # quick check that the function works
#     ret = test_ad_bisection(1500, 1600)
#     # diff with ForwardDiff
#     # Using a single parameter input (derivative)
#     println("Derivative with ForwardDiff")
#     @test begin
#         try
#             # Diff using Enzyme
#             # Forward AD
#             println("Forward jacobian with Enzyme")
#             jacobian(Forward, test_ad_find_kr_inp, [1500.0, 1600.0])
#             println("Forward Derivative with Enzyme")
#             gradient(Forward, test_ad_find_kr_single, 1500.0)
#             # Using an input vector (jacobian)
#             println("Jacobian with ForwardDiff")
#             ForwardDiff.jacobian(test_ad_find_kr_inp, [1500.0, 1600.0])
#             true
#         catch
#             false
#         end
#     end


#     @test begin
#         try
#             # Reverse
#             println("Reverse jacobiant with Enzyme")
#             jacobian(Reverse, test_ad_find_kr_inp, [1500.0, 1600.0])
#             true
#         catch
#             false
#         end
#     end
# end