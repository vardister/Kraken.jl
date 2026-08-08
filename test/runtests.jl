using Test
using TestItems
using TestItemRunner

@testset "Kraken.jl" begin
    # Run TestItems-based tests
    @run_package_tests

    # Include standard @testset-based tests
    include("numerical_methods_tests.jl")
    include("automatic_differentiation_tests.jl")

    # Performance tests (optional, can be slow)
    if get(ENV, "KRAKEN_RUN_PERFORMANCE_TESTS", "false") == "true"
        include("performance_tests.jl")
    end

    # Fortran cross-validation lands in Milestone 3 as test/fortran_reference_tests.jl,
    # driving unmodified kraken.exe from AcousticsToolbox_jll over .env/.mod files.
end
