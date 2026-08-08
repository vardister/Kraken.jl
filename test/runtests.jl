using Test
using TestItems
using TestItemRunner

@testset "Kraken.jl" begin
    # Run TestItems-based tests
    @run_package_tests

    # Include standard @testset-based tests
    include("numerical_methods_tests.jl")
    include("automatic_differentiation_tests.jl")

    # Cross-validation against unmodified kraken.exe from AcousticsToolbox_jll. Self-skipping when
    # no binary resolves, so this is safe to run unconditionally on any platform.
    include("fortran_reference_tests.jl")

    # Performance tests (optional, can be slow)
    if get(ENV, "KRAKEN_RUN_PERFORMANCE_TESTS", "false") == "true"
        include("performance_tests.jl")
    end
end
