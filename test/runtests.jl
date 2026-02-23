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

    # Fortran interface tests require KrakenFortran.jl package
    # Skip if not available
    if get(ENV, "KRAKEN_RUN_FORTRAN_TESTS", "false") == "true"
        include("fortran_interface_tests.jl")
    end
end
