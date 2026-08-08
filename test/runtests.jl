using Test
using TestItems
using TestItemRunner
using Kraken

"""
    check_testing_this_checkout()

Fail loudly if `Kraken` resolved to some *other* copy of the package than the one this `test/`
directory belongs to.

`test/Manifest.toml` is gitignored, so a fresh clone — **or a new git worktree**, which is the easy
one to forget — starts without one. Instantiating the test environment then resolves `Kraken` from
the General registry, where it is a *registered, released* package: the suite runs happily against
code nobody in this working tree wrote. It surfaces as a wall of `UndefVarError` on every symbol
added since that release, which points nowhere near the actual cause.

`Pkg.test()` from the root environment is immune (it always uses the local package), so the two run
modes disagree — and the immune one is not the one to trust when they do.
"""
function check_testing_this_checkout()
    repo = realpath(joinpath(@__DIR__, ".."))
    resolved = realpath(pkgdir(Kraken))
    resolved == repo && return nothing
    return error("""
    Kraken resolved to a different copy of the package than the one under test.

        this checkout : $repo
        resolved to   : $resolved

    The test environment has no dev-link to this working tree, so it took `Kraken` from the
    registry. Every test below would run against the released version. Fix it with:

        julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

    run from $repo — once per clone and once per git worktree. Verify by checking that the
    [[deps.Kraken]] entry in test/Manifest.toml says `path = ".."`.
    """)
end

check_testing_this_checkout()

@testset "Kraken.jl" begin
    # Run TestItems-based tests
    @run_package_tests

    # Include standard @testset-based tests
    include("numerical_methods_tests.jl")
    include("automatic_differentiation_tests.jl")
    include("reverse_ad_tests.jl")

    # Cross-validation against unmodified kraken.exe from AcousticsToolbox_jll. Self-skipping when
    # no binary resolves, so this is safe to run unconditionally on any platform.
    include("fortran_reference_tests.jl")

    # Performance tests (optional, can be slow)
    if get(ENV, "KRAKEN_RUN_PERFORMANCE_TESTS", "false") == "true"
        include("performance_tests.jl")
    end
end
