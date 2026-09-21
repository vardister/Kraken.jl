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

# Where the seventeen minutes go. Measured on a clean full run (1036.6 s end to end), each file timed
# on its own:
#
#   reverse_ad_tests.jl                 ~700 s   68%   Mooncake and Zygote compiling rules
#   fortran_reference_tests.jl           178 s   17%   one kraken.exe subprocess per case × frequency
#   automatic_differentiation_tests.jl   104 s   10%
#   environment_tests.jl                  27 s
#   integration_tests.jl                  13 s
#   numerical_methods_tests.jl            10 s
#
# The reverse-mode AD file is two thirds of it, and almost none of that is Kraken computing anything
# — it is the AD backends compiling. So the Fortran cross-validation, the obvious suspect, is not
# actually what makes the suite slow.
#
# `KRAKEN_SKIP_TESTS` is a comma-separated list of file names to leave out while iterating on code
# they cannot reach. Skipping `reverse_ad_tests.jl` alone takes the suite from ~17 min to ~5.5 min.
# It defaults to empty, CI never sets it, and every skip is announced loudly — a green run that
# quietly skipped two thirds of the suite would be worse than a slow one.
#
#     KRAKEN_SKIP_TESTS=reverse_ad_tests.jl julia --project=. -e 'using Pkg; Pkg.test()'
#
# This is for the edit/test loop only. Run it unset before pushing.
const SKIPPED_SUITES = Set(strip.(split(get(ENV, "KRAKEN_SKIP_TESTS", ""), ','; keepempty=false)))

let unknown = setdiff(
        SKIPPED_SUITES,
        Set([
            "numerical_methods_tests.jl",
            "automatic_differentiation_tests.jl",
            "reverse_ad_tests.jl",
            "fortran_reference_tests.jl",
        ]),
    )
    isempty(unknown) ||
        error("KRAKEN_SKIP_TESTS names files that are not includable suites: $(join(sort(collect(unknown)), ", "))")
end

"""
    include_suite(file)

`include` a test file unless `KRAKEN_SKIP_TESTS` asked for it to be left out, in which case say so.
"""
function include_suite(file)
    if file in SKIPPED_SUITES
        @warn "SKIPPING $file — KRAKEN_SKIP_TESTS is set. This run does NOT show the suite is green."
        return nothing
    end
    return include(file)
end

@testset "Kraken.jl" begin
    # Run TestItems-based tests. These are the cheap ones (~40 s together) and always run.
    @run_package_tests

    # Include standard @testset-based tests
    include_suite("numerical_methods_tests.jl")
    include_suite("automatic_differentiation_tests.jl")
    include_suite("reverse_ad_tests.jl")

    # Cross-validation against unmodified kraken.exe from AcousticsToolbox_jll. Self-skipping when
    # no binary resolves, so this is safe to run unconditionally on any platform.
    include_suite("fortran_reference_tests.jl")

    # Performance tests (optional, can be slow)
    if get(ENV, "KRAKEN_RUN_PERFORMANCE_TESTS", "false") == "true"
        include("performance_tests.jl")
    end
end
