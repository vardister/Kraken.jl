using Documenter
using Kraken

# The build guarantees that the docs see the real module, that every exported name is documented
# (checkdocs=:exports turns a missing docstring into a build failure), and that the code on the pages
# still runs: `guide.md`, `environments.md`, `ad.md` and `validation.md` are `@example` blocks
# throughout, executed here, with their output rendered as it came out.
#
# `ad.md`'s blocks are the reason that matters most: the docs build is what verifies that
# the gradient-based inversion actually converges. That is deliberate — an inversion example that
# has silently stopped converging is worse than no example. It is also why the docs environment
# depends on Zygote and ForwardDiff.
#
# `validation.md` executes too, and is the most expensive page (~1 min): it runs Fortran KRAKEN and
# `field.exe` through the test suite's `test/reference/` harness and tabulates the agreement, so the
# numbers on it are measured by this build rather than copied in. That is why the docs environment
# depends on AcousticsToolbox_jll, LinearAlgebra and Markdown.
#
# Every docstring example is a `jldoctest`, so the build also checks that the API reference's examples
# still run and still print what they say. They share this one setup.
DocMeta.setdocmeta!(Kraken, :DocTestSetup, :(using Kraken); recursive=true)

makedocs(;
    sitename="Kraken.jl",
    modules=[Kraken],
    authors="Ariel Vardi",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vardister.github.io/Kraken.jl",
        # Named explicitly: Documenter otherwise asks `git remote` for the default branch, which fails
        # (with a warning) wherever the checkout has no reachable remote.
        edit_link="master",
    ),
    pages=[
        "Home" => "index.md",
        "Getting started" => "guide.md",
        "Environments" => "environments.md",
        "Automatic differentiation" => "ad.md",
        "Validation" => "validation.md",
        "API reference" =>
            ["api.md", "api/environments.md", "api/solver.md", "api/field.md", "api/pekeris.md", "api/internals.md"],
    ],
    # Every exported name currently has a docstring, so this is enforced rather than warned about:
    # adding an export without documenting it fails the build.
    checkdocs=:exports,
)
