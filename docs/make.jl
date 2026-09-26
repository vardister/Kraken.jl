using Documenter
using Kraken
# `import`, not `using`: `[density](@ref)` links resolve in Main, and Makie exports a `density` too.
using CairoMakie: CairoMakie

# Figures as PNG files next to each page. By default Makie also offers every figure as `text/html` (the
# same PNG, base64-inlined), which Documenter prefers, and which would put the images inside the HTML,
# over `example_size_threshold` and against each page's size limit. Disabling only `text/html` is not
# enough: CairoMakie then advertises `image/svg+xml` as well, Documenter prefers that too, and the TL
# heatmaps (~80k cells each) take minutes apiece to write as vector graphics.
CairoMakie.activate!(; type="png")
CairoMakie.disable_mime!("text/html", "image/svg+xml")

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
# The figures on those pages are drawn by Kraken's own plotting functions (ext/KrakenMakieExt.jl),
# so the docs environment depends on CairoMakie, and the build is also what exercises the extension:
# the test suite does not load Makie.
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

# Publishes to the `gh-pages` branch, from which GitHub Pages serves the site: `dev/` from every push
# to master, and a versioned copy (`v0.4.0/`, `stable/`) from every release tag. Anywhere else,
# including a local build and a push to any other branch, it logs why it is skipping and does nothing,
# so running this file by hand is still just a build. On GitHub Actions it authenticates with the
# workflow's own GITHUB_TOKEN (see .github/workflows/Docs.yml), so no deploy key is needed.
deploydocs(; repo="github.com/vardister/Kraken.jl.git", devbranch="master", push_preview=false)
