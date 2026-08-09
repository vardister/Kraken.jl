using Documenter
using Kraken

# Mostly scaffolding — Milestone 8 writes the remaining guide content. What this build has to
# guarantee today is that it *runs*, that it sees the real module, and that every exported name has
# somewhere to be documented (checkdocs=:exports turns a missing docstring into a build failure).
#
# `ad.md` is the exception: its `@example` blocks execute, so the docs build is what verifies that
# the gradient-based inversion actually converges. That is deliberate — an inversion example that
# has silently stopped converging is worse than no example. It is also why the docs environment
# depends on Zygote and ForwardDiff.
makedocs(;
    sitename="Kraken.jl",
    modules=[Kraken],
    authors="Ariel Vardi",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://vardister.github.io/Kraken.jl"
    ),
    pages=["Home" => "index.md", "Automatic differentiation" => "ad.md", "API reference" => "api.md"],
    # Every exported name currently has a docstring, so this is enforced rather than warned about:
    # adding an export without documenting it fails the build.
    checkdocs=:exports,
)
