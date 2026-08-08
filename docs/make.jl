using Documenter
using Kraken

# Scaffolding only — Milestone 8 writes the actual guide content. What this build has to guarantee
# today is that it *runs*, that it sees the real module, and that every exported name has somewhere
# to be documented (checkdocs=:exports turns a missing docstring into a build failure).
makedocs(;
    sitename="Kraken.jl",
    modules=[Kraken],
    authors="Ariel Vardi",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://vardister.github.io/Kraken.jl"
    ),
    pages=["Home" => "index.md", "API reference" => "api.md"],
    # Every exported name currently has a docstring, so this is enforced rather than warned about:
    # adding an export without documenting it fails the build.
    checkdocs=:exports,
)
