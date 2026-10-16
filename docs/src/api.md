# API reference

```@meta
CurrentModule = Kraken
```

Everything Kraken.jl exports, generated from the docstrings in `src/`, so these pages cannot drift
from the code. Every example in them is a doctest: the docs build runs it and fails if the output
has changed.

The reference is split by what you are doing:

- [Environments](@ref api-environments): profiles, boundary conditions, attenuation units, and the
  standard test environments.
- [Solver](@ref api-solver): [`kraken_jl`](@ref), and the finite-difference pipeline underneath it,
  for driving the steps yourself.
- [Field and transmission loss](@ref api-field): pressure, transmission loss and broadband
  synthesis from a solved set of modes, plus plotting.
- [Pekeris closed form](@ref api-pekeris): the analytic two-layer model the solver is checked
  against.
- [Internals](@ref api-internals): documented but unexported helpers, for reading the source.

`makedocs` runs with `checkdocs=:exports`, so an exported name that none of these pages lists fails
the build instead of being left out of the docs.

```@index
Pages = ["api/environments.md", "api/solver.md", "api/field.md", "api/pekeris.md"]
```
