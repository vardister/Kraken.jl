# API reference

```@meta
CurrentModule = Kraken
```

Everything Kraken.jl exports. Generated from the docstrings in `src/`, so this page cannot drift
from the code.

```@index
```

`:constant` is in the list because Milestone 5 exports two — `ATTENUATION_UNIT_CHARS` and
`DEFAULT_ATTENUATION_UNITS`. `makedocs` runs with `checkdocs=:exports`, so an exported name with
nowhere to land here fails the build rather than quietly vanishing from the docs.

```@autodocs
Modules = [Kraken]
Order = [:type, :function, :constant]
```
