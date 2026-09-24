# [Solver](@id api-solver)

```@meta
CurrentModule = Kraken
```

## Solving

```@docs
kraken_jl
```

## The finite-difference pipeline

What [`kraken_jl`](@ref) does on each mesh, in order. Call these directly only to work on a single
mesh, or to study one step.

```@docs
AcousticProblemProperties
AcousticProblemCache
finite_difference_coefficients
get_g
det_sturm
bisection
find_kr(::UnderwaterEnv, ::AcousticProblemProperties, ::AcousticProblemCache)
solve_for_kr
inverse_iteration
mode_eigenvector
normalize_mode
```

## Attenuation

```@docs
add_attenuation
modal_attenuation
complex_soundspeed
```
