# [Environments](@id api-environments)

```@meta
CurrentModule = Kraken
```

## The environment

```@docs
UnderwaterEnv
UnderwaterEnvFORTRAN
is_lossy
```

## Profiles

Depths are in metres, positive downwards from the surface.

```@docs
SampledSSP
SampledDensity
SampledAttenuation
soundspeed
maxsoundspeed
density
attenuation
SSP_INTERPOLATION_CHARS
DEFAULT_SSP_INTERPOLATION
```

## Boundary conditions

```@docs
BoundaryCondition
PressureRelease
RigidBoundary
AcousticHalfspace
```

## Attenuation units

```@docs
attenuation_nepers_per_m
ATTENUATION_UNIT_CHARS
DEFAULT_ATTENUATION_UNITS
```

## Standard environments

Each returns the `(ssp, layers, sspHS)` triple that [`UnderwaterEnv`](@ref) takes.

```@autodocs
Modules = [Kraken]
Pages = ["kraken_standard_environments.jl"]
Private = false
```
