# Kraken.jl

[![CI](https://github.com/vardister/Kraken.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/vardister/Kraken.jl/actions/workflows/CI.yml)
[![Format](https://github.com/vardister/Kraken.jl/actions/workflows/Format.yml/badge.svg)](https://github.com/vardister/Kraken.jl/actions/workflows/Format.yml)

**❗Documentation is currently under construction.**

KRAKEN.jl is a Normal-Mode based simulation package for underwater acoustic propagation. It's heavily based on [Michael Porter's KRAKEN Fortran code located in the Acoustics Toolbox](https://oalib-acoustics.org/models-and-software/normal-modes/) and [UnderwaterAcoustics.jl](https://github.com/org-arl/UnderwaterAcoustics.jl).

This reimplementation is fully written in Julia, and is designed to be more user-friendly, and easier to extend. It is also designed to be more efficient, and to take advantage of Julia's parallelization capabilities.

## Features

- Normal-Mode based simulation for underwater acoustic propagation fully written in Julia
- User-friendly and easy to extend
- Differentiable code (❗ currently only using [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl))
- Re-using existing environmental data files from the Acoustics Toolbox
- Optional plotting via a package extension — `using CairoMakie` (or `GLMakie`) enables `plot_modes` and `plot_ssp` at no cost to anyone who does not

## Missing features
- [ ] Compressional wave attenuation in environment
- [ ] Inclusion of shear wave properties in environment
- [ ] Cross-validation against the unmodified Fortran `kraken.exe` (in progress — it will be a test-only harness driving `AcousticsToolbox_jll`, not a shared library shipped with this package)

## Installation

```julia
using Pkg
Pkg.add("Kraken")
```

## Usage

```julia
using Kraken

# Load the environment
ssp, layers, sspHS = pekeris_env() # Similar structure to environment files from the Acoustics Toolbox
env = UnderwaterEnv(ssp, layers, sspHS)

# Run the simulation
freq = 100.0
sol = kraken_jl(env, freq)

# Access the results
modes = sol.modes
wavenumbers = sol.kr
zn = vcat(sol.props.zn_vec...)
```

### Calculating group speeds
Group speeds are defined as the derivative of the angular frequency $\omega = 2\pi f$ with respect to the wavenumbers $k_{r,m}$. This is written as

$$ c_g = \frac{\partial \omega}{\partial k_{r,m}} $$

As such, to calculate the group speeds using Kraken.jl we make use of automatic differentiation capabilities using
_ForwardDiff.jl_ and differentiate directly.

```julia
using ForwardDiff
using Kraken
using Roots


function calculate_kr_pekeris(freq)
    ssp, layers, sspHS = pekeris_env()
    env = UnderwaterEnv(ssp, layers, sspHS)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    wavenumbers = find_kr(env, props, cache; method=Roots.A42())
    return wavenumbers
end

freq = 100.0
group_speeds = 2pi ./ ForwardDiff.derivative(calculate_kr_pekeris, freq)
```


## More examples

More examples can be accessed in the `examples` folder. ❗Most of them do not currently run: they
call an `EnvKRAKEN` API that was removed when the Fortran sources moved out of this repository.
They are being rewritten against the current API.

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes. For major changes, please open an issue first to discuss what you would like to change.
