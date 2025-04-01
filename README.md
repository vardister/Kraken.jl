# Kraken.jl

**❗Documentation is currently under construction.**

KRAKEN.jl is a Normal-Mode based simulation package for underwater acoustic propagation. It's heavily based on [Michael Porter's KRAKEN Fortran code located in the Acoustics Toolbox](https://oalib-acoustics.org/models-and-software/normal-modes/) and [UnderwaterAcoustics.jl](https://github.com/org-arl/UnderwaterAcoustics.jl).

This reimplementation is fully written in Julia, and is designed to be more user-friendly, and easier to extend. It is also designed to be more efficient, and to take advantage of Julia's parallelization capabilities.

Access to the Fortran code is also available through Julia calls to the shared library, which can be useful for comparison and validation purposes.


## Features

- Normal-Mode based simulation for underwater acoustic propagation fully written in Julia
- User-friendly and easy to extend
- Differential code (❗ currently only using [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl))
- Access to the Fortran code through Julia calls to the shared library
- Re-using existing environmental data files from the Acoustics Toolbox

## Missing features
- [ ] Compressional wave attenuation in environment
- [ ] Inclusion of shear wave properties in environment

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
Group speeds are defined as the derivative of the wavenumbers $k_{r,m}$ with respect to the angular frequency $\omega$.
As such, to calculate the group speeds using Kraken.jl we make use of automatic differentiation capabilities using
_ForwardDiff.jl_ and differentiate directly.

```julia
using ForwardDiff
using Kraken
using Roots

# Load the environment
ssp, layers, sspHS = pekeris_env() # Similar structure to environment files from the Acoustics Toolbox
env = UnderwaterEnv(ssp, layers, sspHS)

function calculate_kr_pekeris(freq)
    ssp, layers, sspHS = pekeris_env() # Similar structure to environment files from the Acoustics Toolbox
    env = UnderwaterEnv(ssp, layers, sspHS)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    wavenumbers = find_kr(env, props, cache; method=Roots.A42())
    return wavenumbers
end

freq = 100.0
group_speeds = ForwardDiff.derivative(calculate_kr_pekeris, freq)
```


## More examples
More examples can be accessed in the `examples` folder.

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes. For major changes, please open an issue first to discuss what you would like to change.
