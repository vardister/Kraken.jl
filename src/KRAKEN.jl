module KRAKEN

### New Julia implementation
include("kraken_core.jl")

include("kraken_pekeris.jl")
export PekerisUnderwaterEnv

### Legacy Fortran implementation
include("kraken_fortran.jl")
export kraken
export EnvKRAKEN
export env_builder
export pf_adiabatic
export pf_adiabatic_signal
export pf_signal

end