module Kraken

### New Julia implementation
include("kraken_core.jl")
include("kraken_basic.jl")
include("kraken_pekeris.jl")

include("kraken_standard_environments.jl")

### Legacy Fortran implementation
include("kraken_fortran.jl")
export kraken
export EnvKRAKEN
export env_builder
export pf_adiabatic
export pf_adiabatic_signal
export pf_signal
export pressure_field_fortran

end