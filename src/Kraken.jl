module Kraken

### New Julia implementation
include("kraken_core.jl")
include("kraken_ad.jl")
include("kraken_pekeris.jl")
include("kraken_field.jl")

include("kraken_standard_environments.jl")

### Plotting
# Stubs and docstrings; the implementations are in ext/KrakenMakieExt.jl.
include("kraken_plotting.jl")

end
