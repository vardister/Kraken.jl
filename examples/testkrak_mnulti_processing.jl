"""
Same as testkrak1.jl but with multi-threading enabled. This is a test to see if the multi-threading is working correctly with Kraken and if it is faster than the single-threaded version. Also, how to use it with `ccall`?
"""

using Distributed # to get `pmap` and `@everywhere`
addprocs(4) # Add 4 processes
@everywhere using Kraken


using CairoMakie
using BenchmarkTools

@everywhere function main(freq)

nm = 41
hw = 71.0

ssp, b, sspHS = env_builder(hw=hw)

nsr = 1
zsr = 18.0
nrc = 70
zrc = [0.0, hw]

env_temp = EnvKRAKEN(
    ssp,
    b,
    sspHS,
    zrc,
    zsr,
)

return kraken(env_temp, freq; n_modes=nm)
end

freqs = fill(500.0, 100)
rest = Vector{Any}(undef, length(freqs))
res = Vector{Any}(undef, length(freqs))

@benchmark for (i, freq) in enumerate(freqs)
    rest[i] = main(freq)
end

@benchmark res .= pmap($main, $freqs)

# We can see an increase in speed! And a good one as well!
