using DrWatson

using NamedArrays
# include("src_normalmodes.jl")

### Standard Pekeris
# function pekeris_test_env()

#     c = 1500.0
#     cb = 1600.0
#     ρ = 1000.0
#     ρb = 1500.0
#     depth = 100.0
#     freq = 100.0

#     return PekerisEnv(freq, c, cb, ρ, ρb, depth)
# end

function pekeris_test_dict()
    # Water column
    c = 1500.0
    ρ = 1000.0
    # Bottom half-space
    cb = 1600.0
    ρb = 1500.0
    # other
    depth = 100.0
    freq = 100.0

    return Dict(
    "c" => c, "cb" => cb, "ρ" => ρ, "ρb" => ρb, "depth" => depth, "freq" => freq)
end


### Standard 1-layer sediment model with constant sound speeds
function one_layer_test_dict()
    # Water column
    c0 = 1500.0
    ρ0 = 1000.0
    h0 = 100.0
    # sediment layer
    c1 = 1550.0
    ρ1 = 15002.0
    h1 = 20.0
    # bottom half-space
    cb = 1600.0
    ρb = 2000.0
    # other
    depth = h0 + h1;
    freq = 100.0

    return @dict c0 c1 cb ρ0 ρ1 ρb h0 h1 depth freq
end


function one_layer_test_dict_KRAKEN()
    # Water column
    c0 = 1500.0
    ρ0 = 1000.0
    h0 = 100.0
    α0 = 0.0
    # sediment layer
    c1 = 1550.0
    ρ1 = 1500.0
    h1 = 20.0
    α1 = 0.0
    # bottom half-space
    cb = 1600.0
    ρb = 2000.0
    αb = 0.0
    # other
    freq = 100.0
    z0 = h0
    z1 = h0 + h1

    ssp = NamedArray(
            [0.0 c0 ρ0 α0
             z0  c0 ρ0 α0
             z0+eps(z0)  c1 ρ1 α1
             z1  c1 ρ1 α1], dimnames = ("i", "Paramater"))
    setnames!(ssp, ["z", "c", "ρ", "α"], 2)
    ssp = [0.0 c0 ρ0 α0
         z0  c0 ρ0 α0
         z0+eps(z0)  c1 ρ1 α1
         z1  c1 ρ1 α1]

    layers = [0.0 0.0 z0
              0.0 0.0 z1]

    sspHS = [0.0 343.0 0.00121 0.0
             z1  cb    ρb      αb]



    env_dict = Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
    return env_dict
end

