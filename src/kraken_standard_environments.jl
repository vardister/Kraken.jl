using NamedArrays


# There are 4 overall cases to test
export pekeris_test_dict_KRAKEN, one_layer_test_dict_KRAKEN, one_layer_slope_test_dict_KRAKEN, two_layer_slope_test_dict_KRAKEN


### Standard Pekeris
function pekeris_test_dict_KRAKEN(;c0::Real=1500.0, cb::Real=1600.0, ρ0::Real=1000.0, ρb::Real=1500.0, h0::Real=100.0)
    # Input validation
    # for (param_name, param_value) in zip([:c0, :cb, :ρ0, :ρb, :h0], [c0, cb, ρ0, ρb, h0])
    #     if !isfinite(param_value)
    #         throw(DomainError(param_name, "Parameter must be a finite real number"))
    #     end
    #     if param_value <= 0
    #         throw(DomainError(param_name, "Parameter must be positive"))
    #     end
    # end

    # Water column
    α0 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    freq = 100.0
    z0 = h0

    # ssp = NamedArray(
    #     [0.0 c0 0.0 ρ0 α0 0.0
    #      h0 c0 0.0 ρ0 α0 0.0];
    #     dimnames = ("i", "Paramater")
    # )
    # setnames!(ssp, ["z", "cp", "cs", "ρ", "αp", "αs"], 2)
    T = promote_type(typeof(c0), typeof(cb), typeof(ρ0), typeof(ρb), typeof(h0))
    ssp = zeros(T, 2, 6)
    ssp[1, :] = [0.0 c0 0.0 ρ0 α0 0.0]
    ssp[2, :] = [h0 c0 0.0 ρ0 α0 0.0]

    # ssp = [0.0 c0 0.0 ρ0 α0 0.0
           # h0  c0 0.0 ρ0 α0 0.0]

    layers = [0.0 0.0 h0]

    T2 = promote_type(typeof(h0), typeof(cb), typeof(ρb), typeof(αb))
    sspHS = zeros(T2, 2, 6)
    sspHS[1, :] = [0.0 343.0 0.0 0.00121 0.0 0.0]
    sspHS[2, :] = [h0 cb 0.0 ρb αb 0.0]

    # env_dict = Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
    return ssp, layers, sspHS
end


### Standard 1-layer sediment model with constant sound speeds
function one_layer_test_dict_KRAKEN(
        ;c0 = 1500.0, c1 = 1550.0, cb = 1600.0, ρ0 = 1000.0, ρ1 = 1500.0, ρb = 2000.0, h0 = 100.0, h1 = 20.0)
    # Water column
    α0 = 0.0
    # sediment layer
    α1 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    freq = 100.0
    z0 = h0
    z1 = h0 + h1

    ssp = NamedArray(
        [0.0 c0 0.0 ρ0 α0 0.0
         z0 c0 0.0 ρ0 α0 0.0
         z0+eps(z0) c1 0.0 ρ1 α1 0.0
         z1 c1 0.0 ρ1 α1 0.0];
        dimnames = ("i", "Paramater")
    )
    setnames!(ssp, ["z", "cp", "cs", "ρ", "αp", "αs"], 2)
    ssp = [0.0 c0 0.0 ρ0 α0 0.0
           z0 c0 0.0 ρ0 α0 0.0
           z0+eps(z0) c1 0.0 ρ1 α1 0.0
           z1 c1 0.0 ρ1 α1 0.0]

    layers = [0.0 0.0 z0
              0.0 0.0 z1]

    sspHS = [0.0 343.0 0.0 0.00121 0.0 0.0
             z1 cb 0.0 ρb αb 0.0]

    env_dict = Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
    return env_dict
end

### Standard 1-layey model with slope in sound speed
function one_layer_slope_test_dict_KRAKEN(;c0 = 1500.0, c1_1 = 1550.0, c1_2 = 1580.0, cb = 1600.0,
        ρ0 = 1000.0, ρ1 = 1500.0, ρb = 2000.0, h0 = 100.0, h1 = 20.0)
    # Water column
    α0 = 0.0
    # sediment layer
    α1 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    freq = 100.0
    z0 = h0
    z1 = h0 + h1

    ssp = NamedArray(
        [0.0 c0 0.0 ρ0 α0 0.0
         z0 c0 0.0 ρ0 α0 0.0
         z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
         z1 c1_2 0.0 ρ1 α1 0.0];
        dimnames = ("i", "Paramater")
    )
    setnames!(ssp, ["z", "cp", "cs", "ρ", "αp", "αs"], 2)

    ssp = [0.0 c0 0.0 ρ0 α0 0.0
           z0 c0 0.0 ρ0 α0 0.0
           z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
           z1 c1_2 0.0 ρ1 α1 0.0]

    layers = [0.0 0.0 z0
              0.0 0.0 z1]

    sspHS = [0.0 343.0 0.0 0.00121 0.0 0.0
             z1 cb 0.0 ρb αb 0.0]

    return Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
end

### Standard 2-layer model with slope in sound speed
function two_layer_slope_test_dict_KRAKEN(;
        c0 = 1500.0, c1_1 = 1550.0, c1_2 = 1580.0, c2_1 = 1600.0, c2_2 = 1650.0, cb = 1800.0,
        ρ0 = 1000.0, ρ1 = 1500.0, ρ2 = 1600.0, ρb = 2000.0, h0 = 100.0, h1 = 20.0, h2 = 20.0)
    # Water column
    α0 = 0.0
    # sediment layer
    α1 = 0.0
    # sediment layer
    α2 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    freq = 100.0
    z0 = h0
    z1 = h0 + h1
    z2 = h0 + h1 + h2

    ssp = NamedArray(
        [0.0 c0 0.0 ρ0 α0 0.0
         z0 c0 0.0 ρ0 α0 0.0
         z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
         z1 c1_2 0.0 ρ1 α1 0.0
         z1+eps(z1) c2_1 0.0 ρ2 α2 0.0
         z2 c2_2 0.0 ρ2 α2 0.0];
        dimnames = ("i", "Paramater")
    )
    setnames!(ssp, ["z", "cp", "cs", "ρ", "αp", "αs"], 2)

    ssp = [0.0 c0 0.0 ρ0 α0 0.0
           z0 c0 0.0 ρ0 α0 0.0
           z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
           z1 c1_2 0.0 ρ1 α1 0.0
           z1+eps(z1) c2_1 0.0 ρ2 α2 0.0
           z2 c2_2 0.0 ρ2 α2 0.0]

    layers = [0.0 0.0 z0
              0.0 0.0 z1
              0.0 0.0 z2]

    sspHS = [0.0 343.0 0.0 0.00121 0.0 0.0
             z2 cb 0.0 ρb αb 0.0]

    return Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
end
