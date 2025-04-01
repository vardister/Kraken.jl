using NamedArrays

export pekeris_env
export one_layer_env
export one_layer_slope_test_dict_KRAKEN
export two_layer_slope_test_dict_KRAKEN
export munk_test_dict_KRAKEN

### Standard Pekeris
function pekeris_env(;
    c0::Real=1500.0, cb::Real=1600.0, ρ0::Real=1000.0, ρb::Real=1500.0, depth::Real=100.0
)

    # Water column
    α0 = 0.0
    # bottom half-space
    αb = 0.0

    T = promote_type(typeof(c0), typeof(cb), typeof(ρ0), typeof(ρb), typeof(depth))
    ssp = zeros(T, 2, 6)
    ssp[1, :] = [0.0 c0 0.0 ρ0 α0 0.0]
    ssp[2, :] = [depth c0 0.0 ρ0 α0 0.0]

    # ssp = [0.0 c0 0.0 ρ0 α0 0.0
    # depth  c0 0.0 ρ0 α0 0.0]

    layers = [0.0 0.0 depth]

    T2 = promote_type(typeof(depth), typeof(cb), typeof(ρb), typeof(αb))
    sspHS = zeros(T2, 2, 6)
    sspHS[1, :] = [0.0 343.0 0.0 0.00121 0.0 0.0]
    sspHS[2, :] = [depth cb 0.0 ρb αb 0.0]

    # env_dict = Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
    return ssp, layers, sspHS
end

### Standard 1-layer sediment model with constant sound speeds
function one_layer_env(;
    c0=1500.0, c1=1550.0, cb=1600.0, ρ0=1000.0, ρ1=1500.0, ρb=2000.0, h0=100.0, h1=20.0, g=0.0
)
    # Water column
    α0 = 0.0
    # sediment layer
    α1 = 0.0
    # bottom half-space
    αb = 0.0
    # other
    z0 = h0
    z1 = h0 + h1

    
    ssp = [
        0.0 c0 0.0 ρ0 α0 0.0
        z0 c0 0.0 ρ0 α0 0.0
        z0+eps(z0) c1 0.0 ρ1 α1 0.0
        z1 c1*g 0.0 ρ1 α1 0.0
    ]

    layers = [
        0.0 0.0 z0
        0.0 0.0 z1
    ]

    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        z1 cb 0.0 ρb αb 0.0
    ]

    return ssp, layers, sspHS
end



function munk_test_dict_KRAKEN()
    function c(z)
        ϵ = 0.00737
        zhat = 2 * (z - 1300.0) / 1300
        return 1500.0 * (1.0 + ϵ * (zhat - 1.0 + exp(-zhat)))
    end

    # Water column
    ρ0 = 1000.0
    α0 = 0.0
    # bottom half-space
    αb = 0.0
    ρb = 1500.0
    cb = 1600.0

    # other
    zvec = 0:100:5000
    cvec = c.(zvec)

    ssp = transpose(hcat([[x[1], x[2], 0.0, ρ0, α0, 0.0] for x in zip(zvec, cvec)]...))
    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        zvec[end] cb 0.0 ρb αb 0.0
    ]
    layers = [0.0 0.0 zvec[end]]

    return ssp, layers, sspHS
end
