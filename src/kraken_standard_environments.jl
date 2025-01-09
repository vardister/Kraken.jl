using DrWatson
using NamedArrays

### Standard Pekeris
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

  return Dict("c" => c, "cb" => cb, "ρ" => ρ, "ρb" => ρb, "depth" => depth, "freq" => freq)
end

function pekeris_test_dict_KRAKEN()
  d = pekeris_test_dict()
  # Water column
  c0 = d["c"]
  ρ0 = d["ρ"]
  h0 = d["depth"]
  α0 = 0.0
  # bottom half-space
  cb = d["cb"]
  ρb = d["ρb"]
  αb = 0.0
  # other
  freq = 100.0
  z0 = h0
  depth = h0

  ssp = NamedArray(
    [
      0.0 c0 0.0 ρ0 α0 0.0
      h0 c0 0.0 ρ0 α0 0.0
    ];
    dimnames=("i", "Paramater"),
  )
  setnames!(ssp, ["z", "cp", "cs", "ρ", "αp", "αs"], 2)
  ssp = [
    0.0 c0 0.0 ρ0 α0 0.0
    h0 c0 0.0 ρ0 α0 0.0
  ]

  layers = [0.0 0.0 h0]

  sspHS = [
    0.0 343.0 0.0 0.00121 0.0 0.0
    z0 cb 0.0 ρb αb 0.0
  ]

  env_dict = Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
  return env_dict
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
  depth = h0 + h1
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
    [
      0.0 c0 0.0 ρ0 α0 0.0
      z0 c0 0.0 ρ0 α0 0.0
      z0+eps(z0) c1 0.0 ρ1 α1 0.0
      z1 c1 0.0 ρ1 α1 0.0
    ];
    dimnames=("i", "Paramater"),
  )
  setnames!(ssp, ["z", "cp", "cs", "ρ", "αp", "αs"], 2)
  ssp = [
    0.0 c0 0.0 ρ0 α0 0.0
    z0 c0 0.0 ρ0 α0 0.0
    z0+eps(z0) c1 0.0 ρ1 α1 0.0
    z1 c1 0.0 ρ1 α1 0.0
  ]

  layers = [
    0.0 0.0 z0
    0.0 0.0 z1
  ]

  sspHS = [
    0.0 343.0 0.0 0.00121 0.0 0.0
    z1 cb 0.0 ρb αb 0.0
  ]

  env_dict = Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
  return env_dict
end

### Standard 1-layey model with slope in sound speed
function one_layer_slope_test_dict()
  # Water column
  c0 = 1500.0
  ρ0 = 1000.0
  h0 = 100.0
  # sediment layer
  c1_1 = 1550.0
  c1_2 = 1580.0
  ρ1 = 1500.0
  h1 = 20.0
  # bottom half-space
  cb = 1600.0
  ρb = 2000.0
  # other
  depth = h0 + h1
  freq = 100.0

  return @dict c0 c1 cb ρ0 ρ1 ρb h0 h1 depth freq
end

function one_layer_slope_test_dict_KRAKEN()
  # Water column
  c0 = 1500.0
  ρ0 = 1000.0
  h0 = 100.0
  α0 = 0.0
  # sediment layer
  c1_1 = 1550.0
  c1_2 = 1580.0
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
    [
      0.0 c0 0.0 ρ0 α0 0.0
      z0 c0 0.0 ρ0 α0 0.0
      z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
      z1 c1_2 0.0 ρ1 α1 0.0
    ];
    dimnames=("i", "Paramater"),
  )
  setnames!(ssp, ["z", "cp", "cs", "ρ", "αp", "αs"], 2)

  ssp = [
    0.0 c0 0.0 ρ0 α0 0.0
    z0 c0 0.0 ρ0 α0 0.0
    z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
    z1 c1_2 0.0 ρ1 α1 0.0
  ]

  layers = [
    0.0 0.0 z0
    0.0 0.0 z1
  ]

  sspHS = [
    0.0 343.0 0.0 0.00121 0.0 0.0
    z1 cb 0.0 ρb αb 0.0
  ]

  return Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
end

### Standard 2-layer model with slope in sound speed
function two_layer_slope_test_dict()
  # Water column
  c0 = 1500.0
  ρ0 = 1000.0
  h0 = 100.0
  # sediment layer
  c1_1 = 1550.0
  c1_2 = 1580.0
  ρ1 = 1500.0
  h1 = 20.0
  # sediment layer
  c2_1 = 1600.0
  c2_2 = 1650.0
  ρ2 = 1600.0
  h2 = 20.0
  # bottom half-space
  cb = 1800.0
  ρb = 2000.0
  # other
  depth = h0 + h1 + h2
  freq = 100.0

  return @dict c0 c1_1 c1_2 c2_1 c2_2 cb ρ0 ρ1 ρ2 ρb h0 h1 h2 depth freq
end

function two_layer_slope_test_dict_KRAKEN()
  # Water column
  c0 = 1500.0
  ρ0 = 1000.0
  h0 = 100.0
  α0 = 0.0
  # sediment layer
  c1_1 = 1550.0
  c1_2 = 1580.0
  ρ1 = 1500.0
  h1 = 20.0
  α1 = 0.0
  # sediment layer
  c2_1 = 1600.0
  c2_2 = 1650.0
  ρ2 = 1600.0
  h2 = 20.0
  α2 = 0.0
  # bottom half-space
  cb = 1800.0
  ρb = 2000.0
  αb = 0.0
  # other
  freq = 100.0
  z0 = h0
  z1 = h0 + h1
  z2 = h0 + h1 + h2

  ssp = NamedArray(
    [
      0.0 c0 0.0 ρ0 α0 0.0
      z0 c0 0.0 ρ0 α0 0.0
      z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
      z1 c1_2 0.0 ρ1 α1 0.0
      z1+eps(z1) c2_1 0.0 ρ2 α2 0.0
      z2 c2_2 0.0 ρ2 α2 0.0
    ];
    dimnames=("i", "Paramater"),
  )
  setnames!(ssp, ["z", "cp", "cs", "ρ", "αp", "αs"], 2)

  ssp = [
    0.0 c0 0.0 ρ0 α0 0.0
    z0 c0 0.0 ρ0 α0 0.0
    z0+eps(z0) c1_1 0.0 ρ1 α1 0.0
    z1 c1_2 0.0 ρ1 α1 0.0
    z1+eps(z1) c2_1 0.0 ρ2 α2 0.0
    z2 c2_2 0.0 ρ2 α2 0.0
  ]

  layers = [
    0.0 0.0 z0
    0.0 0.0 z1
    0.0 0.0 z2
  ]

  sspHS = [
    0.0 343.0 0.0 0.00121 0.0 0.0
    z2 cb 0.0 ρb αb 0.0
  ]

  return Dict(:ssp => ssp, :layers => layers, :sspHS => sspHS, :freq => freq)
end
