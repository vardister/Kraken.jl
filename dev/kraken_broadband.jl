using Kraken
using DataInterpolations
using NamedArrays
import NaNMath as nm
using UnPack

function pressure_field(mod_sol, r, zs, zr; t0_offset = 0.1, max_modes = 41)
  freq = mod_sol.props.freq
  max_modes = min(size(mod_sol.modes, 2), max_modes)
  krs = mod_sol.kr[1:max_modes]
  modes = mod_sol.modes[:, 1:max_modes]
  zn = prepend!(vcat(mod_sol.props.zn_vec...), 0.0)
  modes_interpolated = [QuadraticInterpolation(mode, zn) for mode in eachcol(modes)]
  ϕ_zr = [interp(zr) for interp in modes_interpolated]
  ϕ_zs = [interp(zs) for interp in modes_interpolated]

  c_max = maximum(mod_sol.env.c.c)
  t0 = r / c_max - t0_offset  # align window correctly in time
  ρw = mod_sol.env.ρ.ρ[1]
  Q = 1im * exp(-1im * π / 4) / (ρw * sqrt(8π * r))
  pf = @. Q * ϕ_zs * ϕ_zr * exp(-im * krs * r) / sqrt(krs)
  pf = sum(pf)
  pf_shifted = pf * exp(2im * π * freq * t0)
  return pf_shifted
end

function pressure_field(mod_sol, r_vec::Vector, zs, zr_vec::Vector; t0_offset = 0.1,
    min_modes = 1, max_modes = typemax(Int), TL = false)
  freq = mod_sol.props.freq
  max_modes = min(size(mod_sol.modes, 2), max_modes)
  krs = mod_sol.kr[min_modes:max_modes]
  modes = mod_sol.modes[:, min_modes:max_modes]
  zn = vcat(mod_sol.props.zn_vec...)
  modes_interpolated = [QuadraticInterpolation(mode, zn) for mode in eachcol(modes)]

  if !TL
    T = promote_type(eltype(modes), ComplexF64)
  elseif TL
    T = promote_type(eltype(modes), Float64)
  end

  pf_field = NamedArray(
    zeros(T, length(zr_vec), length(r_vec)); names = (zr_vec, r_vec ./ 1000), dimnames = (
      "Depth(m)", "Range(km)")
  )
  c_max = maximum(mod_sol.env.c.c)
  ρw = mod_sol.env.ρ.ρ[1]
  ϕ_zs = [interp(zs) for interp in modes_interpolated]

  for ii in eachindex(zr_vec)
    # @infiltrate
    ϕ_zr = [interp(zr_vec[ii]) for interp in modes_interpolated]
    for jj in eachindex(r_vec)
      r = r_vec[jj]
      t0 = r / c_max - t0_offset  # align window correctly in time
      if !TL
        Q = 1im * exp(-1im * π / 4) / (ρw * sqrt(8π * r))
        pf = @. Q * ϕ_zs * ϕ_zr * exp(-im * krs * r) / sqrt(krs)
        pf = sum(pf)
        pf_field[ii, jj] = pf * exp(2im * π * freq * t0)
      elseif TL
        Q = 1 / (ρw * sqrt(r / 2pi))
        pf = @. Q * ϕ_zs * ϕ_zr * exp(-im * krs * r) / sqrt(krs)
        pf = sum(pf)
        pf_field[ii, jj] = -20 * log10.(abs.(pf * exp(2im * π * freq * t0)))
      end
    end
  end
  return pf_field
end

"""
    pressure_field_pekeris(env::PekerisUnderwaterEnv, freq, r, zs, zr; t0=0.1, n_points=20_000, max_modes=Inf)

Calculate the pressure at the `freq` frequency for a Pekeris waveguide defined by `env`.
The pressure is calculated at a range `r` for a source depth `zs` and a receiver depth `zr`.

# Optional Arguments
- `t0::Float64=0.1`: Time offset for the pressure field.
- `n_points::Int=2_000`: Number of discreet wavenumber points to find root intervals for the transcdental equation.
"""
function pressure_field_pekeris(
    env::PekerisUnderwaterEnv, freq, r, zs, zr; t0 = 0.1, n_points = 2_000, max_modes = Inf)
  if freq == 0
    return 0.0 + 0.0im
  end
  ω = 2π * freq
  p = [env.c1, env.cb, env.ρ1, env.ρb, env.depth]
  kr = find_kr(env, freq, p; n_points = n_points)
  if length(kr) == 0
    return 0.0 + 0.0im
  end
  @unpack c1, cb, ρ1, ρb, depth = env
  kzw = sqrt.((ω / c1)^2 .- kr .^ 2) # vertical wavenumber in the water column
  kzb = sqrt.(kr .^ 2 .- (ω / cb)^2) # vertical wavenumber in the bottom half-space
  # The amplitude `amplitude` of the modes
  amplitude = nm.sqrt.(
    (4 .* kzb .* kzw .* ρ1 .* ρb) ./
    (2 .* kzw .* ρ1 .* sin.(depth .* kzw) .^ 2 .-
     kzb .* ρb .* (-2 .* depth .* kzw .+ sin.(2 .* depth .* kzw)))
  )
  # The mode function
  function Ψ(z, mode)
    return amplitude[mode] * sin(kzw[mode] * z) * (z .<= depth) +
           amplitude[mode] * sin(kzw[mode] * depth) * exp(-kzb[mode] * (z - depth)) *
           (z .> depth)
  end

  t_offset = r / min(c1, cb) - t0  # align window correctly in time
  n_modes = Int(min(max_modes, length(kr)))
  pf = -1im * exp(-1im * π / 4) / (ρ1 * sqrt(8π * r)) *
       sum([Ψ(zs, mode) * Ψ(zr, mode) * exp(-1im * kr[mode] * r) / sqrt(kr[mode])
            for mode in range(1, n_modes)])
  return pf * exp(2im * π * freq * t_offset)
end
