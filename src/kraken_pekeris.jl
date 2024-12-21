using NonlinearSolve
using UnPack


struct PekerisUnderwaterEnv{T1 <: Real}
  c1::T1
  cb::T1
  ρ1::T1
  ρb::T1
  depth::T1
end

using NaNMath; nm=NaNMath
function find_kr(env::PekerisUnderwaterEnv, freq; n_points=5_000, method = Ridder())
  @unpack c1, cb, ρ1, ρb, depth = env
  ω = 2π*freq

  kr_max = max(ω / c1, ω / cb)
  kr_min = min(ω / c1, ω / cb)

  func(kr) = @. nm.tan(nm.sqrt((ω / c1)^2 - kr^2) * depth) +
     (ρb / ρ1) * nm.sqrt((ω / c1)^2 - kr^2) / nm.sqrt(kr^2 - (ω / cb)^2)
  func(kr, p) = func(kr)


  kr_try = range(kr_min + eps(kr_min), kr_max - eps(kr_min); length = n_points)
  func_try = func.(kr_try)
  idxes = [i for i = 1:length(func_try)-1 if (func_try[i] > 0 && func_try[i+1] < 0)]
  krspans = [(kr_try[i], kr_try[i+1]) for i in idxes]
  solutions = Vector{eltype(kr_try)}(undef, length(krspans))
  for (ii, krspan) in enumerate(krspans)
    prob = IntervalNonlinearProblem(func, krspan)
    solutions[ii] = solve(prob, method).u
  end
  return reverse(solutions)
end