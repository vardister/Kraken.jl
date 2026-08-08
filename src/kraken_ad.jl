### Reverse-mode automatic differentiation
#
# Forward mode (ForwardDiff) works straight through the solver and needs nothing from this file.
# Reverse mode does: the solver mutates its cache in place, brackets roots by integer mode counting,
# and refines them with an iterative solver — none of which a tape can or should follow.
#
# Instead of making the solver traceable, this file attaches `ChainRulesCore` rules at the two seams
# identified at the top of `kraken_core.jl`, which makes everything behind them opaque:
#
#   * `solve_for_kr`      — implicit function theorem on the Sturm determinant  (this file)
#   * `inverse_iteration` — eigenvector adjoint                                 (task 4.3)
#
# One set of rules serves Zygote, Mooncake, and Enzyme-via-ChainRules at once.

using ChainRulesCore

# `bisection` returns wavenumber brackets derived from integer mode counts. The counts are piecewise
# constant in the parameters, so the brackets are too — and the root that `solve_for_kr` converges to
# does not depend on where the bracket's endpoints sit anyway, only on which root they enclose.
@non_differentiable bisection(::Any, ::Any, ::Any)

# The mesh point count is an integer and piecewise constant in the parameters. It has to be declared:
# the `ceil` inside it is an LLVM intrinsic with no derivative, so reverse mode errors out on it
# rather than noticing that an `Int` cannot carry one. `Δz_vec = h ./ Nz_vec` around it *is*
# differentiable and stays traced — see the seam note in `kraken_core.jl`.
@non_differentiable n_mesh_points(::Any, ::Any, ::Any)

# `maxsoundspeed` exists to feed the `maxsoundspeed(env.c) < env.cb` assertion in `get_Nz_vec`, and
# nothing on the seam consumes its value. Left traced, its `maximum` pullback reaches the assertion's
# `Bool` and reverse mode fails trying to accumulate one into a thunk.
@non_differentiable maxsoundspeed(::Any)

"""
    rrule(::typeof(layer_mesh), z_top, z_bot, Nz)

Reverse-mode rule for one layer's mesh.

[`layer_mesh`](@ref) is `collect(range(z_top, z_bot, Nz))`, and `range` is the one thing on the seam
that cannot simply be rewritten: its `StepRangeLen` carries `Base.TwicePrecision` fields that
reverse mode has no adjoint for, while any hand-rolled replacement shifts interior points by a few
ulp — enough, measured, to break the cross-validation tolerances the solver is held to. So the
primal keeps `range` and the derivative is stated here, where it is exact anyway: the mesh is a
linear interpolation between its endpoints, `z_j = (1 - t_j) z_top + t_j z_bot` with
`t_j = (j - 1) / (Nz - 1)`.
"""
function ChainRulesCore.rrule(::typeof(layer_mesh), z_top, z_bot, Nz)
    t = range(0.0, 1.0, Nz)
    function layer_mesh_pullback(Δz)
        Δ = unthunk(Δz)
        return (NoTangent(), sum((1 .- t) .* Δ), sum(t .* Δ), NoTangent())
    end
    return layer_mesh(z_top, z_bot, Nz), layer_mesh_pullback
end

"""
    sturm_sensitivities(kr, env, props, cache)

Partial derivatives of the (rescaled) Sturm determinant `D` with respect to everything it depends on,
computed at `kr` in one O(N) forward sweep plus one O(N) reverse sweep.

Returns a `NamedTuple`:

| Field | Meaning |
|---|---|
| `D` | the determinant itself — identical to `first(det_sturm(kr, env, props, cache))` |
| `dkr` | `∂D/∂kr`, including the path through [`get_g`](@ref) |
| `da`, `de`, `dλ` | `∂D/∂cache.a_vec`, `.e_vec`, `.λ_scaling` — one entry per mesh point |
| `dcb`, `dρb`, `dfreq` | `∂D/∂` the half-space parameters, via `g` only |

# Why this is written out rather than obtained by nesting AD

The rule below needs `∂D/∂θ` for *every* entry of the three cache vectors, because that is how the
environment's influence reaches the determinant. Forward-mode over `det_sturm` would cost one sweep
per entry — O(N²) for an O(N) problem — which would defeat the point of reverse mode entirely. The
Sturm recurrence is a product of 2×2 transfer matrices, so its adjoint is a backward sweep of the
same length as the forward one, and that is what this does.

# Why the rescaling is harmless

`det_sturm` multiplies the sequence by `scale_const`'s factor `s` whenever it would overflow or
underflow. `s` is piecewise constant in both `kr` and the parameters, so it contributes no derivative
and cancels exactly in the implicit-function ratio: `-(∂(S·D)/∂θ)/(∂(S·D)/∂kr) = -D_θ/D_kr`. The
forward sweep here records each `s` and the reverse sweep treats it as the constant it is — the
resulting `D` is the *rescaled* determinant, matching `det_sturm` bit for bit, and every derivative
carries the same `S` factor so the ratio is unaffected.
"""
function sturm_sensitivities(kr, env, props, cache)
    a = cache.a_vec
    e = cache.e_vec
    L = cache.λ_scaling
    N = sum(props.Nz_vec)
    T = promote_type(typeof(kr), eltype(a), eltype(e), eltype(L))

    q = 2pi * props.freq / env.cb        # horizontal wavenumber at the bottom cutoff
    β = sqrt(kr^2 - q^2)                 # vertical decay rate in the half-space
    g = β / env.ρb                       # == get_g(kr, env, props)

    # --- Forward sweep -------------------------------------------------------------------------
    # State at the start of step k is `v_k = (p0, p1)`; the step is `v_{k+1} = s_k * M_k * v_k` with
    # `M_k = [0 1; -e_k^2  (λ_k - α_k)]`. Record `v_k` and `s_k` for the reverse sweep.
    P0 = Vector{T}(undef, N)
    P1 = Vector{T}(undef, N)
    S = Vector{T}(undef, N)
    p0 = zero(T)
    p1 = one(T)
    D = zero(T)
    for k in 1:N
        P0[k] = p0
        P1[k] = p1
        # The last mesh point carries the bottom half-space boundary condition, exactly as in
        # `det_sturm`: `α = 0.5 * a[N] - g` instead of `α = a[k]`.
        α = k == N ? 0.5 * a[k] - g : a[k]
        p2 = (kr^2 * L[k] - α) * p1 - e[k]^2 * p0
        s = scale_const(p1, p2)
        S[k] = s
        p1 *= s
        p2 *= s
        p0 = p1
        p1 = p2
        D = p2
    end

    # --- Reverse sweep -------------------------------------------------------------------------
    # `(w0, w1) = ∂D/∂v_k`. `D` is the second component of `v_{N+1}`, so the sweep starts at (0, 1).
    da = zeros(T, N)
    de = zeros(T, N)
    dλ = zeros(T, N)
    dg = zero(T)
    dkr = zero(T)
    w0 = zero(T)
    w1 = one(T)
    for k in N:-1:1
        b = w1 * S[k]  # multiplier on everything inside the `p2` expression at step k
        α = k == N ? 0.5 * a[k] - g : a[k]
        dα = k == N ? T(0.5) : one(T)  # ∂α/∂a[k]

        da[k] = -b * dα * P1[k]
        k == N && (dg += b * P1[k])   # α = 0.5a - g, so ∂/∂g flips sign relative to ∂/∂a
        dλ[k] = b * kr^2 * P1[k]
        de[k] = -2 * b * e[k] * P0[k]
        dkr += 2 * b * kr * L[k] * P1[k]

        w0, w1 = -b * e[k]^2, w0 * S[k] + b * (kr^2 * L[k] - α)
    end

    # --- The `g` path --------------------------------------------------------------------------
    # g = sqrt(kr^2 - q^2) / ρb with q = 2π·freq/cb. Written out rather than differentiated because
    # it is three scalar expressions; `test/reverse_ad_tests.jl` checks them against ForwardDiff.
    dkr += dg * kr / (env.ρb * β)
    dcb = dg * q^2 / (env.cb * env.ρb * β)
    dρb = -dg * g / env.ρb
    dfreq = -dg * q^2 / (props.freq * env.ρb * β)

    return (D=D, dkr=dkr, da=da, de=de, dλ=dλ, dcb=dcb, dρb=dρb, dfreq=dfreq)
end

"""
    rrule(::typeof(solve_for_kr), span, env, props, cache; method, kwargs...)

Reverse-mode rule for a single converged wavenumber, by the implicit function theorem.

The root `kr*` satisfies `D(kr*, θ) = 0`, where `D` is the Sturm determinant and `θ` is everything
the discretized problem depends on. Differentiating that identity gives

    ∂kr*/∂θ = -D_θ / D_kr

with both partials evaluated at the converged root — so the rule never looks inside the bracketing
solver, and its cost is one extra pass over the mesh rather than a tape of the iteration.

`θ` reaches `D` by two routes and the rule emits a cotangent for each:

  * through `cache.a_vec`, `.e_vec`, `.λ_scaling`, which is where the sound-speed profile, the
    density profile, and the mesh spacing enter. The AD backend carries these on through
    [`finite_difference_coefficients`](@ref) — that function is on the differentiable seam precisely
    so this hand-off works.
  * through `env.cb`, `env.ρb`, and `props.freq`, which enter directly via [`get_g`](@ref)'s
    half-space term and never pass through the cache.

`span` gets a `ZeroTangent`: it selects *which* root is found, not where that root sits.
"""
function ChainRulesCore.rrule(::typeof(solve_for_kr), span, env, props, cache; method=ITP(), kwargs...)
    kr = solve_for_kr(span, env, props, cache; method=method, kwargs...)
    s = sturm_sensitivities(kr, env, props, cache)

    function solve_for_kr_pullback(Δkr)
        # ∂kr/∂θ = -D_θ / D_kr, scaled by the incoming cotangent.
        f = -unthunk(Δkr) / s.dkr
        env_tangent = Tangent{typeof(env)}(; cb=f * s.dcb, ρb=f * s.dρb)
        props_tangent = Tangent{typeof(props)}(; freq=f * s.dfreq)
        cache_tangent = Tangent{typeof(cache)}(; a_vec=f .* s.da, e_vec=f .* s.de, λ_scaling=f .* s.dλ)
        return (NoTangent(), ZeroTangent(), env_tangent, props_tangent, cache_tangent)
    end

    return kr, solve_for_kr_pullback
end
