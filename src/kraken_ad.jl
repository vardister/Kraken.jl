### Reverse-mode automatic differentiation
#
# Forward mode (ForwardDiff) works straight through the solver and needs nothing from this file.
# Reverse mode does: the solver mutates its cache in place, brackets roots by integer mode counting,
# and refines them with an iterative solver — none of which a tape can or should follow.
#
# Instead of making the solver traceable, this file attaches `ChainRulesCore` rules at the two seams
# identified at the top of `kraken_core.jl`, which makes everything behind them opaque:
#
#   * `solve_for_kr`      — implicit function theorem on the Sturm determinant
#   * `mode_eigenvector`  — eigenvector adjoint via a bordered tridiagonal solve
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
    h = half_space_sensitivities(dg, kr, env, props)

    return (D=D, dkr=dkr + h.dkr, da=da, de=de, dλ=dλ, dcb=h.dcb, dρb=h.dρb, dfreq=h.dfreq)
end

"""
    half_space_sensitivities(dg, kr, env, props)

Chain a cotangent `dg` on the bottom half-space term [`get_g`](@ref) back onto the four parameters it
is built from, returning a `NamedTuple` `(dkr, dcb, dρb, dfreq)`.

`g = √(kr² - q²) / ρb` with `q = 2π·freq/cb` is three scalar expressions, so both rules in this file
state its derivatives rather than differentiating it — `test/reverse_ad_tests.jl` checks them against
ForwardDiff.
"""
function half_space_sensitivities(dg, kr, env, props)
    q = 2pi * props.freq / env.cb
    β = sqrt(kr^2 - q^2)
    g = β / env.ρb
    return (
        dkr=dg * kr / (env.ρb * β),
        dcb=dg * q^2 / (env.cb * env.ρb * β),
        dρb=(-dg * g / env.ρb),
        dfreq=(-dg * q^2 / (props.freq * env.ρb * β)),
    )
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

"""
    shifted_matrix(kr, env, props, cache)

The finite-difference matrix at wavenumber `kr` — the same `M(kr, θ)` that
[`create_finite_diff_matrix!`](@ref) writes into `cache.A`, assembled here into fresh storage.

The rule below needs this matrix *after* the primal has run, and the primal restores the cache to its
unshifted state on the way out. Rebuilding costs one O(N) allocation and keeps the pullback from
mutating state it shares with the rest of the solve — a pullback runs long after the forward pass,
with no guarantee about what has touched the cache in between.
"""
function shifted_matrix(kr, env, props, cache)
    a, e, λ = cache.a_vec, cache.e_vec, cache.λ_scaling
    N = length(a)
    g = get_g(kr, env, props)
    # The last row carries the bottom half-space boundary condition, exactly as in the primal.
    d = [k == N ? 0.5 * a[N] - kr^2 * λ[N] - g : a[k] - kr^2 * λ[k] for k in 1:N]
    return Tridiagonal(e[2:N], d, e[2:N])
end

"""
    rrule(::typeof(mode_eigenvector), kr, env, props, cache; reltol)

Reverse-mode rule for a mode shape, by eigenvector perturbation theory.

The primal returns the unit-2-norm eigenvector `v` of the finite-difference matrix `M(kr, θ)`
belonging to its smallest-magnitude eigenvalue `μ` — which is zero to machine precision when `kr` is
a converged wavenumber, but is treated as an ordinary simple eigenvalue here, because `kr` and `θ`
are independent arguments and nothing constrains `M` to stay singular when they are perturbed
separately. Differentiating `M v = μ v` together with `vᵀv = 1` gives

    dv = -M⁺ (I - v vᵀ) dM v

with `M⁺` the pseudo-inverse on the complement of `v`, so the pullback of a cotangent `Δv` is a
single solve against that same complement,

    ⟨Δv, dv⟩ = -⟨y, dM v⟩,    y = M⁺ (I - v vᵀ) Δv,

i.e. a cotangent `ΔM = -y vᵀ`, read off on `M`'s tridiagonal sparsity pattern and chained back onto
`cache.a_vec`, `.e_vec`, `.λ_scaling`, `kr`, and — through [`get_g`](@ref) — `env.cb`, `env.ρb` and
`props.freq`. Cost is two tridiagonal solves, O(N), against the iteration's own dozens.

# The bordered solve, and why it is spelled as a projection

`y` is defined by the bordered system

    [ M   v ] [ y ]   [ Δv ]
    [ vᵀ  0 ] [ γ ] = [ 0  ]

whose matrix is nonsingular precisely because `v` spans `M`'s null space. Solving it as written would
give up the tridiagonal structure, so it is solved instead as *project, solve, project*: remove the
`v` component from the right-hand side, solve the (numerically nonsingular, since `μ ≈ 1e-13‖M‖`,
but very ill-conditioned in that one direction) tridiagonal system, and remove the `v` component the
solve's roundoff amplified back into the answer. One step of iterative refinement — again inside the
complement — cleans up what the amplification costs in the remaining directions.

# What the normalization is not

`v` here is normalized in the 2-norm only. The physical energy normalization lives in
[`normalize_mode`](@ref), which is on the differentiable seam and traced by the backend; it depends on
the density profile and the mesh coordinates, and pulling it inside this rule would mean carrying an
adjoint for the density interpolant that the backend already has. The sign fix (`v[1] ≥ 0`) is
piecewise constant in the parameters and contributes nothing — the relations above hold for `-v` just
as they do for `v`.
"""
function ChainRulesCore.rrule(::typeof(mode_eigenvector), kr, env, props, cache; reltol=0.01)
    kr_new, v = mode_eigenvector(kr, env, props, cache; reltol=reltol)

    function mode_eigenvector_pullback(Δ)
        # The primal returns a 2-tuple; either component can come back as a structural zero.
        Δkr_new = unthunk(Δ[1])
        Δv = unthunk(Δ[2])
        # `kr_new` is the iteration's own refinement of `kr` — it differs from it by the shift and by
        # the residual eigenvalue estimate, both at the level of roundoff, so it carries `kr`'s
        # derivative unchanged.
        dkr = Δkr_new isa AbstractZero ? zero(kr) : Δkr_new

        if Δv isa AbstractZero
            return (NoTangent(), dkr, ZeroTangent(), ZeroTangent(), ZeroTangent())
        end

        N = length(v)
        M = shifted_matrix(kr, env, props, cache)
        F = lu(M; check=false)
        P(u) = u .- dot(v, u) .* v          # projection onto the complement of `v`
        y = P(F \ P(Δv))
        y = P(y .+ (F \ P(Δv .- M * y)))    # one step of iterative refinement, in the complement

        # ΔM = -y vᵀ, restricted to the tridiagonal pattern. `e_vec[k+1]` sits at both (k, k+1) and
        # (k+1, k), so it collects both; `e_vec[1]` appears in the matrix nowhere.
        Δd = -y .* v
        Δa = [k == N ? 0.5 * Δd[N] : Δd[k] for k in 1:N]
        Δλ = (-kr^2) .* Δd
        Δe = vcat(zero(eltype(Δd)), -(y[2:N] .* v[1:(N - 1)] .+ y[1:(N - 1)] .* v[2:N]))

        # `d[N] = 0.5 a[N] - kr² λ[N] - g`, so `g` picks up the diagonal's cotangent with a sign flip.
        h = half_space_sensitivities(-Δd[N], kr, env, props)
        dkr += h.dkr - 2 * kr * dot(Δd, cache.λ_scaling)

        env_tangent = Tangent{typeof(env)}(; cb=h.dcb, ρb=h.dρb)
        props_tangent = Tangent{typeof(props)}(; freq=h.dfreq)
        cache_tangent = Tangent{typeof(cache)}(; a_vec=Δa, e_vec=Δe, λ_scaling=Δλ)
        return (NoTangent(), dkr, env_tangent, props_tangent, cache_tangent)
    end

    return (kr_new, v), mode_eigenvector_pullback
end
