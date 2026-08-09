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
#
# A third group of rules is here for a different reason: `soundspeed` and `density` *are* traceable,
# but what `DataInterpolations` supplies for them is silently incomplete — its ChainRules rule holds
# an interpolant's knot vector fixed, so every derivative with respect to a layer thickness comes back
# missing a term. That is a wrong number rather than an error, so the rules below replace it outright.

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

### The sampled-profile interpolants -------------------------------------------------------------
#
# `soundspeed(ssp, z)` and `density(ρ, z)` are `ssp.f(z)`, a `DataInterpolations.LinearInterpolation`
# with `ExtrapolationType.Constant`. `DataInterpolations` ships a ChainRules rule for that call, but
# it differentiates only with respect to the *values* and the query point — an interpolant's knots
# are treated as constants. Here the knots are layer boundaries, and layer boundaries are parameters,
# so the missing term is not small: on `two_layer_slope_env` it made `∂kr₁/∂h0` 15% wrong and flipped
# the sign of `∂kr₁/∂h1` and `∂kr₁/∂h2`, with nothing raised. Unmodified `kraken.exe`, finite-
# differenced across perturbed `.env` files, was the arbiter; the table is in the plan.
#
# The three rules below close that. Linear interpolation is local — one query touches exactly two
# knots — so the whole adjoint is a scatter into two slots per query, O(N) for N queries, and the
# constructor rule reduces to storing its arguments. That last one is also where the milestone's
# performance went: tracing the two `DataInterpolations` constructors inside `UnderwaterEnv` was 93%
# of a reverse gradient's cost, and a rule that never enters them removes it.

"""
    linear_interp_partials(A, z)

The value of the linear interpolant `A` at `z`, together with the partial derivatives of that value
with respect to the two knots it lands between and their two values.

Returns `(val, idx, du1, du2, dt1, dt2, dz)`, where `idx`/`idx+1` index `A.t` and `A.u`, so a pullback
scatters `du1`/`dt1` into slot `idx` and `du2`/`dt2` into slot `idx + 1`.

# Matching the primal rather than re-deriving it

`val` is recomputed here rather than taken from `A(z)` so that the two can be compared: the rule is a
derivative of the primal only if it selected the same interval, and with duplicated knots — which
this package uses deliberately, placing `z0` and `z0 + eps(z0)` to represent a discontinuity — the
interval on either side of a boundary carries a completely different slope. `test/reverse_ad_tests.jl`
asserts `val == A(z)` bit for bit across every standard environment for exactly that reason.

Three cases, all inherited from `DataInterpolations` instead of chosen:

  * **Constant extrapolation** outside `[t₁, t_N]` — the value is the nearer endpoint's, so it carries
    the whole cotangent and the knots and the query point carry none. The comparisons are strict
    (`z < t₁`, `z > t_N`), so *at* the last knot the interior branch runs and `∂val/∂z` is the last
    interval's slope, not zero.
  * **Coincident knots** (`t_{i+1} == t_i`) — `linear_interpolation_parameters` forces the slope to
    zero rather than dividing, making `val = u_i` a constant. Every partial but `∂val/∂u_i` is zero.
  * **Interior** — with `w = (z - t_i)/Δt` and `slope = Δu/Δt`:
    `∂val/∂u_i = 1 - w`, `∂val/∂u_{i+1} = w`, `∂val/∂t_i = -slope(1 - w)`,
    `∂val/∂t_{i+1} = -slope·w`, `∂val/∂z = slope`.

# A query exactly on a knot takes the interval *ending* there

Which of the two adjoining intervals a knot-coincident query belongs to changes no value — both
branches evaluate to the shared knot's own `u`, bit for bit on every mesh point of every standard
environment — but it changes the derivative completely, and here it is not a detail: `get_z_vec`
ends each layer's mesh exactly on that layer's lower boundary, so *every* interface knot is queried
exactly. Taking the interval that *starts* at the knot means taking the `eps(z0)`-wide seam this
package uses to spell a discontinuity, whose slope is ~1e15 and whose two huge partials then have to
cancel through unrelated paths in reverse mode; on `two_layer_slope_env` that left `∂kr₁/∂h0` 8.6%
wrong. Taking the interval that ends there gives the water column's own slope, which is the answer
the physics wants and the one ForwardDiff already produced.

That agreement with ForwardDiff is not a coincidence but it *is* an accident: `searchsortedlast` puts
a tie in the right-hand interval, and ForwardDiff only lands in the left-hand one because its `isless`
on `Dual` breaks value ties on the partials, so a `Float64` query compares below an equal-valued
`Dual` knot. Relying on that would be relying on a comparison ForwardDiff does not document, so the
side is selected explicitly here — `side = :first, idx_shift = -1`, which is `searchsortedfirst - 1`
and identical to `searchsortedlast` at every query that is *not* exactly on a knot.
"""
function linear_interp_partials(A, z)
    t, u = A.t, A.u
    N = length(t)
    T = promote_type(eltype(t), eltype(u), typeof(z))
    if z < first(t)
        # Slots 1 and 2, with all the weight on the first — `u[1]` is the extrapolated value.
        return (first(u) + zero(T), 1, one(T), zero(T), zero(T), zero(T), zero(T))
    elseif z > last(t)
        return (last(u) + zero(T), N - 1, zero(T), one(T), zero(T), zero(T), zero(T))
    end
    idx = DataInterpolations.get_idx(A, z, A.iguesser; side=:first, idx_shift=-1)
    t1, t2 = t[idx], t[idx + 1]
    u1, u2 = u[idx], u[idx + 1]
    Δt = t2 - t1
    if iszero(Δt)
        return (u1 + zero(T), idx, one(T), zero(T), zero(T), zero(T), zero(T))
    end
    slope = (u2 - u1) / Δt
    w = (z - t1) / Δt
    return (u1 + slope * (z - t1), idx, 1 - w, w, -slope * (1 - w), -slope * w, slope)
end

"""
    profile_pullback(A, zs, parts, Δ)

Scatter a cotangent `Δ` on interpolated values back onto `(A.u, A.t, zs)`, given the per-query
partials `parts` that [`linear_interp_partials`](@ref) produced on the forward pass.

Returns `(Δu, Δt, Δz)`, all three the shape of what they are a cotangent for. `Δz` follows `zs`: a
scalar query gets a scalar back, a vector query a vector.
"""
function profile_pullback(A, zs, parts, Δ)
    N = length(A.t)
    T = promote_type(eltype(A.t), eltype(A.u), eltype(zs), eltype(Δ))
    Δu = zeros(T, N)
    Δt = zeros(T, N)
    Δz = zeros(T, length(parts))
    # A pullback's own arithmetic is never re-differentiated, so accumulating in place is free here —
    # the no-mutation discipline is about the *primal*, which is what the tape follows.
    for (k, p) in enumerate(parts)
        d = Δ isa Number ? Δ : Δ[k]
        _, idx, du1, du2, dt1, dt2, dz = p
        Δu[idx] += d * du1
        Δu[idx + 1] += d * du2
        Δt[idx] += d * dt1
        Δt[idx + 1] += d * dt2
        Δz[k] = d * dz
    end
    return Δu, Δt, zs isa Number ? Δz[1] : Δz
end

"""
    profile_partials(prof, z)

The forward half shared by the [`soundspeed`](@ref) and [`density`](@ref) rules: the interpolated
value at `z`, plus the per-query partials the pullback scatters. `z` may be a number or a vector; the
value returned follows it, exactly as the primal's does.
"""
function profile_partials(prof, z::Number)
    p = linear_interp_partials(prof.f, z)
    return first(p), (p,)
end

function profile_partials(prof, z::AbstractVector)
    parts = map(zk -> linear_interp_partials(prof.f, zk), z)
    return map(first, parts), parts
end

# The cotangent goes onto the *struct's* fields — `z` and `c`/`ρ` — and never onto `f`, the
# interpolant: `f` is fully determined by the other two, and routing through it is precisely what
# would hand the derivative back to `DataInterpolations` and lose the knot term again.
#
# Note the sign. The struct stores `z = -depth` while the interpolant it built is keyed on `depth`,
# so a cotangent on a knot enters the struct negated. `test/reverse_ad_tests.jl` pins this against
# ForwardDiff, because getting it backwards is invisible in any environment whose profile is constant
# within each layer — which is most of them, and is why the gap survived 4.2 and 4.3.

function ChainRulesCore.rrule(::typeof(soundspeed), ssp::SampledSSP1D, z)
    val, parts = profile_partials(ssp, z)
    function soundspeed_pullback(Δ)
        Δ = unthunk(Δ)
        # A structural zero short-circuits: `profile_pullback` would otherwise widen its accumulators
        # to `Any` trying to promote against it, and scatter zeros for no reason.
        Δ isa AbstractZero && return (NoTangent(), ZeroTangent(), ZeroTangent())
        Δu, Δt, Δz = profile_pullback(ssp.f, z, parts, Δ)
        return (NoTangent(), Tangent{typeof(ssp)}(; z=(-Δt), c=Δu, f=NoTangent()), Δz)
    end
    return val, soundspeed_pullback
end

function ChainRulesCore.rrule(::typeof(density), prof::SampledDensity1D, z)
    val, parts = profile_partials(prof, z)
    function density_pullback(Δ)
        Δ = unthunk(Δ)
        Δ isa AbstractZero && return (NoTangent(), ZeroTangent(), ZeroTangent())
        Δu, Δt, Δz = profile_pullback(prof.f, z, parts, Δ)
        return (NoTangent(), Tangent{typeof(prof)}(; z=(-Δt), ρ=Δu, f=NoTangent()), Δz)
    end
    return val, density_pullback
end

# Building a profile is one `LinearInterpolation` call plus a negation, and reverse mode has no
# business inside the first of those — the rules above already state every derivative the interpolant
# has. Without these two rules the answers are still right (the tangent on `f` is `NoTangent`, so the
# interpolant is never differentiated), but Zygote still *traces* the constructor to build a pullback
# it will not use, and that tracing was 93% of a reverse gradient's cost.
#
# `f` here is the interpolation constructor, not an interpolant — `SampledSSP(depth, c)` passes
# `DataInterpolations.LinearInterpolation` itself — so it is a `NoTangent` on both sides.

"""
    profile_ctor_pullback(Δ, values_field)

Shared pullback body for the two profile constructors: map a tangent on the built struct back onto
`(depth, values)`. `z = -depth` is the only field that is not stored verbatim.
"""
function profile_ctor_pullback(Δ, values_field::Symbol)
    Δ = unthunk(Δ)
    Δ isa AbstractZero && return (NoTangent(), ZeroTangent(), ZeroTangent(), NoTangent())
    Δz = getproperty(Δ, :z)
    Δdepth = Δz isa AbstractZero ? ZeroTangent() : -Δz
    return (NoTangent(), Δdepth, getproperty(Δ, values_field), NoTangent())
end

function ChainRulesCore.rrule(::Type{SampledSSP1D}, depth, c, f)
    SampledSSP1D_pullback(Δ) = profile_ctor_pullback(Δ, :c)
    return SampledSSP1D(depth, c, f), SampledSSP1D_pullback
end

function ChainRulesCore.rrule(::Type{SampledDensity1D}, depth, ρ, f)
    SampledDensity1D_pullback(Δ) = profile_ctor_pullback(Δ, :ρ)
    return SampledDensity1D(depth, ρ, f), SampledDensity1D_pullback
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
