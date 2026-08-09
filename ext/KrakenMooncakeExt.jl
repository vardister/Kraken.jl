"""
Mooncake support for the reverse-mode rules in `src/kraken_ad.jl`.

Mooncake does not read `ChainRulesCore.rrule`s on its own — it derives its own rules from the IR
unless a signature is declared a primitive. Left to itself it walks straight into the solver's
`range`/`TwicePrecision` arithmetic and stops with a `bitcast` error, and even if it got past that it
would be taping the bracketing solver and the inverse iteration, which is exactly what this milestone
exists to avoid. So every signature that carries a rule in `kraken_ad.jl` is declared a Mooncake
primitive here, and `Mooncake.@from_rrule` reuses the ChainRules rule behind it. Nothing in this file
states a derivative; it is a translation layer, and `test/reverse_ad_tests.jl` checks Mooncake against
ForwardDiff on the same targets as Zygote.

# The tangent bridge

`@from_rrule` converts between Mooncake's tangent representation and ChainRules' on every call, and
the conversion it ships covers scalars, arrays and tuples but not `ChainRulesCore.Tangent`s over
composite types — which is what all three of this package's rules return, since their arguments are
`UnderwaterEnv`, `AcousticProblemProperties`, `AcousticProblemCache` and the two sampled profiles.
The `increment_and_get_rdata!` methods below add that case. They are a mechanical field-by-field
recursion, but they do reach into Mooncake's `FData`/`RData`/`MutableTangent` split, which is
internal — hence the `Mooncake = "0.5"` compat bound in `Project.toml`.

Mooncake splits a tangent in two: `fdata` is the part that is mutated in place during the reverse
pass (arrays), and `rdata` is the part that is returned and accumulated (scalars). A struct therefore
has *both*, keyed by field: `UnderwaterEnv`'s `h_vec` lives in fdata while its `cb` lives in rdata.
A ChainRules `Tangent` has neither split nor, necessarily, every field — `Tangent{E}(; cb=…, ρb=…)`
names only two, and `getproperty` returns `ZeroTangent()` for the rest. So the recursion walks the
*Mooncake* side's field list, pulls the matching ChainRules field, and lets the existing
`ZeroTangent`/`NoTangent` methods absorb everything the rule left out.
"""
module KrakenMooncakeExt

using Kraken
using Mooncake

using ChainRulesCore: ChainRulesCore, Tangent, ZeroTangent
using Mooncake: DefaultCtx, FData, MutableTangent, NoFData, NoRData, RData
import Mooncake: increment_and_get_rdata!

using Kraken:
    AcousticProblemCache,
    AcousticProblemProperties,
    SampledDensity1D,
    SampledSSP1D,
    UnderwaterEnv,
    bisection,
    density,
    layer_mesh,
    mode_eigenvector,
    soundspeed,
    solve_for_kr

### The tangent bridge ----------------------------------------------------------------------------

# A structural zero increments nothing and contributes nothing to the returned rdata. Mooncake ships
# this for `NoTangent` but not for `ZeroTangent`, and the rules here emit both: `NoTangent` for an
# argument that has no derivative at all (the `layer_mesh` point count), `ZeroTangent` for one that
# has a derivative which happens to vanish (the root bracket, which selects *which* root is found and
# not where it sits).
increment_and_get_rdata!(::Any, r, ::ZeroTangent) = r

# Both halves of Mooncake's split carry the struct's full field list, so either can supply the names.
_field_names(::FData{N}, ::RData) where {N} = fieldnames(N)
_field_names(::NoFData, ::RData{N}) where {N} = fieldnames(N)
_field_names(::FData{N}, ::NoRData) where {N} = fieldnames(N)

_fdata_field(f::FData, n::Symbol) = getfield(f.data, n)
_fdata_field(::NoFData, ::Symbol) = NoFData()
_rdata_field(r::RData, n::Symbol) = getfield(r.data, n)
_rdata_field(::NoRData, ::Symbol) = NoRData()

# Immutable struct with rdata fields — `UnderwaterEnv` (`cb`, `ρb`) and `AcousticProblemProperties`
# (`freq`). The rdata has to be rebuilt and returned; the fdata is incremented in place.
function increment_and_get_rdata!(f::Union{FData,NoFData}, r::RData, t::Tangent)
    names = _field_names(f, r)
    vals = map(n -> increment_and_get_rdata!(_fdata_field(f, n), _rdata_field(r, n), getproperty(t, n)), names)
    return RData(NamedTuple{names}(vals))
end

# Immutable struct whose every differentiable field is an array — the two sampled profiles. There is
# no rdata to rebuild, so the recursion is for its in-place effect only.
function increment_and_get_rdata!(f::FData, r::NoRData, t::Tangent)
    for n in _field_names(f, r)
        increment_and_get_rdata!(_fdata_field(f, n), NoRData(), getproperty(t, n))
    end
    return NoRData()
end

# Mutable struct — `AcousticProblemCache`. Its whole tangent is fdata, held as full tangents rather
# than an `FData` split, and every field the rules touch (`a_vec`, `e_vec`, `λ_scaling`) is a
# `Vector{Float64}`, whose tangent is a vector Mooncake increments in place. Iterating the *rule's*
# fields rather than the struct's keeps `A`, the assembled `Tridiagonal`, out of it: it is scratch
# space the rules never differentiate.
function increment_and_get_rdata!(f::MutableTangent, ::NoRData, t::Tangent)
    for n in propertynames(t)
        increment_and_get_rdata!(getfield(f.fields, n), NoRData(), getproperty(t, n))
    end
    return NoRData()
end

### The primitives --------------------------------------------------------------------------------
#
# One declaration per rule in `kraken_ad.jl`, in the same order. Signatures are concrete rather than
# `Any`-typed on purpose: Mooncake's own guidance is to claim a primitive only for the argument types
# the ChainRules rule is known to be correct for, and everything the solver reaches these with is
# `Float64`.

Mooncake.@from_rrule(DefaultCtx, Tuple{typeof(layer_mesh),Float64,Float64,Int})

Mooncake.@from_rrule(DefaultCtx, Tuple{typeof(soundspeed),SampledSSP1D,Float64})
Mooncake.@from_rrule(DefaultCtx, Tuple{typeof(soundspeed),SampledSSP1D,Vector{Float64}})
Mooncake.@from_rrule(DefaultCtx, Tuple{typeof(density),SampledDensity1D,Float64})
Mooncake.@from_rrule(DefaultCtx, Tuple{typeof(density),SampledDensity1D,Vector{Float64}})

Mooncake.@from_rrule(DefaultCtx, Tuple{Type{SampledSSP1D},Vector{Float64},Vector{Float64},Type})
Mooncake.@from_rrule(DefaultCtx, Tuple{Type{SampledDensity1D},Vector{Float64},Vector{Float64},Type})

# The counterpart of `kraken_ad.jl`'s `@non_differentiable bisection(...)`.
Mooncake.@zero_derivative(
    DefaultCtx, Tuple{typeof(bisection),UnderwaterEnv,AcousticProblemProperties,AcousticProblemCache}
)

# `true` for the kwargs: both of these are reached with and without them (`find_kr` forwards
# `method`, `inverse_iteration` forwards `reltol`), and `@from_rrule` emits the rule for each form.
Mooncake.@from_rrule(
    DefaultCtx,
    Tuple{typeof(solve_for_kr),Vector{Float64},UnderwaterEnv,AcousticProblemProperties,AcousticProblemCache},
    true
)

Mooncake.@from_rrule(
    DefaultCtx,
    Tuple{typeof(mode_eigenvector),Float64,UnderwaterEnv,AcousticProblemProperties,AcousticProblemCache},
    true
)

end
