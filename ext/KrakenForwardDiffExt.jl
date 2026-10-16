"""
ForwardDiff support for the wavenumber root: the implicit-function theorem in forward mode.

Without this extension ForwardDiff goes straight through `solve_for_kr`, and the `Dual` it returns
carries whatever partials the bracketing solver's last step left there. Most of ITP's steps are
bisections and projections, which copy an endpoint's partials. The correction therefore comes from a
final secant across a bracket a few ulps wide, where the determinant is ~1e-14 and its slope
estimate is rounding noise. Plan task 4.2 measured this at 1.3% on `∂kr₁/∂cb` and accepted it,
with finite differences as the arbiter. Task 6.11 found the same error in group speeds, which are
what the docs tell people to use ForwardDiff for: up to 1.6% on a shelf, at scattered frequencies.

# The rule

The root satisfies `D(kr, θ) = D₀`. Here `D` is the Sturm determinant and `D₀` is its value at the
converged root with every `Dual` stripped, which is zero up to the root tolerance. One Newton step in
`Dual` arithmetic,

    kr ← kr - (D(kr, θ) - D₀) / ∂D/∂kr(kr, θ),

leaves the value alone (the innermost correction is exactly `0 / ∂D/∂kr`) and sets the first-order
partials to `-D_θ / D_kr`, whatever they were before. That is the formula the reverse rule in
`src/kraken_ad.jl` uses, evaluated at the same point, so the two modes agree by construction. Each
further step doubles the order that is right, so a `Dual` nested `n` deep (a Hessian is 2) gets `n`
steps, which is more than enough.

`D` and `∂D/∂kr` come from `Kraken.sturm_and_slope`, one O(N) pass with no allocation.
"""
module KrakenForwardDiffExt

using Kraken
using ForwardDiff: Dual, value

nesting(::Type{<:Dual{T,V}}) where {T,V} = 1 + nesting(V)
nesting(::Type) = 0

innermost(x::Dual) = innermost(value(x))
innermost(x) = x

Kraken.primal_value(x::Dual) = innermost(x)

function Kraken.implicit_root(::Type{D}, kr, env, props, cache) where {D<:Dual}
    k = convert(D, kr)
    D0 = nothing
    for _ in 1:nesting(D)
        det, slope = Kraken.sturm_and_slope(k, env, props, cache)
        D0 === nothing && (D0 = innermost(det))
        k = k - (det - D0) / slope
    end
    return k
end

end
