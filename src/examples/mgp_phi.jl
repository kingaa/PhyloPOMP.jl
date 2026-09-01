# mgp_phi.jl
# =============================================================================
# Generic KLI compatibility-ratio math: production slots, saturation
# enumeration, and the binomial ratio phi_u.
#
#   This is the mathematical core flagged by M00/M01 as missing entirely:
#   `enumerate_saturations` (KLI Sec.3.4.2) and `phi_u` (KLI Eq.9, the
#   production-slot binomial ratio) currently exist only as eight
#   independent hand-derived, hand-transcribed closed forms (seir_naive.jl,
#   mers_naive.jl, and their _soft/_guided/_hard siblings) plus one
#   hand-worked LaTeX derivation (mers_filter_suite.tex). This file
#   implements the two pieces generically, as ordinary callable/testable
#   functions, from the Population IR (`Event.r`) alone.
#
#   Explicitly OUT of scope for this milestone (M02):
#     - deriving Q_u from genealogy/coloring semantics (M03) -- Q_u is an
#       external argument to `kli_binomial_ratio`, default 1.
#     - the coloring operator (chop/swap/fork) itself -- that consumes a
#       saturation, it is not part of computing one.
#     - reducing/marginalizing Phi_u = sum_m phi_u over the event-indicator
#       m (M04) -- this file stops at phi_u for a single given saturation.
#
# Primary source: King, Lin & Ionides, "Exact phylodynamic likelihood via
#   structured Markov genealogy processes" (StructuredMGPs.pdf), Eq.9 /
#   Sec.3.4.2, as excerpted/cross-checked against this repo's own worked
#   derivation `mers_filter_suite.tex` Sec."Construction of Compatibility
#   Terms" (lines 434-541) and its Complete Reference Table (lines 541-638).
#   No page/equation from the PDF itself is reproduced here -- only the
#   already-vetted restatement given in the M02 task spec.
# =============================================================================

export production_slots, enumerate_saturations, kli_binomial_ratio

"""
    production_slots(event::Event) -> Vector{Int}

Return `event`'s production vector `r_u` (KLI Sec.3.3.1's per-mark
constant), one entry per lineage-carrying deme, in `MGPModel.demes` order.

This is a thin, documented accessor -- `Event.r` was already populated by
the DSL (`mgp_macro.jl`'s `_production`) but, per M00's reconnaissance, was
never read by any downstream code. This is the first consumer.
"""
production_slots(event::Event) = event.r

"""
    enumerate_saturations(r::AbstractVector{<:Integer}, ℓ::AbstractVector{<:Integer})
        -> Vector{Vector{Int}}
    enumerate_saturations(event::Event, ℓ::AbstractVector{<:Integer})
        -> Vector{Vector{Int}}

Algorithmically enumerate the saturation space

    S_u(ℓ) = ∏_{d=1}^{D} {0, 1, ..., min(r_d, ℓ_d)}

for a production vector `r_u` (or an `Event`, from which `r_u =
production_slots(event)`) and a current per-deme lineage count `ℓ`. Each
element of the returned vector is one feasible saturation `s`, itself a
`Vector{Int}` of length `D = length(r)`, with `0 <= s_d <= min(r_d, ℓ_d)`
for every deme `d`.

This is a genuinely generic enumeration -- a nested/product iteration over
the per-deme feasible ranges -- not a per-model/per-event literal list.
`length(r) == length(ℓ)` is required (one entry per lineage-carrying deme);
mismatched lengths throw `ArgumentError`.

Boundary cases (see `handoffs/M02_kli_phi.md` for hand-verified examples):
- `ℓ_d == 0` for some deme `d` collapses that deme's range to `{0}` (a
  deme with no tracked lineages can never have a tracked-outcome slot).
- `r_d == 0` for some deme `d` likewise collapses that deme's range to
  `{0}` (nothing is produced there, so nothing can be saturated there).
"""
function enumerate_saturations(r::AbstractVector{<:Integer}, ℓ::AbstractVector{<:Integer})
    length(r) == length(ℓ) ||
        throw(ArgumentError("r and ℓ must have the same length (one entry " *
                             "per lineage-carrying deme); got length(r)=$(length(r)), " *
                             "length(ℓ)=$(length(ℓ))"))
    D = length(r)
    ranges = ntuple(d -> 0:min(r[d], ℓ[d]), D)
    return [collect(Int, s) for s in vec(collect(Iterators.product(ranges...)))]
end

enumerate_saturations(event::Event, ℓ::AbstractVector{<:Integer}) =
    enumerate_saturations(production_slots(event), ℓ)

"""
    safe_binomial(a::Integer, b::Integer) -> Int

`binomial(a, b)` under the KLI convention `C(a,b) = 0` if `b < 0` or
`b > a`.

Julia's `Base.binomial(a, b)` already returns `0` for `b < 0` and for
`0 <= a < b` (verified below, and by the `@assert`s that follow this
docstring at load time) -- but it does *not* return `0` for `a < 0`; it
instead falls back to the generalized/Pascal extension of the binomial
coefficient to negative upper arguments (e.g. `binomial(-1,1) == -1`),
which is the wrong convention for a combinatorial count. Since the KLI
formula's `a` arguments are population/lineage differences (`n_d - ℓ_d`)
that are only ever non-negative under the model invariant `ℓ_d <= n_d`,
`a < 0` should not arise from valid input -- but `kli_binomial_ratio` is
written defensively against it anyway, via this wrapper.
"""
function safe_binomial(a::Integer, b::Integer)
    (a < 0 || b < 0 || b > a) && return 0
    return binomial(a, b)
end

# Verify Julia's Base.binomial behavior on the boundary cases this module
# relies on, rather than assuming it (per the M02 task spec's explicit
# instruction). These run once at module load time.
@assert binomial(5, 0) == 1
@assert binomial(5, -1) == 0          # b < 0 -> 0 (Julia's own convention)
@assert binomial(5, 6) == 0           # b > a (a >= 0) -> 0 (Julia's own convention)
@assert binomial(0, 0) == 1
@assert binomial(0, 1) == 0
@assert binomial(-1, 1) == -1         # a < 0 -> Julia does NOT give 0 (Pascal extension);
                                       # this is exactly why safe_binomial exists.
@assert safe_binomial(-1, 1) == 0
@assert safe_binomial(5, -1) == 0
@assert safe_binomial(5, 6) == 0
@assert safe_binomial(0, 0) == 1

"""
    kli_binomial_ratio(r, s, ℓ, n; Q = 1) -> Rational{Int}
    kli_binomial_ratio(event::Event, s, ℓ, n; Q = 1) -> Rational{Int}

Compute the KLI production-slot compatibility ratio for a single
saturation `s`:

    φ_u(s) = Q · ∏_{d=1}^{D} C(n_d - ℓ_d, r_d - s_d) / C(n_d, r_d)

where `r` is the production vector (or `production_slots(event)`), `ℓ` the
current per-deme lineage (tracked-branch) counts, `n` the per-deme total
population counts, and `C(a,b)` the binomial coefficient under the
zero-outside-range convention implemented by `safe_binomial`. `Q` is the
external compatibility indicator/weight (KLI Eq.38's `Q_u`) -- this
milestone does not derive `Q_u`, it is supplied by the caller (default 1,
i.e. "assume compatible").

Returned as an exact `Rational{Int}` (not `Float64`), so boundary values
(0, 1, exact fractions) are exactly comparable in tests.

All four vector arguments must have equal length `D` (one entry per
lineage-carrying deme); mismatched lengths throw `ArgumentError`.

Convention for `n_d < r_d` (asking a deme to produce more individuals than
it has): `C(n_d, r_d) = 0` in that case (via `safe_binomial`), which would
make the ratio's denominator zero. Rather than throwing a `DivideError`,
this function treats that deme's factor -- and hence the whole product --
as exactly `0`: the saturation is unreachable/impossible given `n`. This is
a defensive convention documented in `handoffs/M02_kli_phi.md`; the
*model-level* invariant is that `n_d >= r_d` always holds for any reachable
post-event population state (the production already happened), so this
branch should not be exercised by valid M03+ callers, but the function
degrades gracefully rather than crashing if it is.

Convention for `s_d` outside its feasible range `0 <= s_d <= min(r_d,
ℓ_d)` (i.e. an `s` never actually produced by `enumerate_saturations`):
`s_d > r_d` is caught automatically (`r_d - s_d < 0` zeroes the numerator
via `safe_binomial`'s `b < 0` rule). `s_d < 0`, however, makes `r_d - s_d`
*larger* than `r_d`, which is not automatically zero -- so this function
explicitly checks `any(s_d < 0)` up front and returns `0` in that case too,
rather than silently computing a number for a nonsensical (negative-count)
saturation. This was found and fixed during M02's own test-writing (see
`handoffs/M02_kli_phi.md`), not assumed correct in advance.
"""
function kli_binomial_ratio(r::AbstractVector{<:Integer}, s::AbstractVector{<:Integer},
                             ℓ::AbstractVector{<:Integer}, n::AbstractVector{<:Integer};
                             Q::Real = 1)
    D = length(r)
    (length(s) == D && length(ℓ) == D && length(n) == D) ||
        throw(ArgumentError("r, s, ℓ, n must all have the same length (one entry " *
                             "per lineage-carrying deme); got lengths " *
                             "r=$(D), s=$(length(s)), ℓ=$(length(ℓ)), n=$(length(n))"))
    any(sd < 0 for sd in s) && return zero(Rational{Int})
    φ = Rational{Int}(Q)
    for d in 1:D
        denom = safe_binomial(n[d], r[d])
        if denom == 0
            return zero(Rational{Int})
        end
        numer = safe_binomial(n[d] - ℓ[d], r[d] - s[d])
        φ *= numer // denom
    end
    return φ
end

kli_binomial_ratio(event::Event, s::AbstractVector{<:Integer}, ℓ::AbstractVector{<:Integer},
                    n::AbstractVector{<:Integer}; Q::Real = 1) =
    kli_binomial_ratio(production_slots(event), s, ℓ, n; Q = Q)
