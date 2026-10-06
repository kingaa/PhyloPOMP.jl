# Saturations S_u(ℓ) and the production-slot binomial ratio φ_u.

export production_slots, enumerate_saturations, kli_binomial_ratio

"""
    production_slots(event::Event) -> Vector{Int}

Return `event`'s production vector `r_u` (a per-mark constant), one entry
per lineage-carrying deme, in `MGPModel.demes` order.
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

`length(r) == length(ℓ)` is required (one entry per lineage-carrying deme);
mismatched lengths throw `ArgumentError`.

`ℓ_d = 0` or `r_d = 0` collapses deme `d`'s range to `{0}`.
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

`binomial(a, b)` under the convention `C(a,b) = 0` if `b < 0` or `b > a`.

`Base.binomial` extends to negative `a` (`binomial(-1, 1) == -1`).
This returns 0 there, as for `b < 0` and `b > a`.
"""
function safe_binomial(a::Integer, b::Integer)
    (a < 0 || b < 0 || b > a) && return 0
    return binomial(a, b)
end

"""
    kli_binomial_ratio(r, s, ℓ, n; Q = 1) -> Rational{Int}
    kli_binomial_ratio(event::Event, s, ℓ, n; Q = 1) -> Rational{Int}

Compute the production-slot compatibility ratio for a single
saturation `s`:

    φ_u(s) = Q · ∏_{d=1}^{D} C(n_d - ℓ_d, r_d - s_d) / C(n_d, r_d)

where `r` is the production vector (or `production_slots(event)`), `ℓ` the
current per-deme lineage (tracked-branch) counts, `n` the per-deme total
population counts, and `C(a,b)` the binomial coefficient under the
zero-outside-range convention implemented by `safe_binomial`.
`Q` is the compatibility indicator, supplied by the caller.

Returned as an exact `Rational{Int}` (not `Float64`), so boundary values
(0, 1, exact fractions) are exactly comparable in tests.

Returns 0 if `n_d < r_d` for some `d` or if any `s_d < 0`.
Throws `ArgumentError` if the vectors differ in length.
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
