using EnumX: @enumx
using Random: AbstractRNG, default_rng

"""
    @marks Name mark1 mark2 ...

Creates a new module, `Name`, which contains an enumeration of the jump-marks `mark1`, `mark2`, ....
"""
macro marks(name, first, rest...)
    esc(:(@enumx $name $first=1 $(rest...)))
end

"""
    rcateg(p, prob = false; rng = Random.default_rng())

If `p` is a vector of weights, `rcateg(p)` returns a draw from the
categorical distribution on `1:length(p)` with these weights.
It also returns the sum of the weights.

If `prob=true`, the normalized weight of the selected category is returned as a third value.

It is not necessary for the weights to be normalized:
this is accomplished internally.

The random draw is made using `rng`, which defaults to the global RNG;
pass an explicit `AbstractRNG` for reproducible, independently-seeded draws.
"""
rcateg(
    p::AbstractVector{<:Real},
    prob::Bool = false;
    rng::AbstractRNG = default_rng(),
) = begin
    s::Prob = zero(Prob)
    for i ∈ eachindex(p)
        @assert p[i] ≥ 0 "invalid p[$i]=$(p[i]) detected"
        s += Prob(p[i])
    end
    if s > 0
        r = s*rand(rng, Prob)
        k::Int = 1
        n = lastindex(p)
        ## Subtraction can leave a remainder past the last weight.
        ## That remainder belongs to the last category with positive weight.
        while k < n && r > p[k]
            r -= p[k]
            k += 1
        end
        if p[k] ≤ 0
            ## k == n: the remainder ran past a zero last weight.
            ## k < n: the draw was 0 and the first weight is 0.
            ## s > 0, so some category has weight.
            k = (k == n) ? findlast(>(0), p) : findfirst(>(0), p)
        end
        if prob
            k, s, p[k]/s
        else
            k, s
        end
    else
        if prob
            zero(Size), zero(Prob), zero(Prob)
        else
            zero(Size), zero(Prob)
        end
    end
end

"""
    rcateg(p, e, prob = false)

This call returns a draw from the categorical distribution on the enumeration `e`.  `e` should be an enumeration-type created, e.g., by `@demes`.
"""
rcateg(
    p::AbstractVector{<:Real},
    demes::Type{D},
    prob::Bool = false;
    rng::AbstractRNG = default_rng(),
) where {D <: Enum} = begin
    k, s... = rcateg(p, prob; rng)
    d = (k > 0) ? demes(k) : missing
    d, s...
end

"""
    rcateg(p, s, prob = false)

This call returns a draw from the categorical distribution on the set `s`.
The latter may be a `Set` or a `BitSet`.
"""
rcateg(
    p::AbstractVector{<:Real},
    set::Union{BitSet, Set},
    prob::Bool = false;
    rng::AbstractRNG = default_rng(),
) = begin
    k, s... = rcateg(p, prob; rng)
    if k > 0
        v = collect(set)::Vector{Int}
        vk = v[k]
    else
        vk = missing
    end
    vk, s...
end

"""
    rcateg(p, v, prob = false)

This call returns a draw from the categorical distribution on the vector `v`.
"""
rcateg(
    p::AbstractVector{<:Real},
    v::AbstractVector,
    prob::Bool = false;
    rng::AbstractRNG = default_rng(),
) = begin
    k, s... = rcateg(p, prob; rng)
    vk = (k > 0) ? v[k] : missing
    vk, s...
end
