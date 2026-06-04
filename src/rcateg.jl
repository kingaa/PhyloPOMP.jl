"""
    @marks Name mark1 mark2 ...

Creates a new module, `Name`, which contains an enumeration of the jump-marks `mark1`, `mark2`, ....
"""
macro marks(name, first, rest...)
    expr = :(@enumx $name $first=1 $(rest...))
    esc(:@eval $expr)
end

"""
    rcateg(p, prob = false)

If `p` is a vector of weights, `rcateg(p)` returns a draw from the
categorical distribution on `1:length(p)` with these weights.
It also returns the sum of the weights.

If `prob=true`, the normalized weight of the selected category is returned as a third value.

It is not necessary for the weights to be normalized:
this is accomplished internally.
"""
rcateg(
    p::AbstractVector{<:Real},
    prob::Bool = false,
) = begin
    s::Prob = zero(Prob)
    for i ∈ eachindex(p)
        @assert p[i] ≥ 0 "invalid p[$i]=$(p[i]) detected"
        s += Prob(p[i])
    end
    ## NB: this is a non-destructive routine.
    ## If it is permissible to destroy p, one could save the
    ## extra round of subtractions.
    if s > 0
        r = s*rand(Prob)
        k::Int = 1
        while r > p[k]
            r -= p[k]
            k += 1
        end
        if prob
            k, s, p[k]/s
        else
            k, s
        end
    else
        if prob
            one(Int), zero(Prob), zero(Prob)
        else
            one(Int), zero(Prob)
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
    prob::Bool = false,
) where {D <: Enum} = begin
    k,s... = rcateg(p,prob)
    demes(k), s...
end
