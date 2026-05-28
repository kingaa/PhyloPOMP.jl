"""
    @marks Name mark1 mark2 ...

Creates a new module, `Name`, which contains an enumeration of the jump-marks `mark1`, `mark2`, ....
"""
macro marks(name, first, rest...)
    expr = :(@enumx $name $first=1 $(rest...))
    esc(:@eval $expr)
end

rcateg(
    p::AbstractVector{<:Real},
) = begin
    s = zero(Float64)
    for i ∈ eachindex(p)
        @assert p[i]>=0 "invalid p[$i]=$(p[i]) detected"
        s += p[i]
    end
    r = s*rand()
    k::Int = 1
    while r > p[k]
        r -= p[k]
        k += 1
    end
    k, s
end

rcateg(
    p::AbstractVector{<:Real},
    types::Module,
) = begin
    k,s = rcateg(p)
    types.T(k), s
end
