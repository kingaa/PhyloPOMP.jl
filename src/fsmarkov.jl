"""
    FSMarkov

This module contains facilities for constructing and working with finite-state Markov processes in continuous time.
"""
module FSMarkov

export fsmarkov, statdist, transition, relprob, relhaz, demekron

import LinearAlgebra: Diagonal, Symmetric, Transpose, eigen

"""
    FSMarkovProc

Represents a finite-state Markov process in continuous time.
"""
struct FSMarkovProc{F<:AbstractFloat,D<:Enum}
    statdist::Vector{F}
    generator::Matrix{F}
    eigenvals::Vector{F}
    left_trans::Matrix{F}
    right_trans::Matrix{F}
    FSMarkovProc(
        args::Union{Pair{D,<:Real},Pair{Tuple{D,D},<:Real}}...,
    ) where {D <: Enum} = FSMarkovProc(Float64,args...)
    FSMarkovProc(
        F::Type{<:AbstractFloat},
        args::Union{Pair{D,<:Real},Pair{Tuple{D,D},<:Real}}...,
    ) where {D <: Enum} = begin
        idx = Int.(instances(D))
        n = length(idx)
        pi = zeros(F,n)
        Q = zeros(F,n,n)
        for r ∈ args
            if isa(r,Pair{D,<:Real})
                s = F(r.second)
                if s ≤ 0.0
                    error("non-positive stationary probability $r")
                end
                pi[Int(r.first)] = s
            else
                if r.first[1] == r.first[2]
                    error("cannot specify conductance of a deme to itself (here, $(r.first[1]) to $(r.first[2]))")
                elseif Q[Int(r.first[2]),Int(r.first[1])] ≠ zero(F)
                    error("double specification of conductance ($(r.first[1]) to $(r.first[2]))")
                else
                    s = F(r.second)
                    if (s < zero(F))
                        error("negative conductance for $(r.first[1]) to $(r.first[2])")
                    end
                    Q[Int(r.first[1]),Int(r.first[2])] = s
                    Q[Int(r.first[2]),Int(r.first[1])] = s
                end
            end
        end
        for i ∈ idx
            if pi[i] ≤ zero(F)
                error("unspecified stationary probability for $(D(i))")
            end
        end
        pi = pi./sum(pi)        # normalize stationary distribution
        Q = Diagonal(pi)*Q
        for i ∈ idx
            Q[i,i] = -sum(Q[:,i])
        end
        s = sqrt.(pi)
        d = Diagonal(s)
        di = Diagonal(one(F)./s)
        Lambda,U = eigen(Symmetric(di*Q*d),sortby=x->abs(x))
        Lambda[1] = zero(F)
        new{F,D}(pi,Q,Lambda,d*U,Transpose(U)*di)
    end
end

"""
    fsmarkov(...)

Constructs a finite-state Markov process in continuous time, represented as an object of type `FSMarkovProc`.
Each argument is of the form `D => prob` or
`(D1,D2) => cond`, where `D`, `D1`, `D2` are demes, `prob` is a probability, and `cond` is a conductance.
Specifically,
- `D => prob` indicates that the stationary distribution on deme `D` is proportional to `prob`.
- `(D1,D2) => cond` mandates that the conductance between `D1` and `D2` is `cond`.
"""
fsmarkov(args...,) = FSMarkovProc(args...)

statdist(m::FSMarkovProc,) = m.statdist

transition(
    m::FSMarkovProc{F}, s::Real,
) where F = begin
    m.left_trans*(Diagonal(exp.(m.eigenvals.*F(s)))*m.right_trans)
end

transition(
    m::FSMarkovProc{F},
    s::Real,
    X::AbstractArray{<:Real,N},
) where {F,N} = begin
    n = size(X,1)
    if n != length(statdist(m))
        error("size mismatch in 'transition'")
    end
    m.left_trans*(Diagonal(exp.(m.eigenvals.*F(s)))*(m.right_trans*X))
end

relprob(
    m::FSMarkovProc,
    s::Real,
    q::AbstractVector{<:Real},
) = begin
    transition(m,s,q) ./ statdist(m)
end

relhaz(
    m::FSMarkovProc{F}, s::Real, q::AbstractVector{<:Real},
) where F = begin
    Q = m.generator
    Pback = relprob(m,s,q)
    [Q[i,j] > 0 ? Pback[i]/Pback[j] : one(F)
     for i ∈ axes(Q,1), j ∈ axes(Q,2)]
end

demekron(F::Type{<:Real}, d::D,) where {D <: Enum} = begin
    [ i == d ? one(F) : zero(F) for i ∈ instances(D) ]
end

demekron(d::Enum,) = demekron(Float64,d)

end
