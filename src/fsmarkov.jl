"""
    FSMarkov

Module implementing finite-state Markov processes for use in filter guides.
"""
module FSMarkov

export fsmarkov, FSMarkovProc, forward_action, statdist, generator

import Base: show
import LinearAlgebra: Diagonal, Symmetric, Transpose, eigen
import ..PhyloPOMP: Prob

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
    ) where {D <: Enum} = FSMarkovProc(Prob,args...)
    FSMarkovProc(
        F::Type{<:AbstractFloat},
        args::Union{Pair{D,<:Real},Pair{Tuple{D,D},<:Real}}...,
    ) where {D <: Enum} = begin
        pi,Q = make_generator(F,args...)
        s = sqrt.(pi)
        d = Diagonal(s)
        di = Diagonal(one(F)./s)
        Lambda,U = eigen(Symmetric(di*Q*d),sortby=-)
        Lambda[1] = zero(F)
        new{F,D}(pi,Q,Lambda,d*U,Transpose(U)*di)
    end
end

"""
    fsmarkov(...)

Constructs a finite-state Markov process in continuous time, represented as an object of type [`FSMarkovProc`](@ref).
Each argument is of the form `D => prob` or
`(D1,D2) => cond`, where `D`, `D1`, `D2` are demes, `prob` is a probability, and `cond` is a conductance.
Specifically,
- `D => prob` indicates that the stationary distribution on deme `D` is proportional to `prob`.
- `(D1,D2) => cond` mandates that the conductance between `D1` and `D2` is `cond`.
"""
fsmarkov(args...,) = FSMarkovProc(args...)

statdist(m::FSMarkovProc,) = m.statdist

generator(m::FSMarkovProc,) = m.generator

forward_action(
    m::FSMarkovProc{F},
    s::Real,
) where F = begin
    d = exp.(m.eigenvals.*F(s))
    m.left_trans*(Diagonal(d)*m.right_trans)
end

forward_action(
    m::FSMarkovProc{F},
    s::Real,
    X::AbstractArray{<:Real,N},
) where {F,N} = begin
    n = size(X,1)
    if n != length(statdist(m))
        error("size mismatch in 'forward_action'")
    end
    d = exp.(m.eigenvals.*F(s))
    m.left_trans*(Diagonal(d)*(m.right_trans*X))
end

make_generator(
    F::Type{<:AbstractFloat},
    args::Union{Pair{D,<:Real},Pair{Tuple{D,D},<:Real}}...,
) where {D <: Enum} = begin
    idx = Int.(instances(D))
    n = length(idx)
    pi = zeros(F,n)
    Q = zeros(F,n,n)
    for arg ∈ args
        if arg isa Pair{D,<:Real}
            s = F(arg.second)
            if s ≤ 0.0
                error("non-positive stationary probability $arg")
            end
            pi[Int(arg.first)] = s
        else
            if arg.first[1] == arg.first[2]
                error("cannot specify conductance of a deme to itself (here, $(arg.first[1]) to $(arg.first[2]))")
            elseif Q[Int(arg.first[2]),Int(arg.first[1])] ≠ zero(F)
                error("double specification of conductance ($(arg.first[1]) to $(arg.first[2]))")
            else
                s = F(arg.second)
                if (s < zero(F))
                    error("negative conductance for $(arg.first[1]) to $(arg.first[2])")
                end
                Q[Int(arg.first[1]),Int(arg.first[2])] = s
                Q[Int(arg.first[2]),Int(arg.first[1])] = s
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
    pi, Q
end

Base.show(
    io::IO,m::FSMarkovProc;kwargs...,
) = print(io,pretty_string(m;kwargs...))

pretty_string(m::FSMarkovProc,) = begin
    Q = generator(m)
    "<FSMarkovProc with generator = $Q>"
end

end
