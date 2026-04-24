import Base: eachindex, length, getindex, eachindex

struct GuideNode{F<:Real,D<:Enum}
    tbeg::Time
    tend::Time
    lineages::Vector{Name}
    probs::Matrix{F}
end

struct Guide{F<:Real,D<:Enum}
    process::FSMarkovProc{F,D}
    nodes::Vector{GuideNode{F,D}}
    Guide(
        g::Genealogy{D},
        m::FSMarkovProc{F,D},
    ) where {F<:Real,D<:Enum} = begin
        ndeme = length(instances(D))
        probs = zeros(F,ndeme,g.nsample)
        nodes = Vector{GuideNode{F,D}}(undef,length(g))
        tend::Time = g.time
        lins::Set{Name} = Set{Name}()
        for n ∈ reverse(eachindex(g))
            ells = collect(lins)
            nodes[n] = GuideNode{F,D}(g[n].slate,tend,ells,probs[:,ells])
            probs[:,ells] = transition(m,tend-g[n].slate,probs[:,ells])
            parlin = g[n].lineage
            childlins = map(x->g[x].lineage, g[n].children)
            if ismissing(g[n].deme)
                @inbounds assimil!(
                    @view(probs[:,parlin]),
                    @view(probs[:,childlins]),
                    statdist(m)
                )
            else
                probs[:,parlin] = demekron(F,g[n].deme)
            end
            tend = g[n].slate
            union!(setdiff!(lins,childlins),parlin)
        end
        new{F,D}(m,nodes)
    end
end

demekron(F::Type, d::D,) where {D <: Enum} = begin
    [i == d ? one(F) : zero(F) for i ∈ instances(D)]
end

relhaz(
    g::Guide{F,D},
    t::Real,
    n::Integer,
) where {F,D} = begin
    p = transition(g.process,t,g[n].probs) ./ statdist(g.process)
    Dict((g[n].lineages[k],D(j)=>D(i)) => p[i,k]/p[j,k]
         for k ∈ eachindex(g[n].lineages), i ∈ axes(p,1), j ∈ axes(p,1) if i != j)
end

assimil!(
    parprob::AbstractVector{F},
    chilprobs::AbstractMatrix{F},
    pi::AbstractVector{F},
) where {F<:Real} = begin
    @assert size(parprob,1)==size(chilprobs,1)
    if isempty(chilprobs)
        parprob[:] = pi
    else
        for i ∈ axes(chilprobs,1)
            for j ∈ axes(chilprobs,2)
                @inbounds chilprobs[i,j] /= pi[i]
            end
        end
        s = zero(F)
        for i ∈ axes(chilprobs,1)
            @inbounds parprob[i] = pi[i]*prod(chilprobs[i,:])
            @inbounds s += parprob[i]
        end
        for i ∈ eachindex(parprob)
            @inbounds parprob[i] /= s
        end
    end
    nothing
end

Base.getindex(g::Guide, i::Integer) = g.nodes[i]
Base.getindex(g::Guide, i::AbstractVector{<:Integer}) = g.nodes[i]
Base.length(g::Guide) = length(g.nodes)
Base.eachindex(g::Guide) = eachindex(g.nodes)
