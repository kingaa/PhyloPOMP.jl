import Base: eachindex, length, getindex, eachindex

"""
    GuideNode{F,D}

A `GuideNode` contains the information needed to guide a phylodynamic Monte Carlo filter on a particular interval between genealogical events.
`F` is an `AbstractFloat` type and `D` is the enumeration of the demes.

A `GuideNode` contains:
- the time interval to which it pertains
- the type of event (`NodeType`) at the left endpoint of the interval
- the name of the parent lineage (ancestral to the genealogical node at the left endpoint)
- the names of the child lineages (immediately descended from the parent node)
- the names of all lineages present on the interval
- the guide probability for each child lineage at the left interval endpoint
- the target probabilities for all lineages at the right interval endpoint.
"""
struct GuideNode{F<:AbstractFloat,D<:Enum}
    type::NodeType
    tbeg::Time
    tend::Time
    parlin::Name
    chillins::Vector{Name}
    alllins::Vector{Name}
    present::Matrix{F}
    target::Matrix{F}
end

"""
    Guide{F,D}

This type holds a "filter-guide" based on a reverse-time finite-state Markov process.
Construct it with a call to `guide`.
`F` is an `AbstractFloat` type and `D` is the enumeration of the demes.
It contains the `FSMarkovProc` object encoding the Markov process and a sequence of `GuideNode` objects.
"""
struct Guide{F<:AbstractFloat,D<:Enum}
    process::FSMarkovProc{F,D}
    nodes::Vector{GuideNode{F,D}}
    Guide(
        g::Genealogy{D},
        m::FSMarkovProc{F,D},
    ) where {F<:AbstractFloat,D<:Enum} = begin
        probs = zeros(F,length(instances(D)),g.nsample)
        nodes = Vector{GuideNode{F,D}}(undef,length(g))
        tend::Time = g.time
        lins::Set{Name} = Set{Name}()
        for n ∈ reverse(eachindex(g))
            ells = collect(lins)
            parlin = g[n].lineage
            chillins = map(x->g[x].lineage, g[n].children)
            target = probs[:,ells]
            probs[:,ells] = forward_action(m,tend-g[n].slate,probs[:,ells])
            present = probs[:,chillins]
            nodes[n] = GuideNode{F,D}(
                g[n].type,
                g[n].slate,tend,
                parlin,chillins,ells,
                present,target,
            )
            if ismissing(g[n].deme)
                @inbounds assimil!(
                    @view(probs[:,parlin]),
                    @view(probs[:,chillins]),
                    statdist(m)
                )
            else
                probs[:,parlin] = demekron(F,g[n].deme)
            end
            tend = g[n].slate
            union!(setdiff!(lins,chillins),parlin)
        end
        new{F,D}(m,nodes)
    end
end

"""
    guide(g, m)

- `g` is a `Genealogy`.
- `m` is a finite-state Markov process (type `FSMarkovProc`) constructed via a call to `fsmarkov`.
"""
guide(args...) = Guide(args...)

"""
    demekron([F::Type], d)

Kronecker delta for deme `d`.
"""
demekron(F::Type, d::D,) where {D <: Enum} = begin
    [i == d ? one(F) : zero(F) for i ∈ instances(D)]
end

"""
    relhaz(g, t, n)

Compute the relative reverse-time hazards associated with filter guide `g` at time `t` relative to the target at guide-node `n`.
"""
relhaz(
    g::Guide{F,D},
    t::Real,
    d::D,
    n::Integer,
) where {F,D} = begin
    if (t > g[n].tend)
        error("improper `relhaz` call: cannot evaluate at t > $(g[n].tend)")
    end
    j = Int(d)
    p = forward_action(g.process,g[n].tend-t,g[n].target) ./ statdist(g.process)
    Dict((g[n].alllins[k],D(i)) => p[i,k]/p[j,k]
         for k ∈ eachindex(g[n].alllins), i ∈ axes(p,1) if i != j)
end

assimil!(
    parprob::AbstractVector{F},
    chilprobs::AbstractMatrix{F},
    pi::AbstractVector{F},
) where {F<:AbstractFloat} = begin
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
