export Genealogy, GenealNode

const Name = Int64
const Time = Union{Float64,Int64}
@enum NodeType Sample Node

"""
    GenealNode{D,T}

Implements a *genealogical node*.
The `Enum` type `D` enumerates the demes, and `T` is the type of the time.
"""
mutable struct GenealNode{D<:Enum,T<:Time}
    type::NodeType
    name::Name
    slate::T
    deme::Union{Missing,D}
    parent::Union{Nothing,Name}
    children::Vector{Name}
    GenealNode{D}(
        name::Name,
        slate::T,
        deme = missing,
        type::NodeType = Node
    ) where {D <: Enum, T <: Time} = begin
        new{D,T}(type,name,slate,deme,nothing,Name[])
    end
end

"""
    Genealogy{D,T}

Implements a genealogy over a time interval `[t0,t]`.
The `Enum`-type `D` enumerates the demes and `T` is the type of the time variable.
Internally, this is represented as a time-ordered sequence of *genealogical nodes* (represented by [`GenealNode`](@ref) objects).
"""
mutable struct Genealogy{D <: Enum, T <: Time}
    t0::T
    time::T
    safe::Bool
    nodes::Vector{GenealNode{D,T}}
    Genealogy{D}(t0::T) where {D <: Enum, T <: Time} = begin
        new{D,T}(t0,t0,true,GenealNode{D,T}[])
    end
end

"""
    weed!(G)

Removes all dead roots.
The genealogy becomes incorrect and needs to be repaired (see [`repair!`](@ref)).
"""
weed!(G::Genealogy) = begin
    filter!(
        n -> (!isnothing(n.parent) || !isempty(n.children)),
        G.nodes
    )
    G.safe = false
    nothing
end

"""
   repair!(G)

Re-sorts the genealogical nodes.
Renames all nodes and corrects the parent/child relationships.
"""
repair!(G::Genealogy{D,T}) where {D,T} = begin
    if !G.safe
        compare(p,q) = p.slate < q.slate ||
            (p.slate == q.slate &&
            (p==q.parent || (q!=p.parent && p.name < q.name)))
        sort!(G.nodes,lt=compare)
        namemap = Dict((n.name,k) for (k,n) ∈ enumerate(G.nodes))
        for n ∈ G.nodes
            n.name = namemap[n.name]
            n.parent = get(namemap,n.parent,nothing)
            n.children = map(x->namemap[x],n.children)
        end
        G.safe = true
    end
    G
end
