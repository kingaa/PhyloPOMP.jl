export Genealogy, GenealNode

const Name = Int64
const Time = Float64
@enum NodeType Sample Node

"""
    GenealNode{D}

Implements a *genealogical node*.
The `Enum` type `D` enumerates the demes.
"""
mutable struct GenealNode{D<:Enum}
    type::NodeType
    name::Name
    slate::Time
    deme::Union{Missing,D}
    parent::Union{Nothing,Name}
    children::Vector{Name}
    GenealNode{D}(
        name::Name,
        slate::Time,
        deme = missing,
        type::NodeType = Node,
        parent::Union{Nothing,Name} = nothing,
    ) where {D <: Enum} = begin
        new{D}(type,name,slate,deme,parent,Name[])
    end
end

"""
    Genealogy{D}

Implements a genealogy over a time interval `[t0,t]`.
The `Enum`-type `D` enumerates the demes.
Internally, this is represented as a time-ordered sequence of *genealogical nodes* (represented by [`GenealNode`](@ref) objects).
"""
mutable struct Genealogy{D <: Enum}
    t0::Time
    time::Time
    nodes::Vector{GenealNode{D}}
    Genealogy{D}(t0::Real) where {D <: Enum} = begin
        new{D}(Time(t0),Time(t0),GenealNode{D}[])
    end
end

Base.getindex(g::Genealogy, i::Integer) = g.nodes[i]
Base.getindex(g::Genealogy, i::Vector{<:Integer}) = g.nodes[i]
Base.length(g::Genealogy) = length(g.nodes)
Base.eachindex(g::Genealogy) = eachindex(g.nodes)

roots(g::Genealogy) = findall(i->isnothing(g[i].parent),eachindex(g))
tips(g::Genealogy) = findall(i->isempty(g[i].children),eachindex(g))
samples(g::Genealogy) = findall(i->(g[i].type==Sample),eachindex(g))
nodes(g::Genealogy) = findall(i->(g[i].type==Node),eachindex(g))

"""
    repair!(G)

Re-sorts the genealogical nodes.
Renames all nodes and corrects the parent/child relationships.
"""
repair!(G::Genealogy{D}) where {D} = begin
    ## weed out dead roots:
    filter!(
        n -> !(isnothing(n.parent) && isempty(n.children)),
        G.nodes
    )
    ## sort nodes:
    compare(p,q) = p.slate < q.slate ||
        (p.slate == q.slate &&
        (p==q.parent || (q!=p.parent && p.name < q.name)))
    sort!(G.nodes,lt=compare)
    ## repair names:
    namemap = Dict(n.name=>k for (k,n) ∈ enumerate(G.nodes))
    for n ∈ G.nodes
        n.name = namemap[n.name]
        if !isnothing(n.parent)
            n.parent = namemap[n.parent]
        end
        n.children = map(i->namemap[i],n.children)
    end
    G
end
