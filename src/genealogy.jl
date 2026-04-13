import Base: eachindex, length, getindex, eachindex

const Name = Int64
const Time = Float64
## FIXME: inclusion of Root in NodeType introduces some inelegant redundancy
@enum NodeType Root Sample Node

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
    lineage::Union{Missing,Name}
    parent::Union{Nothing,Name}
    children::Vector{Name}
    GenealNode{D}(
        name::Name,
        slate::Time,
        deme = missing,
        type::NodeType = Node,
        parent::Union{Nothing,Name} = nothing,
    ) where {D <: Enum} = begin
        new{D}(type,name,slate,deme,missing,parent,Name[])
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

roots(g::Genealogy) = findall(i->(g[i].type==Root),eachindex(g))
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
    ## repair names and types:
    namemap = Dict(n.name=>k for (k,n) ∈ enumerate(G.nodes))
    for n ∈ G.nodes
        n.name = namemap[n.name]
        if isnothing(n.parent)
            n.type = Root
        else
            n.parent = namemap[n.parent]
        end
        n.children = map(i->namemap[i],n.children)
    end
    trace_lineages!(G)
    G
end

trace_lineages!(G::Genealogy, p::Integer, i::Integer) = begin
    G[p].lineage = i
    q = G[p].parent
    if !isnothing(q) && ismissing(G[q].lineage)
        trace_lineages!(G,q,i)
    end
    nothing
end

trace_lineages!(G::Genealogy) = begin
    i::Int64 = 1
    for p ∈ eachindex(G)
        if G[p].type==Sample
            trace_lineages!(G,p,i)
            i = i+1
        end
    end
    nothing
end
