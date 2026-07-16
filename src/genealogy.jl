import Base: eachindex, length, getindex, eachindex, ==, show
import PartiallyObservedMarkovProcesses: times, timezero

"""
    NodeType

Genealogical nodes (see [`GenealNode`](@ref)) can be one of three types: Root, Node, or Sample.
"""
@enum NodeType Root Sample Node
## FIXME: inclusion of Root in NodeType introduces some inelegant redundancy

"""
    GenealNode{E}

Implements a *genealogical node*.
The type `E` is an enumeration of the demes (see [`@demes`](@ref)).
"""
mutable struct GenealNode{E<:Enum}
    "Type of node: Root, Sample, or Node."
    type::NodeType
    "Unique node name."
    name::Name
    "Node-time."
    slate::Time
    "Deme of the branch ascendant from the node."
    deme::Union{Missing,E}
    "Lineage name."
    lineage::Union{Missing,Name}
    "Name of the parent node."
    parent::Union{Nothing,Name}
    "Names of child nodes."
    children::Vector{Name}
    GenealNode{E}(
        name::Integer,
        slate::Time,
        deme = missing,
        type::NodeType = Node,
        parent::Union{Nothing,Name} = nothing,
    ) where {E <: Enum} = begin
        new{E}(type,Name(name),slate,deme,missing,parent,Name[])
    end
end

"""
    Genealogy{D}

Implements a genealogy over a time interval `[t0,t]`.
The module `D` enumerates the demes and should be constructed
via a call to [`@demes`](@ref).
Internally, a genealogy is represented as a time-ordered sequence
of *genealogical nodes* (represented by [`GenealNode`](@ref) objects).
"""
mutable struct Genealogy{D}
    "Root-time of genealogy."
    t0::Time
    "Time of genealogy, i.e., right endpoint of the interval."
    time::Time
    "Number of samples."
    nsample::Size
    "Vector of genealogical nodes."
    nodes::Vector{GenealNode}
    "Constructor of an empty genealogy on the interval [t0,time]."
    Genealogy{D}(t0::Real, time::Real = t0,) where D = begin
        @isademeset D
        new{D}(Time(t0),Time(time),zero(Size),GenealNode{D.DemeSet}[])
    end
end

Base.getindex(g::Genealogy, i::Integer) = g.nodes[i]
Base.getindex(g::Genealogy, i::AbstractVector{<:Integer}) = g.nodes[i]
Base.length(g::Genealogy) = length(g.nodes)
Base.eachindex(g::Genealogy) = eachindex(g.nodes)

roots(g::Genealogy) = findall(i->(g[i].type==Root),eachindex(g))
tips(g::Genealogy) = findall(i->isempty(g[i].children),eachindex(g))
samples(g::Genealogy) = findall(i->(g[i].type==Sample),eachindex(g))
nodes(g::Genealogy) = findall(i->(g[i].type==Node),eachindex(g))
timezero(g::Genealogy) = g.t0
times(g::Genealogy) = [map(n->n.slate,g.nodes[2:end])...,g.time]
nsample(g::Genealogy) = g.nsample

==(n1::GenealNode{E1}, n2::GenealNode{E2}) where {E1,E2} = begin
    E1 == E2 &&
        n1.type == n2.type &&
        n1.name == n2.name &&
        n1.slate == n2.slate &&
        n1.deme === n2.deme &&
        n1.lineage == n2.lineage &&
        n1.parent == n2.parent &&
        n1.children == n2.children
end

==(g1::Genealogy{D1}, g2::Genealogy{D2}) where {D1,D2} = begin
    D1==D2 &&
        length(g1)==length(g2) &&
        g1.t0 == g2.t0 &&
        g1.time == g2.time &&
        all([g1[i]==g2[i] for i ∈ eachindex(g1)])
end

"""
    repair!(G)

Re-sorts the genealogical nodes.
Renames all nodes and corrects the parent/child relationships.
Traces sample-lineages.
"""
repair!(G::Genealogy) = begin
    ## weed out dead roots:
    filter!(
        n -> !(isnothing(n.parent) && isempty(n.children)),
        G.nodes,
    )
    ## sort nodes:
    compare(p,q) = p.slate < q.slate ||
        (p.slate == q.slate &&
        (p==q.parent || (q!=p.parent && p.name < q.name)))
    sort!(G.nodes,lt=compare)
    ## repair names and types:
    namemap = Dict(n.name=>k for (k,n) ∈ enumerate(G.nodes))
    G.nsample = zero(Size)
    for n ∈ G.nodes
        n.name = namemap[n.name]
        if isnothing(n.parent)
            n.type = Root
        else
            n.parent = namemap[n.parent]
        end
        n.children = sort(map(i->namemap[i],n.children))
        if n.type==Sample
            G.nsample += 1
        end
    end
    trace_lineages!(G)
    nothing
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

Base.show(
    io::IO,
    g::Union{Genealogy,GenealNode,AbstractVector{GenealNode}};
    kwargs...,
) = begin
    print(io, pretty_string(g; kwargs...))
end

pretty_string(g::Genealogy; sigdigits=4) = begin
    "<genealogy on ["*
        "$(round(g.t0;sigdigits)),"*
        "$(round(g.time;sigdigits))]:\n" *
        pretty_string(g.nodes; sigdigits) *
        ">"
end

pretty_string(g::AbstractVector{GenealNode}; sigdigits=4) = begin
    join(
        map(eachindex(g)) do i
            "  $(g[i].name): " *
                pretty_string(g[i]; sigdigits)
        end,
        '\n'
    )
end

pretty_string(p::GenealNode; sigdigits = 4) = begin
    parent_string = if isnothing(p.parent)
        ""
    else
        "parent=$(p.parent) "
    end
    children_string = if isempty(p.children)
        ""
    else
        join(map(i->"$i",p.children),',')
    end
    "<$(lowercase(String(Symbol(p.type)))) " *
        "lineage=$(p.lineage) " *
        (ismissing(p.deme) ? "" : "deme=$(p.deme) ") *
        "time=$(round(p.slate;sigdigits)) " *
        parent_string *
        "children=[$children_string]" * ">"
end
