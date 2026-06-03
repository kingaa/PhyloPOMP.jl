import Base: eachindex, length, getindex, eachindex, show
import PartiallyObservedMarkovProcesses: times, timezero

"""
    NodeType

Genealogical nodes (see [`GenealNode`](@ref)) can be one of three types: Root, Node, or Sample.
"""
@enum NodeType Root Sample Node
## FIXME: inclusion of Root in NodeType introduces some inelegant redundancy

"""
    GenealNode{D}

Implements a *genealogical node*.
The module `D` enumerates the demes (see [`@demes`](@ref)).
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
        name::Integer,
        slate::Time,
        deme = missing,
        type::NodeType = Node,
        parent::Union{Nothing,Name} = nothing,
    ) where {D <: Enum} = begin
        new{D}(type,Name(name),slate,deme,missing,parent,Name[])
    end
end

"""
    Genealogy{D}

Implements a genealogy over a time interval `[t0,t]`.
The module `D` enumerates the demes (see [`@demes`](@ref)).
Internally, this is represented as a time-ordered sequence of *genealogical nodes* (represented by [`GenealNode`](@ref) objects).
"""
mutable struct Genealogy{D <: Enum}
    t0::Time
    time::Time
    nsample::Size
    nodes::Vector{GenealNode{D}}
    Genealogy{D}(t0::Real, time::Real = t0) where {D <: Enum} = begin
        new{D}(Time(t0),Time(time),zero(Size),GenealNode{D}[])
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
        n.children = map(i->namemap[i],n.children)
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


"""
    geneal_convert(g, f, D)

Convert the genealogy `g` to a `Genealogy{D}` by applying the function `f` to each of its nodes.  `f` must have signature `f(type, time, deme)`, where `type` will be the node type, `time` its time, and `deme` its deme.  `f` must return a deme. i.e., an instance of `D`, or `missing`.
"""
geneal_convert(
    g::Genealogy,
    f::Function,
    D::Type{E},
) where {E <: Enum} = begin
    h = Genealogy{D}(g.t0,g.time)
    for m in g.nodes
        d = f(m.type,m.slate,m.deme)::Union{D,Missing}
        n = GenealNode{D}(m.name,m.slate,d,m.type,m.parent)
        n.children = copy(m.children)
        push!(h.nodes,n)
    end
    repair!(h)
    h
end

Base.show(
    io::IO,
    g::Union{Genealogy,GenealNode};
    kwargs...,
) = begin
    print(io,pretty_string(g;kwargs...))
end

pretty_string(g::Genealogy; sigdigits=4) = begin
    "<genealogy on ["*
        "$(round(g.t0,sigdigits=sigdigits)),"*
        "$(round(g.time,sigdigits=sigdigits))]:\n" *
        join(
            map(eachindex(g)) do i
                "  $i: " * pretty_string(g[i],sigdigits=sigdigits)
            end,
            '\n'
        ) * ">"
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
        "deme=$(p.deme) " *
        "time=$(round(p.slate,sigdigits=sigdigits)) " *
        parent_string *
        "children=[$children_string]" * ">"
end
