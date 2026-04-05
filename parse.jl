using AbstractTrees

@enum NodeType begin
    Sample
    Node
end

const Name = Int64
const Time = Union{Float64,Int64}

enum_name(name::E) where {E <: Enum} = lowercase(String(Symbol(name)))

name2enum(::Val{E}) where {E <: Enum} = begin
    i = instances(E)
    Dict(zip(enum_name.(i),i))
end

enum2name(::Val{E}) where {E <: Enum} = begin
    i = instances(E)
    Dict(i,zip(enum_name.(i)))
end

const nodetypemap = name2enum(Val(NodeType))

typemap(s::AbstractString) = get(nodetypemap,lowercase(s),Node)

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

mutable struct Genealogy{D <: Enum, T <: Time}
    t0::T
    time::T
    safe::Bool
    nodes::Vector{GenealNode{D,T}}
    Genealogy{D}(t0::T) where {D <: Enum, T <: Time} = begin
        new{D,T}(t0,t0,true,GenealNode{D,T}[])
    end
end

AbstractTrees.ParentLinks(::Type{<:GenealNode}) = AbstractTrees.StoredParents()
AbstractTrees.SiblingLinks(::Type{<:GenealNode}) = AbstractTrees.ImplicitSiblings()
AbstractTrees.ChildIndexing(::Type{<:GenealNode}) = AbstractTrees.IndexedChildren()
AbstractTrees.NodeType(::Type{<:GenealNode}) = HasNodeType()
AbstractTrees.childindices(tree, idx) = tree.nodes[idx].children
AbstractTrees.parentindex(tree, idx) = tree.nodes[idx].parent
AbstractTrees.nodevalue(tree, idx) = tree.nodes[idx]
Base.IteratorEltype(::Type{<:TreeIterator{<:GenealNode}}) = Base.HasEltype()
Base.eltype(::Type{<:TreeIterator{<:GenealNode{D,T}}}) where {D,T} = GenealNode{D,T}

"""
    cap_tips!(G)

Converts tip-nodes to sample-nodes.  The genealogy remains correct.
"""
cap_tips!(G::Genealogy) = begin
    for n ∈ G.nodes
        if isempty(n.children)
            n.type = Sample
        end
    end
    nothing
end

"""
    drop_zlb!(G)

Removes zero-length branches from Genealogy G as needed.
The genealogy is now incorrect: and needs to be repaired (see repair!).
"""
drop_zlb!(G::Genealogy) = begin
    for n ∈ G.nodes
        if !isnothing(n.parent)
            p = G.nodes[n.parent]
            if n.slate == p.slate && n.deme === p.deme
                setdiff!(p.children,n.name)
                append!(p.children,n.children)
                empty!(n.children)
                n.parent = nothing
                @assert p.type==Node
                p.type = n.type
            end
        end
    end
    G.safe = false
    nothing
end

"""
    weed!(G)

Removes all dead roots.
The genealogy becomes incorrect and needs to be repaired (see repair!).
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

"""
    scan_branch!(input, G, p, mapper)

Parse the branch-string in `input`, appending the corresponding
node to Genealogy G.

Arguments:
- `input`: the string containing the branch information
- `G`: the genealogy to be modified
- `p`: the name of the parent node
- `mapper`: a Dict mapping the deme-specification string to the correct deme.
"""
scan_branch!(
    s::String,
    G::Genealogy{D,T},
    p::Name,
    mapper::Dict{String,D}
) where {D <: Enum, T <: Time} = begin
    m = match(r":(.+)$",s)
    if isnothing(m)
        bl = zero(T)
    else
        bl = parse(T,m.captures[1])
    end
    m = match(r"^.*\[&&PhyloPOMP.*deme=(\w+).*\].*$"i,s)
    if isnothing(m)
        deme = missing
    else
        deme = get(mapper,lowercase(m.captures[1]),missing)
    end
    m = match(r"^.*\[&&PhyloPOMP.*type=(\w+).*\].*$"i,s)
    if isnothing(m)
        type = Node
    else
        type = typemap(m.captures[1])
    end
    q = length(G.nodes)+1
    slate = G.nodes[p].slate + bl
    n = GenealNode{D}(q,slate,deme,type)
    n.parent = p
    push!(G.nodes[p].children,q)
    push!(G.nodes,n)
    if (n.slate > G.time)
        G.time = n.slate
    end
    q
end

"""
    parse_newick(input, Val(DemeType), t0, time)

Parse the Newick-format string `input`.

Arguments:
- `t0` is the assumed root-time.
- `time` is the (optional) final-time.
- `DemeType` (of type Enum) enumerates the demes.
"""
parse_newick(
    input::AbstractString,
    ::Val{D},
    t0::T = zero(T),
    time::Union{Missing,Time} = missing,
) where {D <: Enum, T <: Time} = begin
    nsample = count(')',input)    # number of tips
    dememapper = name2enum(Val(D))
    G = Genealogy{D}(t0)
    sizehint!(G.nodes,2*nsample)
    p = 0
    tf = t0
    open = false
    stack = 0
    sqstack = 0
    f = firstindex(input)
    e = lastindex(input)
    b = e
    if input[b] != ';'
        error("invalid Newick format: no final semicolon.")
    end
    while b >= f
        if input[b]==';'
            if stack != 0
                error("invalid Newick: unbalanced parentheses.")
            end
            p = length(G.nodes)+1
            n = GenealNode{D}(p,t0)
            push!(G.nodes,n)
            e = b = b-1
            open = true
        elseif input[b] == ')'
            if (open)
                q = scan_branch!(input[(b+1):e],G,p,dememapper)
                p = q
            else
                error("invalid Newick: missing comma or semicolon.")
            end
            stack += 1
            e = b = b-1
            open = true
        elseif input[b] == '('
            if open
                scan_branch!(input[b+1:e],G,p,dememapper)
            end
            p = G.nodes[p].parent
            e = b = b-1
            stack -= 1
            open = false
        elseif input[b] == ','
            if stack <= 0
                error("invalid Newick string: misplaced comma or unbalanced parentheses.")
            end
            if open
                scan_branch!(input[b+1:e],G,p,dememapper)
            end
            e = b = b-1
            open = true
        elseif input[b] == ']'
            sqstack += 1
            while b >= f && sqstack > 0
                b = b-1
                if input[b] == ']'
                    sqstack += 1
                elseif input[b] == '['
                    sqstack -= 1
                end
            end
            if sqstack != 0
                error("invalid Newick format: unbalanced square brackets.")
            else
                b = b-1
            end
        elseif input[b] == '['
            error("invalid Newick: unbalanced square brackets.")
        else
            b = b-1
        end
    end
    if stack != 0
        error("invalid Newick format: unbalanced parentheses.")
    end
    if open
        scan_branch!(input[b+1:e],G,p,dememapper)
    end
    if !ismissing(time)
        if G.time > T(time)
            error("final time from data ($tf) exceeds specified final time ($time)")
        else
            G.time = T(time)
        end
    end
    cap_tips!(G)
    drop_zlb!(G)
    weed!(G)
    repair!(G)
    G
end

parse_newick(input::AbstractVector{V},args...) where {V<:AbstractString} =
    parse_newick(join(input),args...)

@enum MERSDemes begin
    Camel
    Human
end

x1 = readlines("data/MERS_274_sCoal_phylopomp.nwk")

g1 = parse_newick(x1,Val(MERSDemes),0.0)

@enum Strains begin
    Strain1
    Strain2
    Strain3
end

x2 = readlines("data/B.1.617.all.nwk")
g2 = parse_newick(x2,Val(Strains),0.0)

parse_newick("():0.1;",Val(Strains),0.0)

try
    parse_newick("(:0.1;",Val(Strains),0.0)
catch e
    if hasproperty(e,:msg)
        @info "error: $(e.msg)"
    end
end

try
    parse_newick("[bob=3 tom=[]:0.1;",Val(Strains),0.0)
catch e
    if hasproperty(e,:msg)
        @info "error: $(e.msg)"
    end
end
