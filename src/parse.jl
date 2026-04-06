export parse_newick

enum_name(name::E) where {E <: Enum} = lowercase(String(Symbol(name)))

name2enum(::Val{E}) where {E <: Enum} = begin
    i = instances(E)
    Dict(zip(enum_name.(i),i))
end

enum2name(::Val{E}) where {E <: Enum} = begin
    i = instances(E)
    Dict(i,zip(enum_name.(i)))
end

const nodetypemap = Dict(
    "node"=>Node,"branch"=>Node,"migration"=>Node,"root"=>Node,
    "sample"=>Sample,
)

typemap(s::AbstractString) = get(nodetypemap,lowercase(s),missing)

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

Removes zero-length branches from Genealogy `G` as needed.
The genealogy is now incorrect: and needs to be repaired (see [`repair!`](@ref)).
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
                if (p.type != Node)
                    @warn "dropping zero-length branch yields multiple samples at one node"
                end
                p.type = n.type
            end
        end
    end
    G.safe = false
    nothing
end

"""
    scan_branch!(input, G, p, mapper)

Parse the branch-string in `input`, appending the corresponding
node to Genealogy `G`.

Arguments:
- `input`: the string, or vector of strings, containing the branch information
- `G`: the genealogy to be modified
- `p`: the name of the parent node
- `mapper`: a `Dict` mapping the deme-specification string to the correct deme.
"""
scan_branch!(
    s::String,
    G::Genealogy{D,T},
    p::Name,
    mapper::Dict{String,D}
) where {D <: Enum, T <: Time} = begin
    m = match(r"^.*:([^:]+)$",s)
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
        if ismissing(deme)
            error("unrecognized deme=$(m.captures[1]).")
        end
    end
    m = match(r"^.*\[&&PhyloPOMP.*type=(\w+).*\].*$"i,s)
    if isnothing(m)
        type = Node
    else
        type = typemap(m.captures[1])
        if ismissing(type)
            error("unrecognized type=$(m.captures[1]).")
        end
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
    parse_newick(input, Val(Demes), t0, time)

Parse the Newick-format string (or vector of strings) `input`.

Arguments:
- `t0` is the assumed root-time.
- `time` is the (optional) final-time.
- `Demes` (of type `Enum`) enumerates the demes.
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
            while b > f && sqstack > 0
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
