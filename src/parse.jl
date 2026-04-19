"""
    parse_newick(input; demes, t0, time)

Parse the Newick-format string (or vector of strings) `input`.

Arguments:
- `demes` is the enumeration of the demes.
- `t0` is the assumed root-time.
- `time` is the (optional) final-time.
"""
parse_newick(
    input::AbstractVector{V};
    args...,
) where {V<:AbstractString} =
    parse_newick(join(input);args...)

parse_newick(
    input::AbstractString;
    demes::Type{D} = Unstructured,
    t0::Real = zero(Time),
    time::Union{Missing,Real} = missing,
) where {D <: Enum} = begin
    nnodes = count(')',input)+count(',',input)+2*count(';',input)
    dememapper = name2enum(demes)
    t0 = Time(t0)
    G = Genealogy{D}(t0)
    sizehint!(G.nodes,nnodes)
    p = 0
    tf = t0
    open = false
    bl = zero(Time)
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
            bl = zero(Time)
        elseif input[b] == ')'
            if (open)
                q = scan_branch!(input[(b+1):e], G, p, dememapper, bl)
                p = q
            else
                error("invalid Newick: missing comma or semicolon.")
            end
            stack += 1
            e = b = b-1
            open = true
            bl = zero(Time)
        elseif input[b] == '('
            if open
                scan_branch!(input[(b+1):e], G, p, dememapper, bl)
            end
            p = G[p].parent
            e = b = b-1
            stack -= 1
            open = false
        elseif input[b] == ','
            if stack <= 0
                error("invalid Newick string: misplaced comma or unbalanced parentheses.")
            end
            if open
                scan_branch!(input[(b+1):e], G, p, dememapper, bl)
            end
            e = b = b-1
            open = true
            bl = zero(Time)
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
        elseif input[b] == ':'
            if open
                bl = scan_length(input[(b+1):e])
                e = b = b-1
            else
                error("invalid Newick format: misplaced colon.")
            end
        else
            b = b-1
        end
    end
    if stack != 0
        error("invalid Newick format: unbalanced parentheses.")
    end
    if open
        scan_branch!(input[(b+1):e], G, p, dememapper, bl)
    end
    if !ismissing(time)
        if G.time > Time(time)
            error("final time from data ($tf) exceeds specified final time ($time)")
        else
            G.time = Time(time)
        end
    end
    cap_tips!(G)
    clip_zlb!(G)
    repair!(G)
    G
end

const nodetypemap = Dict(
    "node"=>Node,"branch"=>Node,"migration"=>Node,"root"=>Node,
    "sample"=>Sample,
)

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

##     clip_zlb!(G)
##
## Isolates zero-length branches from Genealogy `G` as needed.
## The genealogy is now incorrect: and needs to be repaired (see [`repair!`](@ref)).
clip_zlb!(G::Genealogy) = begin
    for n ∈ G.nodes
        if !isnothing(n.parent)
            p = G[n.parent]
            if n.slate == p.slate && n.deme === p.deme
                for c ∈ n.children
                    G[c].parent = p.name
                end
                setdiff!(p.children,n.name)
                append!(p.children,n.children)
                empty!(n.children)
                n.parent = nothing
                if (p.type != Node)
                    @warn "dropping zero-length branch yields multiple samples at one node."
                end
                p.type = n.type
            end
        end
    end
    nothing
end

##     insert_zlb!(G)
##
## Adds zero-length branches where needed.
## The genealogy is now incorrect: and needs to be repaired (see [`repair!`](@ref)).
insert_zlb!(G::Genealogy{D}) where D = begin
    for n ∈ G.nodes
        if n.type==Sample && !isempty(n.children)
            q = length(G)+1
            node = GenealNode{D}(q,n.slate,n.deme,Sample,n.name)
            push!(G.nodes,node)
            push!(n.children,q)
            n.type = Node
        end
    end
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
- `mapper`: a function that maps strings to demes.
"""
scan_branch!(
    input::String,
    G::Genealogy{D},
    p::Name,
    dememapper::Function,
    bl::Time,
) where {D <: Enum} = begin
    m = match(r"^.*\[&&PhyloPOMP.*deme=(\w+).*\].*$"i, input)
    if isnothing(m)
        deme = missing
    else
        deme = dememapper(m.captures[1])
        if ismissing(deme)
            error("unrecognized deme '$(m.captures[1])'.")
        end
    end
    m = match(r"^.*\[&&PhyloPOMP.*type=(\w+).*\].*$"i, input)
    if isnothing(m)
        type = Node
    else
        type = get(nodetypemap,lowercase(m.captures[1]),missing)
        if ismissing(type)
            error("unrecognized type '$(m.captures[1])'.")
        end
    end
    q = length(G.nodes)+1
    slate = G[p].slate + bl
    n = GenealNode{D}(q,slate,deme,type)
    n.parent = p
    push!(G[p].children,q)
    push!(G.nodes,n)
    if (n.slate > G.time)
        G.time = n.slate
    end
    q
end

scan_length(input::String) = begin
    m = match(
        r"^(?:\[.*?\])?([^\[\]]+?)(?:\[.*?\])?$",
        input,
    )
    if isnothing(m)
        @warn "no valid branch-length spec found in '$input', assuming zero branch-length."
        bl = zero(Time)
    else
        bl = parse(Time,m.captures[1])
    end
    if (bl < zero(Time))
        error("negative branch length detected.")
    end
    bl
end
