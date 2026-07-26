"""
    parse_newick(input; demes, t0, time)

Parse the Newick-format string (or vector of strings) `input`.
It is possible to include metadata tags of format `[&&PhyloPOMP ...]` in the label-strings.
All tips are assumed to be samples.
Zero-length branches are dropped.

Optional arguments:
- `demes` is a `Module` enumerating the demes (see [`@demes`](@ref)).
  By default, `demes = Unstructured`.
- `t0` is the assumed root-time. By default, `t0 = 0`.
- `time` is the (optional) final-time.
"""
parse_newick(
    input::AbstractVector{V};
    kwargs...,
) where {V<:AbstractString} =
    parse_newick(join(reverse(input)); kwargs...)

parse_newick(
    input::AbstractString;
    demes::Module = Unstructured,
    t0::Real = zero(Time),
    time::Union{Missing,Real} = missing,
) = begin
    @isademeset demes
    D = demes.DemeSet
    dememapper = name2enum(D)
    tf = t0 = Time(t0)
    bl = zero(Time)
    G = Genealogy{demes}(t0)
    nnodes = count(')',input)+count(',',input)+2*count(';',input)
    sizehint!(G.nodes,nnodes)
    p::Union{Nothing,Name} = nothing
    open::Bool = false
    stack::Int = 0
    sqstack::Int = 0
    f = firstindex(input)
    b = e = lastindex(input)
    @assert input[b] == ';' "invalid Newick: no final semicolon."
    while b >= f
        @assert input[b] != '[' "invalid Newick: unbalanced square brackets."
        if input[b]==';'
            @assert stack == 0 "invalid Newick: unbalanced parentheses."
            if open
                scan_branch!(G, input[(b+1):e], p, dememapper, bl)
            end
            p = length(G.nodes)+1
            n = GenealNode{D}(p,t0)
            push!(G.nodes,n)
            e = b = b-1
            open = true
            bl = zero(Time)
        elseif input[b] == ')'
            @assert open "invalid Newick: missing comma or semicolon."
            scan_branch!(G, input[(b+1):e], p, dememapper, bl)
            p = length(G.nodes)
            stack += 1
            e = b = b-1
            open = true
            bl = zero(Time)
        elseif input[b] == '('
            if open
                scan_branch!(G, input[(b+1):e], p, dememapper, bl)
            end
            p = G[p].parent
            e = b = b-1
            stack -= 1
            open = false
        elseif input[b] == ','
            @assert stack > 0 "invalid Newick: misplaced comma or unbalanced parentheses."
            if open
                scan_branch!(G, input[(b+1):e], p, dememapper, bl)
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
            @assert sqstack == 0 "invalid Newick: unbalanced square brackets."
            b = b-1
        elseif input[b] == ':'
            @assert open "invalid Newick: misplaced colon."
            bl = scan_length(input[(b+1):e])
            e = b = b-1
        else
            b = b-1
        end
    end
    @assert stack == 0 "invalid Newick: unbalanced parentheses."
    if open
        scan_branch!(G, input[(b+1):e], p, dememapper, bl)
    end
    set_time!(G,time)
    cap_tips!(G)     # all tips become samples
    clip_zlb!(G)     # samples with zero-length branches become inline
    repair!(G)
    G
end

const nodetypemap = Dict(
    "node"=>Node,
    "branch"=>Node,
    "migration"=>Node,
    "root"=>Node,
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

"""
    clip_zlb!(G)

Isolates zero-length branches from Genealogy `G` as needed.  The
genealogy is now incorrect: and needs to be repaired (see
[`repair!`](@ref)).
"""
clip_zlb!(G::Genealogy) = begin
    for n ∈ G.nodes
        if !isnothing(n.parent)
            p = G[n.parent]
            if n.slate == p.slate
                if (isnothing(p.parent) && length(n.children)==1) ||
                    (!isnothing(p.parent) && n.deme===p.deme)
                    for c ∈ n.children
                        G[c].parent = p.name
                    end
                    setdiff!(p.children,n.name)
                    append!(p.children,n.children)
                    empty!(n.children)
                    n.parent = nothing
                    @assert p.type == Node "dropping zero-length branch collapses multiple samples."
                    p.type = n.type
                end
            end
        end
    end
    nothing
end

"""
    insert_zlb!(G)

Adds zero-length branches where needed.  The genealogy is now
incorrect: and needs to be repaired (see [`repair!`](@ref)).
"""
insert_zlb!(G::Genealogy{D}) where D = begin
    for n ∈ G.nodes
        if n.type==Sample && !isempty(n.children)
            q = length(G)+1
            node = GenealNode{D.DemeSet}(q,n.slate,n.deme,Sample,n.name)
            push!(G.nodes,node)
            push!(n.children,q)
            n.type = Node
        end
    end
    nothing
end


"""
    scan_branch!(G, input, p, mapper, bl)

Parse the branch-string in `input`, appending the corresponding node
to Genealogy `G`.

Arguments:  
- `G`: the genealogy to be modified
- `input`: the string containing the branch information
- `p`: the name of the parent node
- `mapper`: a function that maps strings to demes.
- `bl`: the branch length
"""
scan_branch!(
    G::Genealogy{D},
    input::AbstractString,
    p::Name,
    dememapper::Function,
    bl::Time,
) where D = begin
    m = match(r"^.*\[&&PhyloPOMP.*deme=(\w+).*\].*$"i, input)
    if isnothing(m)
        deme = missing
    else
        deme = dememapper(m.captures[1])
        @assert !ismissing(deme) "unrecognized deme '$(m.captures[1])'."
    end
    m = match(r"^.*\[&&PhyloPOMP.*type=(\w+).*\].*$"i, input)
    if isnothing(m)
        type = Node
    else
        type = get(nodetypemap,lowercase(m.captures[1]),missing)
        @assert !ismissing(type) "unrecognized type '$(m.captures[1])'."
    end
    q = length(G.nodes)+1
    slate = G[p].slate + bl
    n = GenealNode{D.DemeSet}(q,slate,deme,type)
    n.parent = p
    push!(G[p].children,q)
    push!(G.nodes,n)
    if (n.slate > G.time)
        G.time = n.slate
    end
    nothing
end

scan_length(input::AbstractString) = begin
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
    @assert bl >= zero(Time) "negative branch length detected."
    bl
end
