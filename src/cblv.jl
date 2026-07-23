"""
    ladderize(g)

Ladderize the genealogy `g`.  The children of each internal node are
sorted in order of decreasing subtree height.  Return the ladderized
genealogy and the indices of the root nodes ordered in the same way.
"""
ladderize(g::Genealogy) = begin
    g = deepcopy(g)
    insert_zlb!(g)
    repair!(g)
    heights = Vector{Time}(undef,length(g))
    for n ∈ Iterators.reverse(eachindex(g))
        ladderize!(heights,g,n)
    end
    r = roots(g)
    p = sortperm(map(c -> heights[c], r), rev=true)
    g, r[p]
end

ladderize!(
    heights::Vector{Time},
    g::Genealogy,
    n::Integer,
) = begin
    node = g[n]
    if isempty(node.children)
        heights[n] = Time(node.slate)
    else
        h = map(c -> heights[c], node.children)
        heights[n] = maximum(h)
        p = sortperm(h, rev=true)
        node.children = node.children[p]
    end
    nothing
end

## The distance along the genealogy to the first already-visited node.
added_branch_length(
    g::Genealogy,
    memo::Vector{Bool},
    n::Integer,
) = begin
    p = g[n].parent
    if isnothing(p)
        zero(Time)
    elseif memo[p]
        Time(g[n].slate - g[p].slate)
    else
        Time(g[n].slate - g[p].slate) + added_branch_length(g,memo,p)
    end
end

"""
    visit!(x, y, memo, g, n)

Recursively perform an in-order walk over the subtree of genealogy `g`
rooted at node `n`.  Push the length of the new branch for each tip
node into `x` and for each internal node into `y`.  The Boolean vector
`memo` is used to track which nodes have been visited.
"""
visit!(
    x::Vector{Time},
    y::Vector{Time},
    memo::Vector{Bool},
    g::Genealogy,
    n::Integer,
) = begin
    node = g[n]
    @assert !memo[n]
    if isempty(node.children)
        push!(x,added_branch_length(g,memo,n))
        memo[n] = true
    else
        visit!(x,y,memo,g,node.children[1])
        for c ∈ Base.rest(node.children,2)
            push!(y,node.slate-timezero(g))
            memo[n] = true
            visit!(x,y,memo,g,c)
        end
    end
    nothing
end

"""
    cblv(g)

Returns a compact, bijective, ladderized vector (CBLV) representation
of the genealogy `g`.  See also [`parse_cblv`](@ref).
"""
cblv(g::Genealogy) = begin
    g, r = ladderize(g)
    x = Time[]
    y = Time[]
    sizehint!(x,nsample(g))
    sizehint!(y,nsample(g))
    memo = fill(false,length(g))
    for n ∈ r
        visit!(x,y,memo,g,n)
        push!(y,zero(Time))
        memo[n] = true
    end
    if !all(memo)
        @warn "dropping $(sum(.!memo)) inline nodes"
    end
    x, y
end

"""
    parse_cblv(x, y; demes, t0, time)

Parses the compact, bijective, ladderized vector (CBLV) representation
contained in `x`, `y` into a `Genealogy`.  See also [`cblv`](@ref).

Optional arguments:
- `demes` is a `Module` enumerating the demes (see [`@demes`](@ref)).
  By default, `demes = Unstructured`.
- `t0` is the assumed root-time. By default, `t0 = 0`.
- `time` is the (optional) final-time.
"""
parse_cblv(
    x::Vector{F},
    y::Vector{F};
    demes::Module = Unstructured,
    t0::Real = zero(Time),
    time::Union{Missing,Real} = missing,
) where {F <: Real} = begin
    @assert length(x)==length(y) "length mismatch"
    @assert length(x)>0 "empty CBLV representation"
    @assert maximum(x) == x[1] "improper CBLV representation"
    @assert maximum(y) <= x[1] "improper CBLV representation"
    @assert minimum(x) >= 0 "improper CBLV representation"
    @assert minimum(y) == 0 "improper CBLV representation"
    @isademeset demes
    E = demes.DemeSet
    G = Genealogy{demes}(Time(t0),Time(t0)+Time(x[1]))
    sizehint!(G.nodes,2*length(x))
    p::Name = 0
    name::Name = 0
    for k ∈ eachindex(x)
        if p == 0
            node = GenealNode{E}(
                name += 1, Time(t0),
                missing, Root, nothing
            )
            p = node.name
            push!(G.nodes,node)
        end
        tip = GenealNode{E}(
            name += 1, G[p].slate+Time(x[k]),
            missing, Sample, p
        )
        push!(G[p].children,tip.name)
        push!(G.nodes,tip)
        t = t0+Time(y[k])
        if t > t0
            i = p
            j = tip.name
            @assert G[j].slate >= t "invalid CBLV"
            while !isnothing(i) && G[i].slate > t
                j = i
                i = G[i].parent
            end
            @assert !isnothing(i)
            node = GenealNode{E}(
                name += 1, t,
                missing, Node, i
            )
            push!(node.children,j)
            setdiff!(G[i].children,j)
            push!(G[i].children,node.name)
            G[j].parent = node.name
            push!(G.nodes,node)
            p = node.name
        else
            p = 0
        end
    end
    @assert p == 0 "invalid CBLV: last value of `y` is nonzero."
    set_time!(G,time)
    cap_tips!(G)
    clip_zlb!(G)
    repair!(G)
    G
end

"""
    parse_cblv((x, y); kwargs...)

Equivalent to `parse_cblv(x, y; kwargs...)`.
"""
parse_cblv((x, y,); kwargs...,) = parse_cblv(x, y; kwargs...)
