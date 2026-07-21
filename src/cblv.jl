"""
    ladderize(g)

Ladderize the genealogy `g`.  The children of each internal node are
sorted in order of decreasing subtree height.  Return the ladderized
genealogy and the indices of the root nodes ordered in the same way.
"""
ladderize(g::Genealogy) = begin
    g = deepcopy(g)
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
            push!(y,node.slate)
            memo[n] = true
            visit!(x,y,memo,g,c)
        end
    end
    nothing
end

"""
    cblv(g)

Returns a compact, bijective, ladderized vector (CBLV) representation
of the genealogy `g`.
"""
cblv(g::Genealogy) = begin
    x = Time[]
    y = Time[]
    sizehint!(x,nsample(g))
    sizehint!(y,nsample(g))
    memo = fill(false,length(g))
    g, r = ladderize(g)
    for n ∈ r
        visit!(x,y,memo,g,n)
        push!(y,zero(Time))
        memo[n] = true
    end
    @assert all(memo)
    x, y
end
