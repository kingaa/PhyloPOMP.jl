"""
    genealogyplot(g::Genealogy; ladderize = true, points = true, palette, title = "")

Draw `g` as a tree, time running left to right.
Returns a Makie `Figure`.
Samples are dots, colored by deme when the genealogy records one.

- `ladderize`: order the branches below each node by clade size.
- `points`: draw a dot at each sample.
- `palette`: colors for the demes, indexed by deme number.
- `title`: title of the axis.
- `figure`, `axis`: `NamedTuple`s of keywords passed to `Figure` and `Axis`.

Requires Makie or one of its backends, such as `CairoMakie`.
Load it with `import CairoMakie` or `using CairoMakie`.
"""
function genealogyplot end

"""
    treelayout(g; ladderize = true)

Where to draw each part of `g` when time runs along the horizontal axis.
Returns a `NamedTuple`:

- `y`: vertical position of each node. Tips are `1, 2, ...`, in drawing order.
  Other nodes sit halfway between their lowest and highest child.
- `branches`: one `(x0, x1, y)` for each node with a parent, the horizontal line
  from the parent's time to the node's time.
- `connectors`: one `(x, y0, y1)` for each node with children, the vertical line
  joining its children.
"""
treelayout(g::Genealogy; ladderize::Bool = true) = begin
    n = length(g)
    kids = [Int.(g[i].children) for i ∈ 1:n]
    ntip = ones(Int, n)
    ## children come after their parents, so go backwards
    for i ∈ n:-1:1
        isempty(kids[i]) || (ntip[i] = sum(ntip[c] for c ∈ kids[i]))
    end
    ladderize && foreach(k -> sort!(k; by = c -> ntip[c]), kids)
    y = zeros(n)
    stack = reverse(Int.(roots(g)))
    tip = 0
    while !isempty(stack)
        i = pop!(stack)
        if isempty(kids[i])
            y[i] = (tip += 1)
        else
            append!(stack, reverse(kids[i]))
        end
    end
    for i ∈ n:-1:1
        isempty(kids[i]) || (y[i] = (minimum(y[kids[i]]) + maximum(y[kids[i]])) / 2)
    end
    branches = NTuple{3,Float64}[]
    connectors = NTuple{3,Float64}[]
    for i ∈ 1:n
        nd = g[i]
        isnothing(nd.parent) ||
            push!(branches, (g[Int(nd.parent)].slate, nd.slate, y[i]))
        isempty(kids[i]) ||
            push!(connectors, (nd.slate, minimum(y[kids[i]]), maximum(y[kids[i]])))
    end
    (y = y, branches = branches, connectors = connectors)
end
