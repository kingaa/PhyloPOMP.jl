using AbstractTrees
import OrderedCollections: OrderedDict
import NewickTree: readnw, getheights, getleaves, prewalk, Node, NewickData, distance, nv
export parse_newick, NodeType

@enum NodeType begin
    Root
    Sample
    Branch
end

repair_roots(x::Node{I, T}, i::J) where {I <: Integer, J <: Integer, D, S, T <: NewickData{D, S}} = begin
    bl = distance(x)
    if bl != 0
        m = Node(I(i + 1), NewickData(d = zero(D)))
        push!(m, x)
    else
        x
    end
end

repair_roots(x::Node{I, T}) where {I, T} = begin
    n = nv(x)
    y = similar(x.children)
    for k in eachindex(y)
        y[k] = repair_roots(pop!(x), n += 1)
    end
    for k in eachindex(y)
        push!(x, y[k])
    end
    x
end

parse_newick(
    input::Vector{<:AbstractString},
    time::Real;
    t0::Real = 0.0,
) = begin
    x = "(" * join(map(n -> match(r"^(.+);\w*$", n).captures[1], input), ",") * ");" |> readnw |> repair_roots
    ids = [n.id for n in prewalk(x)]
    n = OrderedDict(zip(ids, eachindex(ids) .- 1))
    sort!(n, byvalue = true)
    h = getheights(x)
    ids = sort(ids, lt = (i1, i2) -> ((h[i1] < h[i2]) || (h[i1] == h[i2] && n[i1] < n[i2])))
    children = Dict([
        n.id => if isdefined(n, :children)
            map(n -> n.id, n.children)
        else
            UInt16[]
        end
        for n in prewalk(x)]
                    )
    parent = Dict([
        n.id => if isdefined(n, :parent)
            n.parent.id
        else
            missing
        end
        for n in prewalk(x)]
                  )
    labels = Dict(
        [n.id => match(r"([bsgmo])_(\d+)_(\d+)", n.data.name)
         for n in prewalk(x)]
    )
    typemap = Dict("b" => Sample, "s" => Sample, "g" => Branch, "m" => Root)
    type = Dict([
        if !isnothing(v)
            k => typemap[v.captures[1]]
        else
            k => 0
        end for (k, v) in labels
            ]
                )
    deme = Dict([
        if !isnothing(v)
            k => parse(Int, v.captures[2])
        else
            k => missing
        end for (k, v) in labels
            ]
                )
    names = Dict([
        if !isnothing(v)
            k => parse(Int, v.captures[3])
        else
            k => missing
        end for (k, v) in labels
            ]
                 )
    lineages = Dict([
        n.id => minimum(map(k -> names[k.id], getleaves(n)))
        for n in prewalk(x)]
                    )
    dur = Dict(zip(ids, diff([getindex.(Ref(h), ids)..., time])))
    seq = [
        (
            n = n[id],
            id = id,
            lineage = lineages[id],
            parent = begin
                p = getindex(n, parent[id])
                p == 0 ? n[id] : p
            end,
            children = map(i -> n[i], children[id]),
            time = t0 + h[id],
            duration = dur[id],
            type = type[id],
            deme = deme[id],
            saturation = length(children[id]),
        ) for id in Base.rest(ids, 2)
            ]
end
