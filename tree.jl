using AbstractTrees
import NewickTree: readnw, getheights, getleaves, prewalk, Node, NewickData, distance, nv

@enum NodeType begin
    Root
    Sample
    Branch
end

x = readlines("data/MERS_274_sCoal_phylopomp.nwk")

x = [
    "(((((b_1_1:0.1,s_0_2:0.1)g_0_1:0.2,b_0_3:0.3)g_0_2:0.4,b_1_4:0.5)g_0_3:0.6)b_0_8:0.2,b_0_5:0.7)m_0_4:0.8;",
    "(b_0_6:0.1,b_0_7:0.1)m_0_5:0.5;"
]

repairRoots(x::Node{I,T}, i::J) where {I<:Integer,J<:Integer,D,S,T<:NewickData{D,S}} = begin
    bl = distance(x)
    if bl != 0
        m = Node(I(i + 1), NewickData(d=zero(D)))
        push!(m, x)
    else
        x
    end
end

repairRoots(x::Node{I,T}) where {I,T} = begin
    n = nv(x)
    y = similar(x.children)
    for k in eachindex(y)
        y[k] = repairRoots(pop!(x), n += 1)
    end
    for k in eachindex(y)
        push!(x, y[k])
    end
    x
end

parseNewick(
    input::Vector{<:AbstractString},
    time::Real;
    t0::Real=0.0,
) = begin
    x = "(" * join(map(n -> match(r"^(.+);\w*$", n).captures[1], input), ",") * ");" |> readnw |> repairRoots
    ids = [n.id for n in prewalk(x)]
    n = sort(Dict(zip(ids, eachindex(ids) .- 1)), byvalue=true)
    h = getheights(x)
    ids = sort(ids, lt=(i1, i2) -> ((h[i1] < h[i2]) || (h[i1] == h[i2] && n[i1] < n[i2])))
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
    seq = [(n=n[id],
        id=id,
        lineage=lineages[id],
        parent=begin
            p = getindex(n, parent[id])
            p == 0 ? n[id] : p
        end,
        children=map(i -> n[i], children[id]),
        time=t0 + h[id],
        duration=dur[id],
        type=type[id],
        deme=deme[id],
        saturation=length(children[id]),
    ) for id in Base.rest(ids, 2)]
end

parseNewick(x, 5.0)
