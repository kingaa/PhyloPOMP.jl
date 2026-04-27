"""
    newick(g; extended = true, sigdigits = 6)

Generates a `String` containing a Newick-format description of the [`Genealogy`](@ref) `g`.
If `extended = true` (the default), then the **PhyloPOMP**-flavored Newick extension is used (i.e., metadata tags are included).
"""
newick(
    g::Genealogy;
    extended::Bool = true,
    sigdigits::Integer = 6
) = begin
    if extended
        extnewick(g,sigdigits)
    else
        plainnewick(g,sigdigits)
    end
end

plainnewick(
    g::Genealogy,
    sigdigits::Integer = 6,
) = begin
    g = deepcopy(g)
    insert_zlb!(g)
    nstr = similar(Array{String},length(g))
    for i ∈ reverse(eachindex(g))
        n = g[i].name
        bstr = if isnothing(g[i].parent)
            "$n:0;"
        else
            bl = round(g[i].slate-g[g[i].parent].slate,sigdigits=sigdigits)
            "$n:$bl"
        end
        nstr[i] = if isempty(g[i].children)
            bstr
        else
            c = join(map(j -> nstr[j], g[i].children),',')
            "("*c*")"*bstr
        end
    end
    nstr[roots(g)]
end

extnewick(
    g::Genealogy{D},
    sigdigits::Integer = 6,
) where D = begin
    typenodemap = Dict(Node=>"node",Sample=>"sample",Root=>"root")
    nstr = similar(Array{String},length(g))
    dememap = enum2name(Val(D))
    for i ∈ reverse(eachindex(g))
        t = typenodemap[g[i].type]
        d = get(dememap,g[i].deme,missing)
        dstr = if ismissing(d)
            ""
        else
            " deme=$d"
        end
        n = g[i].name
        bstr = if isnothing(g[i].parent)
            "[&&PhyloPOMP type=root]$n:0;"
        else
            bl = round(g[i].slate-g[g[i].parent].slate,sigdigits=sigdigits)
            "[&&PhyloPOMP type=$t$dstr]$n:$bl"
        end
        nstr[i] = if isempty(g[i].children)
            bstr
        else
            c = join(map(j -> nstr[j], g[i].children),',')
            "("*c*")"*bstr
        end
    end
    nstr[roots(g)]
end
