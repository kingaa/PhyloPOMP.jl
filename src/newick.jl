export newick

newick(g::Genealogy{D}) where D = begin
    typenodemap = Dict(Node=>"node",Sample=>"sample")
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
            "[&&PhyloPOMP type=$t$dstr]$n:0;"
        else
            bl = g[i].slate - g[g[i].parent].slate
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
