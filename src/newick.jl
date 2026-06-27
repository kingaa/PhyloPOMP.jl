"""
    newick(g; extended = true, sigdigits = 6, demenames = true)

Generates a `String` containing a Newick-format description of the
[`Genealogy`](@ref) `g`.  If `extended = true` (the default), then the
**PhyloPOMP**-flavored Newick extension is used (i.e., metadata tags
are included).  If `demenames = false`, then demes are indicated using
integers instead of names (in the extended format only).
"""
newick(
    g::Genealogy;
    extended::Bool = true,
    sigdigits::Integer = 6,
    demenames = true,
) = begin
    if extended
        extnewick(g,sigdigits,demenames)
    else
        plainnewick(g,sigdigits)
    end
end

"""
    plainnewick(g; sigdigits = 6)

Generates a `String` containing a plain Newick representation of the [`Genealogy`](@ref) `g`.  All deme information is discarded.
"""
plainnewick(
    g::Genealogy,
    sigdigits::Integer = 6,
) = begin
    g = deepcopy(g)
    insert_zlb!(g)
    nstr = Array{String}(undef,length(g))
    for i ∈ reverse(eachindex(g))
        bstr = if isnothing(g[i].parent)
            ";"
        else
            bl = round(g[i].slate-g[g[i].parent].slate,sigdigits=sigdigits)
            ":$bl"
        end
        nstr[i] = if isempty(g[i].children)
            bstr
        else
            "("*join(map(j -> nstr[j], g[i].children),',')*")"*bstr
        end
    end
    nstr[roots(g)]
end

"""
    extnewick(g; sigdigits = 6, demenames = true)

Generates a `String` containing a representation of the [`Genealogy`](@ref) `g` in the special, **PhyloPOMP**-flavored, extended Newick format.
"""
extnewick(
    g::Genealogy{D},
    sigdigits::Integer = 6,
    demenames::Bool = true,
) where D = begin
    typenodemap = Dict(Node=>"node",Sample=>"sample",Root=>"root")
    nstr = Array{String}(undef,length(g))
    dememap = enum2name(D.DemeSet)
    for i ∈ reverse(eachindex(g))
        t = typenodemap[g[i].type]
        dstr = if ismissing(g[i].deme)
            ""
        elseif demenames
            " deme=$(get(dememap,g[i].deme,missing))"
        else
            " deme=$(Int(g[i].deme))"
        end
        bstr = if isnothing(g[i].parent)
            ";"
        else
            bl = round(g[i].slate-g[g[i].parent].slate,sigdigits=sigdigits)
            "[&&PhyloPOMP type=$t$dstr]:$bl"
        end
        nstr[i] = if isempty(g[i].children)
            bstr
        else
            "("*join(map(j -> nstr[j], g[i].children),',')*")"*bstr
        end
    end
    nstr[roots(g)]
end
