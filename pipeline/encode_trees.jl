# M5: encode trees.tsv (from generate_trees.jl) into the CNN-CBLV input of Voznica et al., as phylodeep's
# `encode_into_most_recent` does (used by encode_generated_trees.py). Python and phylodeep are a reference only.
#
#   julia -t 16 --project=pipeline pipeline/encode_trees.jl <run_dir>
#
# Writes <run_dir>/cblv.h5 with datasets `cblv` (402 or 1002 rows: the encoded tree with p_sample), `ss` (99 rows:
# the summary statistics of sumstats.jl, then p_sample), `rescale` (phylodeep's rescale factor, the same for both
# encodings), `targets` (the target columns of trees.tsv: those between `tree` and `p_sample`), `target_names`, and
# `tree` (row of trees.tsv). One tree per column. Columns are found by their header names.

using PhyloPOMP
using PhyloPOMP: Root
using HDF5
include(joinpath(@__DIR__, "sumstats.jl"))

"""
    phylodeep_cblv(g, p_sample) -> (v, rescale)

The CBLV vector of a single-root genealogy `g` as phylodeep's `encode_into_most_recent` builds it.
The root edge is dropped, branch lengths are divided by `rescale`, their mean, and the encoding
`[tip₁, node₁, tip₂, node₂, …, tipₙ]` is padded with zeros to 399 entries (fewer than 200 tips) or
999, then reordered with two copies of `p_sample` as phylodeep does: internal-node half first, tip
half second. `length(v)` is 402 or 1002.
"""
phylodeep_cblv(g::Genealogy, p::Real) = begin
    x, y = cblv(g)
    n = length(x)
    length(collect(roots(g))) == 1 || throw(ArgumentError("one root expected"))
    pop!(y)                              # the zero cblv appends for the root
    ## cblv measures depths from the root; phylodeep from the first branch point.
    l0 = minimum(y)
    x[1] -= l0
    y .-= l0
    rescale = branch_mean(g, l0)
    e = Vector{Float64}(undef, 2n - 1)
    e[1:2:end] .= x ./ rescale
    e[2:2:end] .= y ./ rescale
    maxlen = n < 200 ? 399 : 999
    full = zeros(maxlen + 3)
    full[1:length(e)] .= e
    full[maxlen + 2] = p                 # 0-based positions maxlen+1 and maxlen+2
    full[maxlen + 3] = p
    ints = [maxlen; 1:2:(maxlen - 4); maxlen + 2; maxlen - 2]
    tips = [0:2:(maxlen - 3); maxlen + 1; maxlen - 1]
    full[[ints; tips] .+ 1], rescale
end

## Mean branch length as phylodeep's `rescale_tree`: over every node of the tree below the root edge,
## each node's branch to its parent. The first branch point keeps its own branch, the root edge `l0`,
## because ete3's `detach` leaves `dist` unchanged; without it the encoding differs by up to 0.18.
branch_mean(g::Genealogy, l0) = begin
    s = l0; k = 1
    for i in eachindex(g)
        nd = g[i]
        nd.type == Root && continue
        par = nd.parent
        (isnothing(par) || g[par].type == Root) && continue   # the root edge, counted above
        s += nd.slate - g[par].slate
        k += 1
    end
    s / k
end

main(args) = begin
    length(args) == 1 || error("usage: encode_trees.jl <run_dir>")
    dir = args[1]
    lines = eachline(joinpath(dir, "trees.tsv"))
    header = split(first(lines), '\t')
    col(name) = findfirst(==(name), header)::Int
    ip, inw, itree = col("p_sample"), col("newick"), col("tree")
    tcols = (itree + 1):(ip - 1)
    rows = [split(l, '\t') for l in lines]
    n = length(rows)
    enc = Vector{Vector{Float64}}(undef, n)
    ss = Vector{Vector{Float64}}(undef, n)
    resc = Vector{Float64}(undef, n)
    Threads.@threads :dynamic for i in 1:n
        r = rows[i]
        g, p = parse_newick(String(r[inw])), parse(Float64, r[ip])
        enc[i], resc[i] = phylodeep_cblv(g, p)
        ss[i], s2 = phylodeep_sumstats(g, p)
        isapprox(s2, resc[i]; rtol = 1e-12) || error("tree $(r[1]): rescale factors differ ($s2, $(resc[i]))")
    end
    widths = unique(length.(enc))
    length(widths) == 1 || error("trees of both sizes in one run: widths $widths")
    h5open(joinpath(dir, "cblv.h5"), "w") do f
        f["cblv"] = Float32.(reduce(hcat, enc))
        f["ss"] = reduce(hcat, ss)                  # Float64; standardized at training time
        f["rescale"] = resc
        f["targets"] = [parse(Float64, r[c]) for c in tcols, r in rows]
        f["target_names"] = String.(header[tcols])
        f["tree"] = [parse(Int, r[itree]) for r in rows]
    end
    println("$n trees, $(only(widths)) columns -> $(joinpath(dir, "cblv.h5"))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
