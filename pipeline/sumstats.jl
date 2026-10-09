# M5: the 98 summary statistics of Voznica et al. (2022) plus p_sample, as phylodeep's
# `encode_into_summary_statistics` computes them (phylodeep/sumstats.py, Python reference only).
# Each function below transcribes the phylodeep function of the same name, including its quirks
# (noted where they matter); pipeline/test/runtests.jl checks the result against phylodeep.

using PhyloPOMP
using PhyloPOMP: Root
using Statistics

## The tree as ete3 holds it after `prune(..., preserve_branch_length = true)`: nodes with one child
## removed (their branch added to the child's), so the root is the first branch point and keeps the
## root edge as its own `dist`. Node 1 is the root; `children` keeps the genealogy's child order.
struct ETree
    parent::Vector{Int}
    children::Vector{Vector{Int}}
    dist::Vector{Float64}
end

etree(g::Genealogy) = begin
    r = only(i for i in eachindex(g) if g[i].type == Root)
    parent = Int[]; children = Vector{Int}[]; dist = Float64[]
    ## follow single-child chains from genealogy node `n`, reached by a branch of length `len`
    add(n, p, len) = begin
        while length(g[n].children) == 1
            c = only(g[n].children)
            len += g[c].slate - g[n].slate
            n = c
        end
        push!(parent, p); push!(children, Int[]); push!(dist, len)
        k = length(parent)
        p > 0 && push!(children[p], k)
        for c in g[n].children
            add(c, k, g[c].slate - g[n].slate)
        end
    end
    add(r, 0, 0.0)
    ETree(parent, children, dist)
end

isleaf(t::ETree, i) = isempty(t.children[i])
levelorder(t::ETree) = begin
    out = [1]; k = 1
    while k ≤ length(out)
        append!(out, t.children[out[k]]); k += 1
    end
    out
end
nleaves(t::ETree) = begin
    n = zeros(Int, length(t.parent))
    for i in reverse(levelorder(t))
        n[i] = isleaf(t, i) ? 1 : sum(n[c] for c in t.children[i])
    end
    n
end

## numpy's var (ddof = 0) and median
pvar(x) = var(x; corrected = false)

"""
    phylodeep_sumstats(g, p_sample) -> (v, rescale)

The 99 inputs of FFNN-SS (98 statistics, then `p_sample`) for genealogy `g`, and the rescale factor,
as phylodeep's `encode_into_summary_statistics`.
"""
phylodeep_sumstats(g::Genealogy, p::Real) = begin
    t = etree(g)
    N = length(t.parent)
    rescale = mean(t.dist)                               # rescale_tree: every node, root edge included
    d = t.dist ./ rescale
    lo = levelorder(t)
    leaves = [i for i in lo if isleaf(t, i)]
    internal = [i for i in lo if !isleaf(t, i)]
    nl = nleaves(t)
    nleaf = length(leaves)

    ## add_depth_and_get_max: children of the root have depth 1; the maximum skips depth-1 nodes
    depth = zeros(Int, N); maxdepth = 0
    for i in lo[2:end]
        depth[i] = t.parent[i] == 1 ? 1 : depth[t.parent[i]] + 1
        t.parent[i] != 1 && (maxdepth = max(maxdepth, depth[i]))
    end
    dtr = zeros(N)                                       # add_dist_to_root: the root is at 0
    for i in lo[2:end]
        dtr[i] = dtr[t.parent[i]] + d[i]
    end
    ## add_ladder
    ladder = fill(-1, N)
    for i in lo[2:end]
        isleaf(t, i) && continue
        c1, c2 = t.children[i][1], t.children[i][2]
        if t.parent[i] == 1
            ladder[i] = (isleaf(t, c1) || isleaf(t, c2)) ? 0 : -1
        elseif isleaf(t, c1) && isleaf(t, c2)
            ladder[i] = 0
        elseif isleaf(t, c1) || isleaf(t, c2)
            ladder[i] = ladder[t.parent[i]] + 1
        else
            ladder[i] = 0
        end
    end
    ## add_height: leaves 0, internal nodes 1 + highest child
    height = zeros(Int, N)
    for i in reverse(lo)
        isleaf(t, i) || (height[i] = 1 + maximum(height[c] for c in t.children[i]))
    end

    s = Float64[]
    ## tree_height
    th = [dtr[i] for i in leaves]
    append!(s, (maximum(th), minimum(th)))
    ## branches
    dall, dext = d[lo], d[leaves]
    append!(s, (mean(dall), median(dall), pvar(dall), mean(dext), median(dext), pvar(dext)))
    ## piecewise_branches: internal nodes by distance to the root, in thirds of the tree height
    allmax, emean, emed, evar = s[1], s[6], s[7], s[8]
    parts = (Float64[], Float64[], Float64[])
    for i in internal
        if dtr[i] < allmax / 3
            push!(parts[1], d[i])
        elseif dtr[i] < 2allmax / 3
            push!(parts[2], d[i])
        elseif dtr[i] > 2allmax / 3                     # a node exactly at 2/3 is in no part
            push!(parts[3], d[i])
        end
    end
    for x in parts
        append!(s, isempty(x) ? zeros(6) :
                (mean(x), median(x), pvar(x), mean(x) / emean, median(x) / emed, pvar(x) / evar))
    end
    ## colless (at a polytomy phylodeep replaces the running sum by the polytomy's mean; kept)
    col = 0.0
    for i in internal
        ch = t.children[i]
        if length(ch) == 2
            col += abs(nl[ch[1]] - nl[ch[2]])
        else
            col = mean(abs(nl[ch[j]] - nl[ch[k]]) for j in eachindex(ch) for k in (j + 1):length(ch))
        end
    end
    push!(s, col)
    ## sackin
    push!(s, sum(depth[i] for i in leaves))
    ## wd_ratio_delta_w (phylodeep compares width[0] with width[-1], the last entry; kept)
    width = zeros(maxdepth + 1)
    for i in lo[2:end]
        width[depth[i] + 1] += 1
    end
    δw = 0.0
    for j in 0:(length(width) - 2)
        prev = j == 0 ? width[end] : width[j]
        δw = max(δw, abs(width[j + 1] - prev))
    end
    append!(s, (maximum(width) / maxdepth, δw))
    ## max_ladder_il_nodes (internal nodes, root included with ladder -1)
    append!(s, (max(0, maximum(ladder[i] for i in internal)) / nleaf,
                count(i -> ladder[i] > 0, internal) / (nleaf - 1)))
    ## staircaseness
    ratios = [min(nl[t.children[i][1]], nl[t.children[i][2]]) / max(nl[t.children[i][1]], nl[t.children[i][2]])
              for i in internal]
    append!(s, (count(i -> nl[t.children[i][1]] != nl[t.children[i][2]], internal) / (nleaf - 1), mean(ratios)))
    ## ltt_plot: sampling (-1) and branching (+1) events sorted by time; at equal times phylodeep's sort on
    ## the float bits as Int64 puts -1 first
    ev = Tuple{Float64,Float64}[]
    for i in lo
        if isleaf(t, i)
            push!(ev, (dtr[i], -1.0))
        else
            foreach(_ -> push!(ev, (dtr[i], 1.0)), 2:length(t.children[i]))
        end
    end
    sort!(ev)
    tm = first.(ev); kind = last.(ev)
    lin = cumsum(kind) .+ 1                              # lineages after each event, starting from 1
    ## ltt_plot_comput
    jmax = argmax(lin)                                   # first index of the maximum
    slope(x, y) = (xm = mean(x); sum((x .- xm) .* (y .- mean(y))) / sum((x .- xm) .^ 2))   # linregress slope
    s1 = slope(lin[1:jmax-1], tm[1:jmax-1])
    s2 = slope(lin[jmax:end], tm[jmax:end])
    tmax = tm[end]
    samp = tm[kind .== -1]
    br = [tm[k] for k in eachindex(tm) if kind[k] == 1]
    br1 = filter(x -> x < tmax / 3, br)
    br2 = filter(x -> tmax / 3 ≤ x < 2tmax / 3, br)
    br3 = filter(x -> x ≥ 2tmax / 3, br)
    mdiff(x) = length(x) > 1 ? mean(diff(x)) : 0.0
    append!(s, (lin[jmax], tm[jmax], s1, s2, s1 / s2, mean(diff(samp)), mdiff(br1), mdiff(br2), mdiff(br3)))
    ## coordinates_comp: 20 bins of events, mean time then mean lineages
    n = length(ev)
    b = floor.(Int, range(0, n, length = 21))
    append!(s, [mean(tm[b[j]+1:b[j+1]]) for j in 1:20])
    append!(s, [mean(lin[b[j]+1:b[j+1]]) for j in 1:20])
    ## number of tips
    push!(s, nleaf)
    ## compute_chain_stats (order 4): from each node of height > 3, follow the shortest child branch
    chains = Float64[]
    for i in lo
        height[i] > 3 || continue
        c = Float64[]; k = i
        while length(c) < 4
            ch = t.children[k]
            j = argmin(d[ch])                            # first child among equals, as list.index(min)
            push!(c, d[ch[j]])
            k = ch[j]
            isleaf(t, k) && break
        end
        length(c) == 4 && push!(chains, sum(c))
    end
    if length(chains) > 1
        append!(s, (length(chains), mean(chains)))
        append!(s, quantile(chains, 0:0.1:1))            # numpy percentile, linear interpolation
        push!(s, pvar(chains))
    else
        append!(s, zeros(14))
    end
    push!(s, p)
    length(s) == 99 || error("expected 99 values, got $(length(s))")
    s, rescale
end
