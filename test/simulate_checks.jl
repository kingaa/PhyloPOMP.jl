"""
Checks shared by `seir_simulate.jl`, `mers_simulate.jl`, and `sir_simulate.jl`.

A finite filter likelihood does not pin down internal-node times.
`simulate_checked` runs the simulator's own loop and checks each event as it happens.
`structure_problems` checks the finished tree: positive edge lengths, root degree 1, internal degree 2, sample degree at most `sample_max_children`.
`roundtrip_problems` checks that Newick and CBLV keep those times.
`run_invariants` runs all three over many seeds, including forests.
"""
module SimulateChecks

using Test
using Random: MersenneTwister
using PhyloPOMP
using PhyloPOMP: Root, Node, Sample, Name, apply_delta, BIRTH, SAMPLE

export simulate_checked, structure_problems, roundtrip_problems,
    total_branch_length, run_invariants

function simulate_checked(
    model, θ; x0, graft, tmax, rng,
    demeset::Module = PhyloPOMP.Unstructured, samplemap = nothing,
)
    ndeme = length(model.demes)
    bad = String[]
    ## The state just before the current event, saved by the previous call.
    n0 = Ref(0)
    before = Ref([Name[] for _ ∈ 1:ndeme])
    xprev = Ref(x0)
    check(G, inv, ev, t, x) = begin
        if !isnothing(ev)
            x == apply_delta(xprev[], ev) ||
                push!(bad, "$(ev.name): state $x is not the previous state plus Δ")
            made = length(G.nodes) - n0[]
            if ev.type == BIRTH || ev.type == SAMPLE
                made == 1 || push!(bad, "$(ev.name): created $made nodes, expected 1")
                nn = G.nodes[end]
                nn.slate == t || push!(bad, "$(ev.name): slate $(nn.slate) ≠ firing time $t")
                par = G.nodes[findfirst(m -> m.name==nn.parent, G.nodes)]
                par.slate < t || push!(bad, "$(ev.name): parent slate $(par.slate) not < $t")
                nn.name ∈ par.children || push!(bad, "$(ev.name): new node not among parent's children")
                if ev.type == BIRTH
                    nn.type == Node || push!(bad, "$(ev.name): birth node has type $(nn.type)")
                    ## the acted-on lineage (the new node's parent) is gone and
                    ## exactly r[j] copies of the new node are open in each deme j
                    for j ∈ 1:ndeme
                        count(==(nn.name), inv[j]) == ev.r[j] ||
                            push!(bad, "$(ev.name): $(count(==(nn.name), inv[j])) products in deme $j, r says $(ev.r[j])")
                        extra = count(==(nn.parent), inv[j]) - count(==(nn.parent), before[][j]) +
                            (j == ev.from ? 1 : 0)
                        extra == 0 || push!(bad, "$(ev.name): parent lineage not replaced in deme $j")
                    end
                else
                    nn.type == Sample || push!(bad, "$(ev.name): sample node has type $(nn.type)")
                    destructive = ev.r[ev.from] == 0
                    (nn.name ∈ inv[ev.from]) == !destructive ||
                        push!(bad, "$(ev.name): lineage continuation disagrees with the event's production vector")
                end
            else
                made == 0 || push!(bad, "$(ev.name): created $made nodes, expected 0")
            end
            for (i,sym) ∈ enumerate(model.demes)
                length(inv[i]) == getproperty(x,sym) || push!(bad, "inventory/population mismatch in $sym")
            end
        end
        n0[] = length(G.nodes)
        before[] = [copy(inv[d]) for d ∈ 1:ndeme]
        xprev[] = x
    end
    ## The simulator's own loop, with `check` called after every event.
    G = PhyloPOMP._simulate(model, θ, x0, graft, 0.0, tmax, rng, demeset, samplemap, check)
    G, bad
end

function structure_problems(g; sample_max_children::Integer = 1)
    bad = String[]
    byname = Dict(g[i].name => g[i] for i ∈ eachindex(g))
    length(byname) == length(g) || push!(bad, "duplicate node names")
    for i ∈ eachindex(g)
        n = g[i]
        for c ∈ n.children
            haskey(byname,c) || (push!(bad, "child $c of $(n.name) missing"); continue)
            byname[c].parent == n.name || push!(bad, "child $c does not point back to $(n.name)")
            byname[c].slate > n.slate || push!(bad, "non-positive edge $(n.name)->$c")
        end
        if !isnothing(n.parent)
            haskey(byname,n.parent) || push!(bad, "parent of $(n.name) missing")
            haskey(byname,n.parent) && n.name ∉ byname[n.parent].children &&
                push!(bad, "$(n.name) not listed among parent's children")
        end
        deg = length(n.children)
        n.type == Root   && deg != 1 && push!(bad, "Root $(n.name) has degree $deg")
        n.type == Node   && deg != 2 && push!(bad, "Node $(n.name) has degree $deg")
        n.type == Sample && deg > sample_max_children && push!(bad, "Sample $(n.name) has degree $deg")
        isnothing(n.parent) == (n.type == Root) || push!(bad, "Root/parent mismatch at $(n.name)")
        (g.t0 ≤ n.slate ≤ g.time) || push!(bad, "slate $(n.slate) outside [$(g.t0),$(g.time)]")
    end
    issorted([g[i].slate for i ∈ eachindex(g)]) || push!(bad, "nodes not time-sorted")
    bad
end

total_branch_length(g) = sum(
    g[i].slate - g[g[i].parent].slate for i ∈ eachindex(g) if !isnothing(g[i].parent);
    init = 0.0,
)

internal_times(g) = sort([g[i].slate for i ∈ nodes(g)])

function roundtrip_problems(g, D::Module)
    bad = String[]
    g2 = parse_newick(newick(g); demes = D, time = g.time)
    length(g2) == length(g) || push!(bad, "newick: node count $(length(g)) -> $(length(g2))")
    nsample(g2) == nsample(g) || push!(bad, "newick: nsample $(nsample(g)) -> $(nsample(g2))")
    length(roots(g2)) == length(roots(g)) || push!(bad, "newick: root count changed")
    isapprox(internal_times(g2), internal_times(g); atol = 1e-3) ||
        push!(bad, "newick: internal-node times changed")
    isapprox(total_branch_length(g2), total_branch_length(g); rtol = 1e-4) ||
        push!(bad, "newick: total branch length changed")
    Set(g2[i].deme for i ∈ samples(g2)) == Set(g[i].deme for i ∈ samples(g)) ||
        push!(bad, "newick: sample deme set changed")
    if length(roots(g)) == 1 && nsample(g) ≥ 2      # cblv is single-tree only
        xy = cblv(g)
        g3 = parse_cblv(xy...; demes = D, time = g.time)
        nsample(g3) == nsample(g) || push!(bad, "cblv: nsample changed")
        isapprox(internal_times(g3), internal_times(g); atol = 1e-3) ||
            push!(bad, "cblv: internal-node times changed")
        xy2 = cblv(g3)
        (xy2[1] ≈ xy[1] && xy2[2] ≈ xy[2]) || push!(bad, "cblv: vectors not reproduced")
    end
    bad
end

"""
    run_invariants(model, θ, cases; ntree, seed, demeset, samplemap, sample_max_children)

`cases` is a vector of `(x0, graft)` pairs (e.g. one single-root, one forest).
Runs `ntree` seeds per case through all three layers, with `@test`s, and
returns `(nonempty, inline)`: the number of trees with ≥1 sample and the number
of degree-1 (inline) Sample nodes seen, so callers can assert coverage.
"""
function run_invariants(
    model, θ, cases; ntree::Integer = 100, seed::Integer = 20261003,
    demeset::Module = PhyloPOMP.Unstructured, samplemap = nothing,
    sample_max_children::Integer = 1, tmax::Real = 15.0,
)
    rng = MersenneTwister(seed)
    nonempty = 0; inline = 0
    for (x0, graft) ∈ cases, _ ∈ 1:ntree
        s = rand(rng, UInt32)
        g, bad = simulate_checked(model, θ; x0, graft, tmax, rng = MersenneTwister(s), demeset, samplemap)
        @test isempty(bad)
        g0 = simulate(model, θ; x0, graft, tmax, rng = MersenneTwister(s), demeset, samplemap)
        @test newick(g) == newick(g0)   # the check draws nothing, so the tree is the same
        nsample(g0) == 0 && continue
        nonempty += 1
        inline += count(length(g0[i].children) == 1 for i ∈ samples(g0))
        @test isempty(structure_problems(g0; sample_max_children))
        @test isempty(roundtrip_problems(g0, demeset))
        @test g0.t0 == 0.0 && g0.time == tmax
        @test nsample(g0) == length(samples(g0))
    end
    (nonempty = nonempty, inline = inline)
end

end # module
