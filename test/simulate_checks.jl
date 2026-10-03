"""
Shared checks for the generic forward simulator (`src/simulate.jl`), used by
`seir_simulate.jl`, `mers_simulate.jl`, and `sir_simulate.jl`.

Why these exist: the simulate-then-filter "finite logLik" tests cannot see
internal-node timing (pulling every internal node 90% of the way toward its
parent still gives a finite likelihood), and the two timing bugs found in
`check_milestone3.md`/`check_milestone4.md` were therefore invisible to them.
The checks below were run against the pre-fix simulator (commit `e5d6b30`)
and fail on it; they pass on the current code.

Three layers:

1. `simulate_checked` replays `simulate`'s loop using the package's own
   `apply_event!` and asserts the timing contract after EVERY event: BIRTH and
   SAMPLE append exactly one node, stamped with the firing time, whose parent
   is strictly earlier and on which every lineage the event produced now sits;
   MIGRATION/DEATH/NEUTRAL append none. Its output must equal `simulate`'s for
   the same seed, which is also checked.
2. `structure_problems` inspects a finished genealogy: parent strictly earlier
   than child, no zero-length edges, Root degree 1 / Node degree 2 / Sample
   degree ≤ `sample_max_children`, consistent parent/child links, time-sorted.
3. `roundtrip_problems` requires Newick and CBLV round trips to preserve the
   internal-node times and total branch length, not just node counts.

`run_invariants` drives all three over many seeds (single root and forest).
"""
module SimulateChecks

using Test
using Random: MersenneTwister
using PhyloPOMP
using PhyloPOMP: Root, Node, Sample, Time, SimInventory, push_node!, add!,
    apply_delta, apply_event!, prune!, repair!, rcateg, BIRTH, SAMPLE

export simulate_checked, structure_problems, roundtrip_problems,
    total_branch_length, run_invariants

function simulate_checked(
    model, θ; x0, graft, tmax, rng,
    demeset::Module = PhyloPOMP.Unstructured, samplemap = nothing,
)
    ndeme = length(model.demes)
    G = Genealogy{demeset}(Time(0.0))
    inv = SimInventory(ndeme)
    for (d,n) ∈ enumerate(graft), _ ∈ 1:n
        r = push_node!(G,Time(0.0),Root,nothing)
        c = push_node!(G,Time(0.0),Node,r)
        push!(G[r].children,c)
        add!(inv,d,c)
    end
    x = x0
    t = Time(0.0)
    haz = zeros(length(model.events))
    bad = String[]
    while true
        for (k,ev) ∈ enumerate(model.events)
            haz[k] = ev.hazard(x,θ)
        end
        total = sum(haz)
        total ≤ 0 && break
        dt = -log(rand(rng))/total
        t + dt ≥ tmax && break
        t += dt
        k, _ = rcateg(haz; rng)
        ev = model.events[k]
        x = apply_delta(x,ev)
        n0 = length(G.nodes)
        apply_event!(G,inv,ev,model,t,rng,samplemap)
        made = length(G.nodes) - n0
        if ev.type == BIRTH || ev.type == SAMPLE
            made == 1 || push!(bad, "$(ev.name): created $made nodes, expected 1")
            nn = G.nodes[end]
            nn.slate == t || push!(bad, "$(ev.name): slate $(nn.slate) ≠ firing time $t")
            par = G.nodes[findfirst(m -> m.name==nn.parent, G.nodes)]
            par.slate < t || push!(bad, "$(ev.name): parent slate $(par.slate) not < $t")
            nn.name ∈ par.children || push!(bad, "$(ev.name): new node not among parent's children")
            if ev.type == BIRTH
                nn.type == Node || push!(bad, "$(ev.name): birth node has type $(nn.type)")
                ## the continuing lineage and every new lineage sit on the new node
                nn.name ∈ inv[ev.from] || push!(bad, "$(ev.name): continuing lineage not on new node")
                for j ∈ ev.into
                    nn.name ∈ inv[j] || push!(bad, "$(ev.name): new lineage in deme $j not on new node")
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
    G.time = Time(tmax)
    prune!(G); repair!(G)
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
        @test newick(g) == newick(g0)
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
