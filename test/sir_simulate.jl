"""
Tests for the generic forward simulator (`src/simulate.jl`) on the SIR model
declared in `src/examples/mgp_sir.jl`.

SIR has no hand-coded filter, so there is no simulate-then-filter check here.
The finite-logLik checks in `seir_simulate.jl`/`mers_simulate.jl` are blind to
internal-node timing (see `check_milestone4.md`), so this file instead tests the
properties directly:

  * every BIRTH/SAMPLE event creates exactly one node, stamped with the
    firing time, whose parent is strictly earlier; MIGRATION/DEATH/NEUTRAL
    create none (a driver that replays the loop of `simulate`);
  * every finished tree has parent < child in time, no zero-length edges,
    Root degree 1, internal Node degree 2, Sample degree 0 or 1, consistent
    parent/child links;
  * Newick and CBLV round trips preserve internal-node times and total
    branch length, not just counts;
  * `prune!` on hand-built genealogies gives the expected tree.
"""
module SIRSimulateTest

## `h1`/`h2` come from runtests.jl; fall back when run standalone.
const h1 = isdefined(Main, :h1) ? Main.h1 : identity
const h2 = isdefined(Main, :h2) ? Main.h2 : identity

@info h1("SIR forward simulator")

using Test
using Random: MersenneTwister
using PhyloPOMP
using PhyloPOMP: Root, Node, Sample, Time, SimInventory, push_node!, add!,
    apply_delta, apply_event!, prune!, repair!, rcateg, BIRTH, SAMPLE

## Replays the loop in `simulate` (src/simulate.jl) and checks the timing
## contract after every event. Returns the finished genealogy and a vector of
## violation messages.
function simulate_checked(model, θ; x0, graft, tmax, rng)
    ndeme = length(model.demes)
    G = Genealogy{PhyloPOMP.Unstructured}(Time(0.0))
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
        apply_event!(G,inv,ev,model,t,rng,nothing)
        made = length(G.nodes) - n0
        if ev.type == BIRTH || ev.type == SAMPLE
            made == 1 || push!(bad, "$(ev.name): created $made nodes")
            nn = G.nodes[end]
            nn.slate == t || push!(bad, "$(ev.name): slate $(nn.slate) ≠ firing time $t")
            par = G.nodes[findfirst(m -> m.name==nn.parent, G.nodes)]
            par.slate < t || push!(bad, "$(ev.name): parent not strictly earlier")
        else
            made == 0 || push!(bad, "$(ev.name): created $made nodes")
        end
        for (i,sym) ∈ enumerate(model.demes)
            length(inv[i]) == getproperty(x,sym) || push!(bad, "inventory mismatch")
        end
    end
    G.time = Time(tmax)
    prune!(G); repair!(G)
    G, bad
end

## Structural problems in a finished genealogy.
function structure_problems(g)
    bad = String[]
    byname = Dict(g[i].name => g[i] for i ∈ eachindex(g))
    for i ∈ eachindex(g)
        n = g[i]
        for c ∈ n.children
            byname[c].parent == n.name || push!(bad, "child $c does not point to $(n.name)")
            byname[c].slate > n.slate || push!(bad, "non-positive edge $(n.name)->$c")
        end
        deg = length(n.children)
        n.type == Root   && deg != 1 && push!(bad, "Root $(n.name) degree $deg")
        n.type == Node   && deg != 2 && push!(bad, "Node $(n.name) degree $deg")
        n.type == Sample && deg > 1  && push!(bad, "Sample $(n.name) degree $deg")
        isnothing(n.parent) == (n.type == Root) || push!(bad, "Root/parent mismatch at $(n.name)")
    end
    issorted([g[i].slate for i ∈ eachindex(g)]) || push!(bad, "nodes not time-sorted")
    bad
end

total_branch_length(g) = sum(
    g[i].slate - g[g[i].parent].slate for i ∈ eachindex(g) if !isnothing(g[i].parent);
    init = 0.0,
)

@testset verbose=true "SIR forward simulator" begin

    θ = (β=3.0, γ=1.0, ψ=0.3, N=100.0)
    x0 = (S=99, I=1, R=0)

    @info h2("model table")
    @test PhyloPOMP.SIR isa PhyloPOMP.MGPModel
    @test PhyloPOMP.SIR.demes == [:I]
    @test [ev.type for ev ∈ PhyloPOMP.SIR.events] ==
        [PhyloPOMP.BIRTH, PhyloPOMP.DEATH, PhyloPOMP.SAMPLE]
    ## sampling is non-destructive: no negative Δ on its own deme
    @test isempty(PhyloPOMP.SIR.events[3].Δ)

    @info h2("argument validation")
    @test_throws ArgumentError simulate(PhyloPOMP.SIR, θ; x0=x0, graft=[0], tmax=10.0)
    @test_throws ArgumentError simulate(PhyloPOMP.SIR, θ; x0=x0, graft=[1,0], tmax=10.0)

    @info h2("per-event timing and structure over many trees")
    rng = MersenneTwister(20261003)
    ntree = 0; nbirth_trees = 0; ninline = 0
    for (graft, x) ∈ (([1], x0), ([4], (S=96, I=4, R=0)))
        for _ ∈ 1:100
            seed = rand(rng, UInt32)
            g, bad = simulate_checked(PhyloPOMP.SIR, θ; x0=x, graft=graft, tmax=15.0,
                                      rng=MersenneTwister(seed))
            @test isempty(bad)
            ## the driver must reproduce `simulate` exactly
            g0 = simulate(PhyloPOMP.SIR, θ; x0=x, graft=graft, tmax=15.0,
                          rng=MersenneTwister(seed))
            @test newick(g) == newick(g0)
            nsample(g0) == 0 && continue
            ntree += 1
            @test isempty(structure_problems(g0))
            @test g0.t0 == 0.0 && g0.time == 15.0
            @test nsample(g0) == length(samples(g0))
            ninline += count(length(g0[i].children)==1 for i ∈ samples(g0))
            length(nodes(g0)) > 0 && (nbirth_trees += 1)

            ## round trips preserve times, not only counts
            g2 = parse_newick(newick(g0); time=g0.time)
            @test length(g2) == length(g0)
            @test sort([g2[i].slate for i ∈ nodes(g2)]) ≈
                sort([g0[i].slate for i ∈ nodes(g0)]) atol=1e-3
            @test total_branch_length(g2) ≈ total_branch_length(g0) rtol=1e-4
            if length(roots(g0)) == 1 && nsample(g0) ≥ 2
                xy = cblv(g0)
                g3 = parse_cblv(xy...; time=g0.time)
                @test nsample(g3) == nsample(g0)
                @test sort([g3[i].slate for i ∈ nodes(g3)]) ≈
                    sort([g0[i].slate for i ∈ nodes(g0)]) atol=1e-3
            end
        end
    end
    @test ntree > 100
    @test nbirth_trees > 50
    ## non-destructive sampling must produce degree-1 (inline) samples
    @test ninline > 0

    @info h2("scaling: a full N=30,000 epidemic must not be quadratic in nodes")
    ## Before the O(1) node lookup in `apply_event!` this took 17.8 s
    ## (`findfirst` over all nodes at every event); after, 0.14 s. The bound is
    ## loose on purpose so slow CI machines do not fail it, but a return to
    ## quadratic behaviour (~100x slower) will.
    let N = 30_000
        θb = (β=2.0, γ=1.0, ψ=0.05, N=Float64(N))
        xb = (S=N-10, I=10, R=0)
        simulate(PhyloPOMP.SIR, θb; x0=xb, graft=[10], tmax=1.0, rng=MersenneTwister(1)) # compile
        tsec = @elapsed gb = simulate(PhyloPOMP.SIR, θb; x0=xb, graft=[10], tmax=40.0, rng=MersenneTwister(7))
        @test nsample(gb) > 500
        @test tsec < 5.0
    end

    @info h2("prune!: unsampled side branch + lineage alive at tmax")
    G = Genealogy{PhyloPOMP.Unstructured}(Time(0.0))
    r  = push_node!(G, 0.0, Root, nothing)
    c  = push_node!(G, 0.0, Node, r);    push!(G[r].children, c)
    p1 = push_node!(G, 1.0, Node, c);    push!(G[c].children, p1)   # birth; side lineage dies unsampled
    p2 = push_node!(G, 2.0, Node, p1);   push!(G[p1].children, p2)  # birth
    s2 = push_node!(G, 2.5, Sample, p2); push!(G[p2].children, s2)  # sampled, alive at tmax
    s1 = push_node!(G, 3.0, Sample, p2); push!(G[p2].children, s1)  # sampled, alive at tmax
    G.time = 5.0
    prune!(G); repair!(G)
    @test length(G) == 4
    @test [G[i].type for i ∈ 1:4] == [Root, Node, Sample, Sample]
    @test [G[i].slate for i ∈ 1:4] == [0.0, 2.0, 2.5, 3.0]
    @test G[1].children == [2] && G[2].parent == 1
    @test G[2].children == [3,4] && G[3].parent == 2 && G[4].parent == 2
    @test isempty(structure_problems(G))

    @info h2("prune!: unsampled sub-clade and unsampled extant lineage")
    G = Genealogy{PhyloPOMP.Unstructured}(Time(0.0))
    r  = push_node!(G, 0.0, Root, nothing)
    c  = push_node!(G, 0.0, Node, r);    push!(G[r].children, c)
    p1 = push_node!(G, 1.0, Node, c);    push!(G[c].children, p1)
    q1 = push_node!(G, 1.2, Node, p1);   push!(G[p1].children, q1)  # side clade, never sampled
    q2 = push_node!(G, 1.4, Node, q1);   push!(G[q1].children, q2)
    p2 = push_node!(G, 2.0, Node, p1);   push!(G[p1].children, p2)  # other branch, never sampled
    s1 = push_node!(G, 3.0, Sample, p2); push!(G[p2].children, s1)
    G.time = 5.0
    prune!(G); repair!(G)
    @test length(G) == 2
    @test [G[i].type for i ∈ 1:2] == [Root, Sample]
    @test G[2].slate == 3.0 && G[2].parent == 1 && G[1].children == [2]
end

end
