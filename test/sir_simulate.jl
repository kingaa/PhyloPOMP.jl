"""
Tests for the generic forward simulator (`src/simulate.jl`) on the SIR model
declared in `src/examples/mgp_sir.jl`.

SIR has no hand-coded filter, so there is no simulate-then-filter check here.
Instead this file leans on the shared invariants in `simulate_checks.jl`
(per-event timing, structure, time-preserving round trips, forests), adds a
scaling check, and tests `prune!` on hand-built genealogies.
"""
module SIRSimulateTest

## `h1`/`h2` come from runtests.jl; fall back when run standalone.
const h1 = isdefined(Main, :h1) ? Main.h1 : identity
const h2 = isdefined(Main, :h2) ? Main.h2 : identity
isdefined(Main, :SimulateChecks) || Base.include(Main, joinpath(@__DIR__, "simulate_checks.jl"))

@info h1("SIR forward simulator")

using Test
using Random: MersenneTwister
using PhyloPOMP
using PhyloPOMP: Root, Node, Sample, Time, push_node!, prune!, repair!
using Main.SimulateChecks

@testset verbose=true "SIR forward simulator" begin

    θ = (β=3.0, γ=1.0, ψ=0.3, N=100.0)
    x0 = (S=99, I=1, R=0)

    @info h2("model table")
    @test PhyloPOMP.SIR isa PhyloPOMP.MGPModel
    @test PhyloPOMP.SIR.demes == [:I]
    @test [ev.type for ev ∈ PhyloPOMP.SIR.events] ==
        [PhyloPOMP.BIRTH, PhyloPOMP.DEATH, PhyloPOMP.SAMPLE]
    ## sampling is non-destructive: production 1 in deme I, no negative Δ on I
    @test PhyloPOMP.SIR.events[3].r == [1]
    @test isempty(PhyloPOMP.SIR.events[3].Δ)

    @info h2("argument validation")
    @test_throws ArgumentError simulate(PhyloPOMP.SIR, θ; x0=x0, graft=[0], tmax=10.0)
    @test_throws ArgumentError simulate(PhyloPOMP.SIR, θ; x0=x0, graft=[1,0], tmax=10.0)

    @info h2("per-event timing, structure, round trips (single root and forest)")
    cov = run_invariants(PhyloPOMP.SIR, θ, [(x0, [1]), ((S=96, I=4, R=0), [4])];
                         ntree=100, tmax=15.0, sample_max_children=1)
    @test cov.nonempty > 100
    ## non-destructive sampling must produce degree-1 (inline) samples
    @test cov.inline > 0

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

    @info h2("prune!: forest whose first founding lineage dies unsampled")
    G = Genealogy{PhyloPOMP.Unstructured}(Time(0.0))
    r1 = push_node!(G, 0.0, Root, nothing)
    c1 = push_node!(G, 0.0, Node, r1);   push!(G[r1].children, c1)  # dies, never sampled
    r2 = push_node!(G, 0.0, Root, nothing)
    c2 = push_node!(G, 0.0, Node, r2);   push!(G[r2].children, c2)
    s  = push_node!(G, 1.5, Sample, c2); push!(G[c2].children, s)
    G.time = 5.0
    prune!(G); repair!(G)
    @test length(roots(G)) == 1
    @test [G[i].type for i ∈ 1:2] == [Root, Sample] && G[2].slate == 1.5
end

end
