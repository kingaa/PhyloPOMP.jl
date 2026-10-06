"""
Checks for `simulate` on the SIR model (`src/examples/mgp_sir.jl`).

SIR has no hand-written filter, so the checks are the ones in
`simulate_checks.jl`, a scaling bound, and `prune!` on a few hand-built trees.
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
using PhyloPOMP: Root, Node, Sample, Time, push_node!, prune!, repair!, check_simulated!
using PhyloPOMP: MGPModel, Event, BIRTH, MIGRATION, DEATH, SAMPLE, NEUTRAL, @demes, @mgp
using Main.SimulateChecks

## `Division` is `fork(A => B, B)`: the parent leaves its deme.
## `Twins` is `fork(I => E, E)`.
@mgp Division begin
    compartments = (A, B)
    demes = (A, B)
    params = (λ, μ)
    @event divide  rate=λ*A pop=(A=-1, B=+2) move=fork(A => B, B)
    @event sampleB rate=μ*B pop=(B=-1)       move=sample_remove(B) kind=singular
end
## `Leak`: a constant-rate drain of a non-lineage compartment; its hazard does
## not vanish when S is empty, so it must be rejected at run time.
@mgp Leak begin
    compartments = (S, I)
    demes = (I,)
    params = (ρ,)
    @event drain rate=ρ pop=(S=-1) move=none
end
@mgp Twins begin
    compartments = (S, E, I, R)
    demes = (E, I)
    params = (β, σ, γ, ψ, N)
    @event infection   rate=β*S*I/N pop=(S=-1, E=+2, I=-1) move=fork(I => E, E)
    @event progression rate=σ*E     pop=(E=-1, I=+1)       move=swap(E => I)
    @event recovery    rate=γ*I     pop=(I=-1, R=+1)       move=chop(I)
    @event sampling    rate=ψ*I     pop=()                 move=sample(I) kind=singular
end

## is node `i` in the subtree rooted at `root`?
function _descends(g, i, root)
    while !isnothing(i)
        i == root && return true
        i = g[i].parent
    end
    false
end

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

    @info h2("@mgp rejects events whose `move` and `pop` disagree on a deme")
    ## Destructive sampling is decided from `r`. The macro checks that `pop` agrees.
    bad(ev) = :(PhyloPOMP.@mgp Bad begin
        compartments = (S, E, I, R); demes = (E, I); params = (β, σ, γ, ψ, N)
        $ev
    end)
    ok(ev) = macroexpand(@__MODULE__, bad(ev)) isa Expr
    @test ok(:(@event s rate=ψ*I pop=()      move=sample(I)        kind=singular))
    @test ok(:(@event s rate=ψ*I pop=(I=-1)  move=sample_remove(I) kind=singular))
    @test ok(:(@event f rate=β*S*I/N pop=(S=-1, E=+1) move=fork(I => E, I)))
    @test ok(:(@event f rate=β*S*I/N pop=(S=-1, I=+1) move=fork(I => I, I)))
    @test ok(:(@event m rate=σ*E pop=(E=-1, I=+1) move=swap(E => I)))
    @test ok(:(@event d rate=γ*I pop=(I=-1, R=+1) move=chop(I)))
    @test_throws ErrorException ok(:(@event s rate=ψ*I pop=(I=-1) move=sample(I) kind=singular))  # continues, yet pop drops I
    @test_throws ErrorException ok(:(@event s rate=ψ*I pop=()     move=sample_remove(I) kind=singular))
    @test_throws ErrorException ok(:(@event f rate=β*S*I/N pop=(S=-1) move=fork(I => E, I)))      # missing E=+1
    @test_throws ErrorException ok(:(@event m rate=σ*E pop=(I=+1) move=swap(E => I)))             # missing E=-1
    @test_throws ErrorException ok(:(@event d rate=γ*I pop=(R=+1) move=chop(I)))                  # missing I=-1
    @test_throws ErrorException ok(:(@event n rate=γ*I pop=(I=-1) move=none))                     # none cannot change a deme
    @test_throws ErrorException ok(:(@event n rate=γ*I pop=(Q=-1) move=none))                     # undeclared compartment
    @test_throws ErrorException ok(:(@event n rate=γ pop=(S=-1, S=+1) move=none))                 # repeated compartment
    ## births are binary: exactly two products
    @test ok(:(@event f rate=β*S*I/N pop=(S=-1, E=+2, I=-1) move=fork(I => E, E)))              # parent leaves its deme
    @test_throws ErrorException ok(:(@event f rate=β*S*I/N pop=(S=-1, E=+1, I=-1) move=fork(I => E)))
    @test_throws ErrorException ok(:(@event f rate=β*S*I/N pop=(S=-1, E=+1, I=+1) move=fork(I => E, I, I)))

    @info h2("birth forms: parent leaves its deme; both products in another deme")
    ## A -> B + B. After k divisions and s samples, A = A0 - k and B = 2k - s.
    let θd = (λ = 1.0, μ = 0.7)
        covd = run_invariants(Division, θd, [((A=5, B=0), [5, 0]), ((A=2, B=3), [2, 3])];
                              ntree=100, tmax=10.0, sample_max_children=0)
        @test covd.nonempty > 100
        tr = simulate_trajectory(Division, θd; x0=(A=5, B=0), graft=[5, 0], tmax=10.0,
                                 rng=MersenneTwister(3))
        k = count(==(:divide), tr.events); s = count(==(:sampleB), tr.events)
        @test tr.states[end] == (A = 5 - k, B = 2k - s)
        @test k == 5                                      # every A eventually divides
        @test nsample(tr.genealogy) == s
    end
    let θt = (β = 3.0, σ = 1.0, γ = 1.0, ψ = 0.3, N = 100.0)
        covt = run_invariants(Twins, θt, [((S=98, E=1, I=1, R=0), [1, 1]), ((S=95, E=2, I=3, R=0), [2, 3])];
                              ntree=100, tmax=15.0, sample_max_children=1)
        @test covt.nonempty > 100
        @test covt.inline > 0
    end

    @info h2("argument validation")
    @test_throws ArgumentError simulate(PhyloPOMP.SIR, θ; x0=x0, graft=[1], tmax=-1.0)
    @test_throws ArgumentError simulate(PhyloPOMP.SIR, θ; x0=(S=99, I=1), graft=[1], tmax=10.0)   # R missing
    ## a hazard that does not vanish when the compartment it depletes is empty
    @test_throws AssertionError simulate(Leak, (ρ = 1.0,); x0=(S=2, I=1), graft=[1], tmax=100.0)
    ## The loop's own check fires first on any compiled model, so call `apply_event!` directly.
    let G = Genealogy{PhyloPOMP.Unstructured}(Time(0.0)), inv = PhyloPOMP.SimInventory(1)
        push_node!(G, 0.0, Root, nothing)
        ev = PhyloPOMP.SIR.events[findfirst(e -> e.name == :recovery, PhyloPOMP.SIR.events)]
        err = try PhyloPOMP.apply_event!(G, inv, ev, Time(1.0), MersenneTwister(1), nothing); nothing
              catch e; e end
        @test err isa AssertionError && occursin("recovery", err.msg) && occursin("deme 1", err.msg)
    end
    ## prune! precondition: names 1:n with parents before children
    let G = Genealogy{PhyloPOMP.Unstructured}(Time(0.0))
        r = push_node!(G, 0.0, Root, nothing)
        c = push_node!(G, 1.0, Sample, r); push!(G[r].children, c)
        G.nodes[1], G.nodes[2] = G.nodes[2], G.nodes[1]          # child stored before parent
        @test_throws ArgumentError prune!(G)
    end
    @test_throws ArgumentError simulate(PhyloPOMP.SIR, θ; x0=x0, graft=[0], tmax=10.0)
    @test_throws ArgumentError simulate(PhyloPOMP.SIR, θ; x0=x0, graft=[1,0], tmax=10.0)
    ## non-finite hazards are rejected (an Inf hazard gives dt = 0 forever)
    @test_throws AssertionError simulate(PhyloPOMP.SIR, merge(θ, (N=0.0,)); x0=x0, graft=[1], tmax=10.0)
    @test_throws AssertionError simulate(PhyloPOMP.SIR, merge(θ, (β=NaN,)); x0=x0, graft=[1], tmax=10.0)
    @test_throws AssertionError simulate(PhyloPOMP.SIR, merge(θ, (γ=-1.0,)); x0=x0, graft=[1], tmax=10.0)

    @info h2("per-event timing, structure, round trips (single root and forest)")
    cov = run_invariants(PhyloPOMP.SIR, θ, [(x0, [1]), ((S=96, I=4, R=0), [4])];
                         ntree=100, tmax=15.0, sample_max_children=1)
    @test cov.nonempty > 100
    ## non-destructive sampling must produce degree-1 (inline) samples
    @test cov.inline > 0

    @info h2("simulate_trajectory: same tree as simulate, consistent population path")
    let seed = 20261005
        gA = simulate(PhyloPOMP.SIR, θ; x0=(S=96,I=4,R=0), graft=[4], tmax=15.0, rng=MersenneTwister(seed))
        tr = simulate_trajectory(PhyloPOMP.SIR, θ; x0=(S=96,I=4,R=0), graft=[4], tmax=15.0, rng=MersenneTwister(seed))
        @test newick(tr.genealogy) == newick(gA)            # identical random draws
        @test tr.times[1] == 0.0 && tr.states[1] == (S=96,I=4,R=0)
        @test length(tr.states) == length(tr.times) == length(tr.events) + 1
        @test issorted(tr.times) && allunique(tr.times) && tr.times[end] < 15.0
        byname = Dict(ev.name => ev for ev ∈ PhyloPOMP.SIR.events)
        @test all(tr.states[i+1] == PhyloPOMP.apply_delta(tr.states[i], byname[tr.events[i]])
                  for i ∈ eachindex(tr.events))               # each step is exactly one event's Δ
        @test all(x.S + x.I + x.R == 100 for x ∈ tr.states)   # closed population
        @test all(x.I ≥ 0 && x.S ≥ 0 for x ∈ tr.states)
        ## every ψ-sample is kept by prune!, so the tree's sample count is the
        ## trajectory's sampling-event count; the final I bounds the lineages
        ## that could still be open at tmax
        @test nsample(tr.genealogy) == count(==(:sampling), tr.events)
        @test tr.states[end].I ≥ 0
        @test length(tr.events) > 10
    end

    @info h2("distribution: rates, event choice, lineage choice (β = 0, exact laws)")
    ## β = 0, so the I individuals are independent. Each lives Exp(γ) and is
    ## sampled at rate ψ while alive:
    ##   E[nsample] = ψ I0 / γ,  variance I0 (ψ/γ + ψ²/γ²)
    ##   P(recovery) = γ / (γ + ψ)
    ##   recovery times ~ Exp(γ)
    ## Bounds are 4 to 5 standard errors.
    let γ = 1.0, ψ = 0.5, I0 = 50, nseed = 200
        θ0 = (β = 0.0, γ = γ, ψ = ψ, N = 100.0)
        ns = Int[]; nrec = 0; nev = 0; rectimes = Float64[]
        for seed ∈ 1:nseed
            tr = simulate_trajectory(PhyloPOMP.SIR, θ0; x0=(S=50, I=I0, R=0), graft=[I0], tmax=40.0,
                                     rng=MersenneTwister(seed))
            push!(ns, nsample(tr.genealogy))
            nrec += count(==(:recovery), tr.events); nev += length(tr.events)
            append!(rectimes, tr.times[i+1] for i ∈ eachindex(tr.events) if tr.events[i] == :recovery)
            @test tr.states[end].I == 0                   # tmax = 40 ≈ ∞ for Exp(1) lifetimes
        end
        μn = ψ * I0 / γ; sen = sqrt(I0 * (ψ/γ + (ψ/γ)^2) / nseed)
        @test abs(sum(ns)/nseed - μn) < 4sen
        pr = γ / (γ + ψ); sep = sqrt(pr * (1 - pr) / nev)
        @test abs(nrec/nev - pr) < 5sep
        ## one-sample KS against Exp(γ)
        sort!(rectimes); m = length(rectimes)
        D = maximum(max((i/m) - (1 - exp(-γ*rectimes[i])), (1 - exp(-γ*rectimes[i])) - (i-1)/m) for i ∈ 1:m)
        @test m == I0 * nseed
        @test D < sqrt(-log(5e-5) / (2m))                 # Kolmogorov bound at α = 1e-4
    end
    ## Two founders, same rates. Samples should split evenly between the roots.
    let nseed = 400, θ0 = (β = 0.0, γ = 1.0, ψ = 2.0, N = 100.0)
        n1 = 0; n2 = 0
        for seed ∈ 1:nseed
            g = simulate(PhyloPOMP.SIR, θ0; x0=(S=98, I=2, R=0), graft=[2], tmax=40.0,
                         rng=MersenneTwister(seed))
            r = roots(g)
            length(r) == 2 || continue                    # a founder with no samples is dropped
            n1 += count(_descends(g, i, r[1]) for i ∈ samples(g))
            n2 += count(_descends(g, i, r[2]) for i ∈ samples(g))
        end
        ## per founder: mean 2, var 2 + 4 = 6 samples; difference sd = sqrt(2 * 6 * nseed)
        @test abs(n1 - n2) < 4 * sqrt(12 * nseed)
        @test n1 + n2 > 1000
    end

    @info h2("distribution: mean prevalence follows the SIR ODE for a large population")
    ## N = 10,000. Mean I(t) over seeds should track the SIR ODE through the peak.
    let N = 10_000, I0 = 100, β = 2.0, γ = 1.0, nseed = 60, tcheck = (1.0, 2.0, 3.0, 4.0)
        θo = (β = β, γ = γ, ψ = 0.0, N = Float64(N))
        acc = zeros(length(tcheck))
        for seed ∈ 1:nseed
            tr = simulate_trajectory(PhyloPOMP.SIR, θo; x0=(S=N-I0, I=I0, R=0), graft=[I0], tmax=4.01,
                                     rng=MersenneTwister(seed))
            for (j, tc) ∈ enumerate(tcheck)
                acc[j] += tr.states[searchsortedlast(tr.times, tc)].I
            end
        end
        ## RK4 on dS/dt = -βSI/N, dI/dt = βSI/N - γI
        f(u) = (-β*u[1]*u[2]/N, β*u[1]*u[2]/N - γ*u[2])
        u = (Float64(N - I0), Float64(I0)); h = 1e-3; tt = 0.0; ode = Float64[]
        for tc ∈ tcheck
            while tt < tc - 1e-12
                k1 = f(u); k2 = f(u .+ h/2 .* k1); k3 = f(u .+ h/2 .* k2); k4 = f(u .+ h .* k3)
                u = u .+ h/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4); tt += h
            end
            push!(ode, u[2])
        end
        for j ∈ eachindex(tcheck)
            @test abs(acc[j]/nseed - ode[j]) / ode[j] < 0.05
        end
        @test ode[3] > 1000                               # the check spans the epidemic peak
    end

    @info h2("scaling: a full N=30,000 epidemic must not be quadratic in nodes")
    ## N = 100,000 takes about 0.1 s and N = 300,000 about 0.35 s.
    ## A search on the event path makes this quadratic, and 5 s catches that.
    let N = 30_000
        θb = (β=2.0, γ=1.0, ψ=0.05, N=Float64(N))
        xb = (S=N-10, I=10, R=0)
        simulate(PhyloPOMP.SIR, θb; x0=xb, graft=[10], tmax=1.0, rng=MersenneTwister(1)) # compile
        tsec = @elapsed gb = simulate(PhyloPOMP.SIR, θb; x0=xb, graft=[10], tmax=40.0, rng=MersenneTwister(7))
        @test nsample(gb) > 500
        @test tsec < 5.0
    end

    @info h2("simulate: a hand-built model is checked like an @mgp model")
    ## Each model breaks one rule that `@mgp` would have caught.
    bad = [
        (MGPModel(:HalfBirth, [:S,:I], [:I],
            [Event(:b, [:S=>-1], (x,θ)->0.05*x.I*x.S, [1], BIRTH, 1, Int[], true, false)]),
         (S=10, I=1)),
        (MGPModel(:TypoDelta, [:S,:I,:R], [:I],
            [Event(:rec, [:I=>-1, :Rr=>+1], (x,θ)->1.0*x.I, [0], DEATH, 1, Int[], true, false)]),
         (S=0, I=3, R=0)),
        (MGPModel(:DoubleSample, [:I], [:I],
            [Event(:s, Pair{Symbol,Int}[], (x,θ)->1.0*x.I, [2], SAMPLE, 1, Int[], false, true)]),
         (I=2,)),
        (MGPModel(:BadFrom, [:I], [:I],
            [Event(:d, [:I=>-1], (x,θ)->1.0*x.I, [0], DEATH, 2, Int[], true, false)]),
         (I=2,)),
    ]
    for (m, x0) ∈ bad
        @test_throws ArgumentError simulate(m, (;); x0, graft=[x0.I], tmax=5.0, rng=MersenneTwister(1))
    end
    θref = (β=4.0, σ=1.0, γ=1.0, ω=1.0, ψ=0.3, χ=0.0, N=100.0)
    G = simulate(PhyloPOMP.SEIR_REFERENCE, θref; x0=(S=99, E=0, I=1, R=0), graft=[0,1],
                 tmax=5.0, rng=MersenneTwister(2))
    @test isempty(structure_problems(G))

    @info h2("check_simulated!: always-on post-condition rejects malformed output")
    ## Each of these is a tree `check_simulated!` should reject.
    mk() = begin
        G = Genealogy{PhyloPOMP.Unstructured}(Time(0.0))
        r = push_node!(G, 0.0, Root, nothing)
        a = push_node!(G, 1.0, Node, r);   push!(G[r].children, a)
        s = push_node!(G, 2.0, Sample, a); push!(G[a].children, s)
        u = push_node!(G, 3.0, Sample, a); push!(G[a].children, u)
        G.time = 5.0
        G, r, a, s, u
    end
    let (G, _, _, _, _) = mk()                       # well-formed: passes
        repair!(G); @test isnothing(check_simulated!(G))
    end
    let (G, _, a, s, _) = mk()                       # Sample stamped before its parent
        G[s].slate = 0.5
        repair!(G); @test_throws AssertionError check_simulated!(G)
    end
    let (G, _, a, _, u) = mk()                       # degree-1 internal Node
        filter!(!=(u), G[a].children); G[u].parent = G[a].parent
        push!(G[G[a].parent].children, u)
        repair!(G); @test_throws AssertionError check_simulated!(G)
    end
    let (G, r, a, s, _) = mk()                       # stale parent pointer
        G[s].parent = r
        repair!(G); @test_throws AssertionError check_simulated!(G)
    end
    let (G, _, _, _, _) = mk()                       # Sample with two children
        G[2].type = Sample
        repair!(G); @test_throws AssertionError check_simulated!(G)
    end
    let (G, _, _, _, _) = mk()                       # nsample bookkeeping wrong
        repair!(G); G.nsample = 1
        @test_throws AssertionError check_simulated!(G)
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
