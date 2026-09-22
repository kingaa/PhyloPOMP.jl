"""
Acceptance test for the generic forward simulator (`src/simulate.jl`).

The headline check: simulate a genealogy from the `@mgp`-generated `SEIR`
event table (`src/examples/mgp.jl`), then feed the result straight into the
already-tested, hand-coded `NaiveSEIR.filter_pomp` (`src/examples/seir_naive.jl`)
and require a finite log-likelihood. These are two independent encodings of
the same SEIR rate structure; agreement is far stronger evidence of
correctness than a Newick/CBLV round trip alone (which would faithfully
serialize a structurally wrong tree). Round-trip checks are included too,
but only as supporting, non-sufficient evidence.

A second, independently-written filter over the same rate structure --
`GuidedSEIR.filter_pomp(gen, m; ...)` (`src/examples/seir_guided.jl`) -- is
exercised at the end of the file. It is the stricter of the two: it runs
`GuidedSEIR.check` on the genealogy (root degree 1, internal-node degree 2,
sample degree < 2) before building the pomp object, and it needs a guiding
finite-state Markov process (`fsmarkov`) as well as the genealogy.
"""
module SEIRSimulateTest

import ..Main: h1, h2

@info h1("SEIR forward simulator")

using Test
using Random: MersenneTwister, seed!
using PhyloPOMP
using PhyloPOMP: Sample
using PhyloPOMP.NaiveSEIR
using PhyloPOMP.GuidedSEIR
using PhyloPOMP.GuidedSEIR.Demes: Expos, Infec
import PartiallyObservedMarkovProcesses as POMP

@testset verbose=true "SEIR forward simulator" begin

    ## parameters shared between the simulator and the filter -- the point
    ## of this test is that these are the SAME numbers, not fitted/matched
    ## after the fact.
    β, σ, γ, ω, ψ, χ = 4.0, 1.0, 1.0, 1.0, 0.10, 0.0
    pop = 100
    θ = (β=β, σ=σ, γ=γ, ω=ω, ψ=ψ, χ=χ, N=Float64(pop))
    x0 = (S=pop-1, E=0, I=1, R=0)

    @info h2("argument validation")
    @test_throws ArgumentError simulate(
        PhyloPOMP.SEIR, θ; x0=x0, graft=[0,2], tmax=10.0,
    )
    @test_throws ArgumentError simulate(
        PhyloPOMP.SEIR, θ; x0=x0, graft=[0,1,0], tmax=10.0,
    )

    @info h2("simulate a non-degenerate genealogy")
    rng = MersenneTwister(20260813)
    local g
    ok = false
    for _ ∈ 1:100
        g = simulate(PhyloPOMP.SEIR, θ; x0=x0, graft=[0,1], tmax=20.0, rng=rng)
        if nsample(g) ≥ 5
            ok = true
            break
        end
    end
    @test ok
    @test g isa Genealogy{PhyloPOMP.Unstructured}
    @test length(roots(g)) == 1
    @test nsample(g) == length(samples(g))
    @test all(isempty(g[i].children) == (g[i].type==Sample) for i ∈ tips(g))
    @test timezero(g) == 0.0
    @test g.time == 20.0

    @info h2("newick round trip (necessary, not sufficient)")
    s = newick(g)
    @test s isa Vector{String}
    @test length(s) == 1
    g2 = parse_newick(s, time = g.time)
    @test g2 isa Genealogy{PhyloPOMP.Unstructured}
    @test length(g2) == length(g)
    @test nsample(g2) == nsample(g)
    @test times(g2) ≈ times(g) atol=1e-4

    @info h2("cblv round trip (necessary, not sufficient)")
    xy = cblv(g)
    g3 = parse_cblv(xy...; time = g.time)
    @test nsample(g3) == nsample(g)
    xy2 = cblv(g3)
    @test xy2[1] ≈ xy[1] atol=1e-4
    @test xy2[2] ≈ xy[2] atol=1e-4

    @info h2("simulate -> filter: finite log-likelihood (headline check)")
    p = NaiveSEIR.filter_pomp(
        g;
        β=β, σ=σ, γ=γ, ω=ω, ψ=ψ, χ=χ, pop=pop,
        S0=x0.S/pop, E0=x0.E/pop, I0=x0.I/pop, R0=x0.R/pop,
    )
    @test p isa POMP.PompObject
    pf = pfilter(p, Np=1000)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(logLik(pf))

    @testset "simulate -> GuidedSEIR: check passes and logLik finite" begin
        @info h2("simulate -> GuidedSEIR: check passes and logLik finite")
        ## The SEIR `sampling` event (`src/examples/mgp_macro.jl:197`) is
        ## NON-destructive -- `pop=()`, so `simulate` leaves the sampled
        ## lineage open (src/simulate.jl:197-198) and a later transmission
        ## can hang a child off the Sample node itself. Such "inline"
        ## samples are exactly what `GuidedSEIR.check` permits (`< 2`
        ## children, seir_guided.jl:46) and what
        ## `GuidedSEIR.inline_sample!` (seir_guided.jl:140) handles, so the
        ## loop below insists on a realization that contains at least one:
        ## it is the case a destructive-sampling simulator would get wrong.
        ## Independently seeded from the headline `g` above, so this is a
        ## second realization rather than a re-run on one lucky tree.
        rng2 = MersenneTwister(20260901)
        local gs
        ok2 = false
        inline2 = false
        for _ ∈ 1:300
            gs = simulate(PhyloPOMP.SEIR, θ; x0=x0, graft=[0,1], tmax=20.0, rng=rng2)
            inline2 = any(length(gs[i].children)==1 for i ∈ samples(gs))
            ## cap the sample count to keep the pfilter below cheap
            if 5 ≤ nsample(gs) ≤ 40 && inline2
                ok2 = true
                break
            end
        end
        @test ok2
        @test inline2
        ## `check` returns `nothing` and throws on failure, so calling it at
        ## all is the assertion; `=== nothing` records it as a passing test.
        @test GuidedSEIR.check(gs) === nothing
        ## ...and the invariants it asserts, spelled out, so a regression
        ## names the broken one instead of just firing an @assert:
        @test all(length(gs[i].children)==1 for i ∈ roots(gs))
        @test all(length(gs[i].children)<2 for i ∈ samples(gs))
        @test all(length(gs[i].children)==2 for i ∈ nodes(gs))

        ## same numbers as the simulation: `seir_rinit`
        ## (src/examples/seir_guided.jl:298) rescales by pop/(S0+E0+I0+R0),
        ## so fractions of `pop` reproduce x0 exactly. NOTE θ.N ↔ `pop`, and
        ## χ must be 0: the `@mgp SEIR` table has no χ-event at all, its only
        ## sampling event being `rate=ψ*I` (mgp_macro.jl:197).
        pg = GuidedSEIR.filter_pomp(
            gs, fsmarkov(Expos=>0.1, Infec=>1, (Expos,Infec)=>1);
            β=β, σ=σ, γ=γ, ω=ω, ψ=ψ, χ=χ, pop=pop,
            S0=x0.S/pop, E0=x0.E/pop, I0=x0.I/pop, R0=x0.R/pop,
        )
        @test pg isa POMP.PompObject
        ## pfilter draws from the global RNG, so seed for reproducibility
        ## independent of whatever ran earlier in the suite.
        seed!(20260815)
        pfg = pfilter(pg, Np=1000)
        @test pfg isa POMP.PfilterdPompObject
        @test isfinite(logLik(pfg))
    end

end

end
