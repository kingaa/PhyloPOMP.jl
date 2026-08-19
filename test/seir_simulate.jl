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
"""
module SEIRSimulateTest

import ..Main: h1, h2

@info h1("SEIR forward simulator")

using Test
using Random: MersenneTwister
using PhyloPOMP
using PhyloPOMP: Sample
using PhyloPOMP.NaiveSEIR
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

end

end
