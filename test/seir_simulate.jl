"""
`simulate` on `@mgp SEIR`, then `NaiveSEIR.filter_pomp` and `GuidedSEIR.filter_pomp` on the result.
Both should return a finite log likelihood.
`GuidedSEIR.check` also requires root degree 1, internal degree 2, and sample degree below 2.
"""
module SEIRSimulateTest

## `h1`/`h2` come from runtests.jl; fall back when run standalone.
const h1 = isdefined(Main, :h1) ? Main.h1 : identity
const h2 = isdefined(Main, :h2) ? Main.h2 : identity
isdefined(Main, :SimulateChecks) || Base.include(Main, joinpath(@__DIR__, "simulate_checks.jl"))

@info h1("SEIR forward simulator")

using Test
using Random: MersenneTwister, seed!
using PhyloPOMP
using PhyloPOMP: Sample
using Main.SimulateChecks
using PhyloPOMP.NaiveSEIR
using PhyloPOMP.GuidedSEIR
using PhyloPOMP.GuidedSEIR.Demes: Expos, Infec
import PartiallyObservedMarkovProcesses as POMP

@testset verbose=true "SEIR forward simulator" begin

    ## Same numbers for the simulator and the filter.
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

    @info h2("per-event timing, structure, round trips (single root and forest)")
    ## These look at node times. The filter tests below only ask for a finite likelihood.
    cov = run_invariants(PhyloPOMP.SEIR, θ, [(x0, [0,1]), ((S=95, E=2, I=3, R=0), [2,3])];
                         ntree=100, tmax=20.0, sample_max_children=1)
    @test cov.nonempty > 100
    @test cov.inline > 0        # ψ-sampling is non-destructive: inline samples must occur

    @info h2("χ > 0: destructive `culling` samples terminate their lineage")
    ## χ > 0: ψ-samples may have a child, χ-culls are leaves.
    ## NaiveSEIR draws ψ against χ at each sample.
    let θχ = merge(θ, (ψ = 0.05, χ = 0.10))
        covχ = run_invariants(PhyloPOMP.SEIR, θχ, [(x0, [0,1])];
                              ntree=100, seed=20261004, tmax=20.0, sample_max_children=1)
        @test covχ.nonempty > 50
        @test covχ.inline > 0
        rngχ = MersenneTwister(20261004)
        leaf_samples = 0; gχ = g
        for _ ∈ 1:50
            gχ = simulate(PhyloPOMP.SEIR, θχ; x0=x0, graft=[0,1], tmax=20.0, rng=rngχ)
            leaf_samples += count(isempty(gχ[i].children) for i ∈ samples(gχ))
            5 ≤ nsample(gχ) ≤ 40 && break
        end
        @test leaf_samples > 0
        pχ = NaiveSEIR.filter_pomp(gχ; β=β, σ=σ, γ=γ, ω=ω, ψ=θχ.ψ, χ=θχ.χ, pop=pop,
                                   S0=x0.S/pop, E0=x0.E/pop, I0=x0.I/pop, R0=x0.R/pop)
        seed!(20261004)
        @test isfinite(logLik(pfilter(pχ, Np=500)))
        ## With χ = 0 the culling hazard is 0.
        @test all(ev -> ev.hazard((S=90,E=4,I=5,R=1), θ) ≥ 0, PhyloPOMP.SEIR.events)
        @test PhyloPOMP.SEIR.events[end].hazard((S=90,E=4,I=5,R=1), θ) == 0
    end

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
        ## SEIR sampling leaves the host in place, so a later birth can put a
        ## child on the Sample node. `GuidedSEIR.check` allows that.
        ## Separate seed from `g` above.
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
        @test GuidedSEIR.check(gs) === nothing
        @test all(length(gs[i].children)==1 for i ∈ roots(gs))
        @test all(length(gs[i].children)<2 for i ∈ samples(gs))
        @test all(length(gs[i].children)==2 for i ∈ nodes(gs))

        ## `seir_rinit` rescales by `pop / (S0+E0+I0+R0)`, so fractions of `pop` give `x0` back. `θ.N` is `pop`.
        pg = GuidedSEIR.filter_pomp(
            gs, fsmarkov(Expos=>0.1, Infec=>1, (Expos,Infec)=>1);
            β=β, σ=σ, γ=γ, ω=ω, ψ=ψ, χ=χ, pop=pop,
            S0=x0.S/pop, E0=x0.E/pop, I0=x0.I/pop, R0=x0.R/pop,
        )
        @test pg isa POMP.PompObject
        ## pfilter uses the global RNG, so seed it.
        seed!(20260815)
        pfg = pfilter(pg, Np=1000)
        @test pfg isa POMP.PfilterdPompObject
        @test isfinite(logLik(pfg))
    end

end

end
