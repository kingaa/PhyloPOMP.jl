"""
`simulate` on `@mgp MERS`, then `SoftMERS.filter_pomp` and `GuidedMERS.filter_pomp`.
MERS records the deme of a sample, so this file is where `demeset` and `samplemap` get used.
`GuidedMERS` has to be given `demeset = GuidedMERS.Demes`. A `SoftMERS` deme is a different type, and `GuidedMERS` reads `gen[node].deme` directly.
"""
module MERSSimulateTest

## `h1`/`h2` come from runtests.jl; fall back when run standalone.
const h1 = isdefined(Main, :h1) ? Main.h1 : identity
const h2 = isdefined(Main, :h2) ? Main.h2 : identity
isdefined(Main, :SimulateChecks) || Base.include(Main, joinpath(@__DIR__, "simulate_checks.jl"))

@info h1("MERS forward simulator")

using Test
using Random: MersenneTwister, seed!
using PhyloPOMP
using PhyloPOMP: Sample
using Main.SimulateChecks
using PhyloPOMP.SoftMERS
using PhyloPOMP.SoftMERS.Demes: Camel, Human
using PhyloPOMP.GuidedMERS
import PartiallyObservedMarkovProcesses as POMP

## NB: `Camel`/`Human` above are SoftMERS's; GuidedMERS's same-named
## instances belong to a DIFFERENT enum type and are always qualified
## (`GuidedMERS.Camel`, `GuidedMERS.Human`) below.

@testset verbose=true "MERS forward simulator" begin

    ## Parameters shared by simulator and filter. Small N and R0 near 1 make
    ## cross-species, few-sample trees common within the retry budget.
    β_cc, β_ch, β_hc, β_hh = 3.0, 0.5, 0.5, 3.0
    γ_c, γ_h = 1.0, 1.0
    χ_c, χ_h = 0.3, 0.3
    B_c, B_h = 0.0, 0.0
    N_c, N_h = 20, 20
    θ = (
        β_cc=β_cc, β_ch=β_ch, β_hc=β_hc, β_hh=β_hh,
        γ_c=γ_c, γ_h=γ_h, χ_c=χ_c, χ_h=χ_h,
        B_c=B_c, B_h=B_h, N_c=Float64(N_c), N_h=Float64(N_h),
    )
    x0 = (S_c=N_c-1, I_c=1, S_h=N_h, I_h=0)
    samplemap = [Camel, Human]

    ## The same simulation, recording sample demes as GuidedMERS's enum
    ## instances, with the retry loop shared by the two guided blocks below.
    ## Returns `(genealogy, ok)`; `lo` is the minimum sample count.
    guided_sim(rng, lo) = begin
        local gsim
        okl = false
        for _ ∈ 1:300
            gsim = simulate(
                PhyloPOMP.MERS, θ; x0=x0, graft=[1,0], tmax=10.0, rng=rng,
                demeset=GuidedMERS.Demes,
                samplemap=[GuidedMERS.Camel, GuidedMERS.Human],
            )
            demes_sampled = Set(gsim[i].deme for i ∈ samples(gsim))
            if lo ≤ nsample(gsim) ≤ 12 &&
                GuidedMERS.Camel ∈ demes_sampled &&
                GuidedMERS.Human ∈ demes_sampled
                okl = true
                break
            end
        end
        gsim, okl
    end

    ## Filter parameters match the simulation; `mers_rinit` rescales each
    ## species by N/(S0+I0).
    guided_pomp(gen, m) = GuidedMERS.filter_pomp(
        gen, m;
        β_cc=β_cc, β_ch=β_ch, β_hc=β_hc, β_hh=β_hh,
        γ_c=γ_c, γ_h=γ_h, χ_c=χ_c, χ_h=χ_h, B_c=B_c, B_h=B_h,
        S_c0=(N_c-1)/N_c, I_c0=1/N_c, S_h0=1.0, I_h0=0.0,
        N_c=N_c, N_h=N_h,
    )

    @info h2("argument validation")
    @test_throws ArgumentError simulate(
        PhyloPOMP.MERS, θ; x0=x0, graft=[1,0,0], tmax=10.0,
        demeset=SoftMERS.Demes, samplemap=samplemap,
    )
    @test_throws ArgumentError simulate(
        PhyloPOMP.MERS, θ; x0=x0, graft=[0,0], tmax=10.0,
        demeset=SoftMERS.Demes, samplemap=samplemap,
    )
    @test_throws ArgumentError simulate(
        PhyloPOMP.MERS, θ; x0=x0, graft=[1,0], tmax=10.0,
        demeset=SoftMERS.Demes, samplemap=[Camel],
    )

    @info h2("simulate a non-degenerate, cross-species genealogy")
    rng = MersenneTwister(20260813)
    local g
    ok = false
    for _ ∈ 1:300
        g = simulate(
            PhyloPOMP.MERS, θ; x0=x0, graft=[1,0], tmax=10.0, rng=rng,
            demeset=SoftMERS.Demes, samplemap=samplemap,
        )
        ## Require both species sampled (exercises transmission_hc/ch) and at
        ## most 12 samples (keeps pfilter cheap).
        demes_sampled = Set(g[i].deme for i ∈ samples(g))
        if 4 ≤ nsample(g) ≤ 12 && Camel ∈ demes_sampled && Human ∈ demes_sampled
            ok = true
            break
        end
    end
    @test ok
    @test g isa Genealogy{SoftMERS.Demes}
    @test length(roots(g)) == 1
    @test nsample(g) == length(samples(g))
    ## MERS sampling is destructive: Sample nodes have no children.
    @test all(isempty(g[i].children) for i ∈ samples(g))
    @test all(!ismissing(g[i].deme) for i ∈ samples(g))
    @test all(ismissing(g[i].deme) for i ∈ eachindex(g) if g[i].type != Sample)

    @info h2("per-event timing, structure, round trips (single root and forest, demography on)")
    ## Same checks as SEIR. Sampling removes the host, and births are on (`B_c`, `B_h` > 0).
    θd = merge(θ, (B_c=0.5, B_h=0.5))
    cov = run_invariants(PhyloPOMP.MERS, θd,
                         [(x0, [1,0]), ((S_c=18, I_c=2, S_h=19, I_h=1), [2,1])];
                         ntree=100, tmax=10.0, demeset=SoftMERS.Demes, samplemap=samplemap,
                         sample_max_children=0)
    @test cov.nonempty > 100
    @test cov.inline == 0       # sample_remove: a Sample node never has a child

    @info h2("newick round trip (necessary, not sufficient)")
    s = newick(g)
    @test length(s) == 1
    g2 = parse_newick(s, demes=SoftMERS.Demes, time = g.time)
    @test g2 isa Genealogy{SoftMERS.Demes}
    @test length(g2) == length(g)
    @test nsample(g2) == nsample(g)
    @test times(g2) ≈ times(g) atol=1e-4
    @test Set(g2[i].deme for i ∈ samples(g2)) == Set(g[i].deme for i ∈ samples(g))

    @info h2("cblv round trip (necessary, not sufficient)")
    xy = cblv(g)
    g3 = parse_cblv(xy...; demes=SoftMERS.Demes, time = g.time)
    @test nsample(g3) == nsample(g)

    @info h2("simulate -> filter: finite log-likelihood (headline check)")
    m = fsmarkov(Camel=>0.3, Human=>0.7, (Camel,Human)=>1)
    p = SoftMERS.filter_pomp(
        g, m;
        β_cc=β_cc, β_ch=β_ch, β_hc=β_hc, β_hh=β_hh,
        γ_c=γ_c, γ_h=γ_h, χ_c=χ_c, χ_h=χ_h, B_c=B_c, B_h=B_h,
        S_c0=(N_c-1)/N_c, I_c0=1/N_c, S_h0=1.0, I_h0=0.0,
        N_c=N_c, N_h=N_h,
    )
    @test p isa POMP.PompObject
    ## pfilter uses the global RNG, so seed it.
    seed!(20260814)
    pf = pfilter(p, Np=5000)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(logLik(pf))

    ## demeset/samplemap only label nodes and use no randomness, so the same
    ## seed replays g; the newick equality checks that.
    gg, okg = guided_sim(MersenneTwister(20260813), 4)
    @test okg
    @test gg isa Genealogy{GuidedMERS.Demes}
    @test newick(gg) == newick(g)
    mg = fsmarkov(
        GuidedMERS.Camel=>0.5, GuidedMERS.Human=>0.5,
        (GuidedMERS.Camel,GuidedMERS.Human)=>0.05,
    )
    pg = guided_pomp(gg, mg)
    @test pg isa POMP.PompObject
    seed!(20260814)
    pfg = pfilter(pg, Np=2000)
    @test pfg isa POMP.PfilterdPompObject
    @test isfinite(logLik(pfg))

    @testset "simulate -> GuidedMERS: check passes and logLik finite" begin
        ## Independent realization. GuidedMERS.check returns nothing or throws.
        @info h2("simulate -> GuidedMERS: check passes and logLik finite")
        g4, ok4 = guided_sim(MersenneTwister(20260901), 3)
        @test ok4
        @test 3 ≤ nsample(g4) ≤ 12
        @test GuidedMERS.check(g4) === nothing
        ## The properties check asserts, so a failure names the broken one.
        @test all(length(g4[i].children)==1 for i ∈ roots(g4))
        @test all(isempty(g4[i].children) for i ∈ samples(g4))
        @test all(
            g4[i].deme ∈ [GuidedMERS.Camel, GuidedMERS.Human]
            for i ∈ samples(g4)
        )
        @test all(length(g4[i].children)==2 for i ∈ nodes(g4))

        m4 = fsmarkov(
            GuidedMERS.Camel=>0.5, GuidedMERS.Human=>0.5,
            (GuidedMERS.Camel,GuidedMERS.Human)=>0.05,
        )
        p4 = guided_pomp(g4, m4)
        @test p4 isa POMP.PompObject
        seed!(20260815)
        pf4 = pfilter(p4, Np=2000)
        @test pf4 isa POMP.PfilterdPompObject
        @test isfinite(logLik(pf4))
    end

end

end
