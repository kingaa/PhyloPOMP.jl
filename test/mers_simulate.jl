"""
Milestone-2 acceptance test for the generic forward simulator
(`src/simulate.jl`), generalizing the SEIR check in `test/seir_simulate.jl`
to a second, independent model: MERS's two-host (Camel/Human) structure,
already declared via `@mgp MERS` in `src/examples/mgp_mers.jl` -- no new
model-spec front-end code was needed for this milestone, only the
`demeset`/`samplemap` generalization to `simulate` itself (MERS observes
deme identity at the moment of sampling, unlike SEIR, where deme is purely
latent -- see `src/simulate.jl`'s module header).

`NaiveMERS.filter_pomp` cannot be used as the validating filter here: unlike
`NaiveSEIR.filter_pomp`, it takes no genealogy argument at all and always
filters the fixed empirical `mers_tree` (confirmed in `handoff.md`'s
"Follow-up round 4"). `SoftMERS.filter_pomp(gen, m; ...)` (shared machinery
in `src/examples/mers_funs.jl`) DOES accept an arbitrary genealogy, exactly
like the SEIR case, so it is the validating filter used below.
"""
module MERSSimulateTest

import ..Main: h1, h2

@info h1("MERS forward simulator")

using Test
using Random: MersenneTwister, seed!
using PhyloPOMP
using PhyloPOMP: Sample
using PhyloPOMP.SoftMERS
using PhyloPOMP.SoftMERS.Demes: Camel, Human
import PartiallyObservedMarkovProcesses as POMP

@testset verbose=true "MERS forward simulator" begin

    ## parameters shared between the simulator and the filter. Kept modest
    ## (small N, R0 not too far above 1) so a "cross-species, few-sample"
    ## realization is common within the retry budget below, rather than
    ## needing to fish for a rare narrow window between "died out" and "large
    ## outbreak" -- see the diagnostic in the implementation notes.
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
        ## demand at least one sample of EACH species (so the test actually
        ## exercises the cross-deme transmission_hc/transmission_ch events,
        ## not just the single-host case), and cap the sample count so the
        ## tree stays cheap and non-collapse-prone for pfilter below --
        ## an uncapped major outbreak can reach 50+ samples, which needs
        ## far more than a modest Np to reliably avoid -Inf.
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
    ## MERS sampling is always destructive (`sample_remove`): every Sample
    ## node has exactly zero children, matching NaiveMERS/SoftMERS's
    ## `@assert length(n.children)==0` invariant.
    @test all(isempty(g[i].children) for i ∈ samples(g))
    @test all(!ismissing(g[i].deme) for i ∈ samples(g))
    @test all(ismissing(g[i].deme) for i ∈ eachindex(g) if g[i].type != Sample)

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
    ## pfilter draws from the global RNG (unlike `simulate` above, which was
    ## given its own `rng`), so seed explicitly for a reproducible result
    ## independent of what ran earlier in the suite.
    seed!(20260814)
    pf = pfilter(p, Np=5000)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(logLik(pf))

end

end
