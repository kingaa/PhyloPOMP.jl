module GuidedMERSTest

import ..Main: h1, h2

@info h1("MERS model with guided proposals")

using Test
using BenchmarkTools
using Random: seed!
using Distributions: LogNormal
using PhyloPOMP
using PhyloPOMP.GuidedMERS
using PhyloPOMP.GuidedMERS.Demes: Camel, Human
import PartiallyObservedMarkovProcesses as POMP

heavy = occursin(r"y|yes|t|true", get(ENV,"RUN_HEAVY_TESTS","no"))

@testset verbose=true "MERS model with guided proposals" begin

    seed!(2121916527)

    m0 = fsmarkov(Camel=>0.5, Human=>0.5, (Camel,Human)=>0.01)
    g = parse_newick(GuidedMERS.mers_newick, demes=GuidedMERS.Demes)
    p = GuidedMERS.filter_pomp(g, m0)
    @test p isa POMP.PompObject
    pf = pfilter(p,Np=10,trigger=0.2,target=0.8)
    @time pf = pfilter(p,Np=1000,trigger=0.2,target=0.8)
    @test pf isa POMP.PfilterdPompObject
    mf = mif(
        pf,Nmif=3,Np=1000,trigger=0.2,target=0.8,
        cooling=geometric_cooling(1.0),
        perturbations=@perturbn(
            @lognormal(βcc,0.002),
            @lognormal(βhh,0.002),
            @lognormal(βhc,0.002),
            @lognormal(βch,0.002),
            @ivp(@logbarynormal((Sc0,Ic0),0.2)),
            @ivp(@logbarynormal((Sh0,Ih0),0.2)),
        )
    )
    @test mf isa POMP.MifdPompObject

    if heavy
        @time mf = mif(mf,Nmif=100)
        @time mf1 = mif(mf,Nmif=50,cooling=geometric_cooling(0.8))
        @time mf2 = mif(mf1,Nmif=100,cooling=geometric_cooling(0.1))
        @info "estimates: $(map(x->round(x,sigdigits=3),coef(mf2)))"
        @info "logLik estimate = $(logLik(mf2))"
    end

end

end
