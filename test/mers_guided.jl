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
        cooling_schedule=geometric_cooling(1.0),
        perturbations=function(
            s, lag;
            βcc, βhh, βch, βhc,
            Sc0, Ic0, Sh0, Ih0,
            _...,
            )
            βcc = rand(LogNormal(log(βcc),0.002*s))
            βhh = rand(LogNormal(log(βhh),0.002*s))
            βhc = rand(LogNormal(log(βhc),0.002*s))
            βch = rand(LogNormal(log(βch),0.002*s))
            if lag == 0
                Sc0 = rand(LogNormal(log(Sc0),0.2*s))
                Ic0 = rand(LogNormal(log(Ic0),0.2*s))
                Sh0 = rand(LogNormal(log(Sh0),0.2*s))
                Ih0 = rand(LogNormal(log(Ih0),0.2*s))
                m = Sc0 + Ic0
                Sc0 /= m
                Ic0 /= m
                m = Sh0 + Ih0
                Sh0 /= m
                Ih0 /= m
            end
            (;βcc,βhh,βch,βhc,Sc0,Ic0,Sh0,Ih0)
        end,
    )
    @test mf isa POMP.MifdPompObject

    if heavy
        @time mf = mif(mf,Nmif=100)
        @time mf1 = mif(mf,Nmif=50,cooling_schedule=geometric_cooling(0.8))
        @time mf2 = mif(mf1,Nmif=100,cooling_schedule=geometric_cooling(0.1))
        @info "$(coef(mf2))"
    end

end

end
