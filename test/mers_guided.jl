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

heavy = occursin(r"y|yes|t|true", get(ENV,"RUN_HEAVY_TESTS","yes"))

@testset verbose=true "MERS model with guided proposals" begin

    seed!(2121916527)

    m0 = fsmarkov(Camel=>0.5, Human=>0.5, (Camel,Human)=>0.01)
    g = parse_newick(GuidedMERS.mers_newick, demes=GuidedMERS.Demes)
    p = GuidedMERS.filter_pomp(g, m0)
    @time pf = pfilter(p,Np=1000,trigger=0.2,target=0.8)
    @time mf = mif(
        pf,Nmif=100,Np=1000,trigger=0.2,target=0.8,
        cooling_schedule=geometric_cooling(1.0),
        perturbations=function(
            s, lag;
            β_cc, β_hh, β_ch, β_hc,
            S_c0, I_c0, S_h0, I_h0,
            _...,
            )
            β_cc = rand(LogNormal(log(β_cc),0.002*s))
            β_hh = rand(LogNormal(log(β_hh),0.002*s))
            β_hc = rand(LogNormal(log(β_hc),0.002*s))
            β_ch = rand(LogNormal(log(β_ch),0.002*s))
            if lag == 0
                S_c0 = rand(LogNormal(log(S_c0),0.2*s))
                I_c0 = rand(LogNormal(log(I_c0),0.2*s))
                S_h0 = rand(LogNormal(log(S_h0),0.2*s))
                I_h0 = rand(LogNormal(log(I_h0),0.2*s))
                m = S_c0 + I_c0
                S_c0 /= m
                I_c0 /= m
                m = S_h0 + I_h0
                S_h0 /= m
                I_h0 /= m
            end
            (;β_cc,β_hh,β_ch,β_hc,S_c0,I_c0,S_h0,I_h0)
        end,
    )
    @time mf = mif(mf,Nmif=20)
    @time mf = mif(mf,Nmif=50,cooling_schedule=geometric_cooling(0.8))

end
