module GuidedSEIRTest

import ..Main: h1, h2

@info h1("SEIR model with guided proposals")

using Test
using BenchmarkTools
using Random: seed!
using PhyloPOMP
using PhyloPOMP.GuidedSEIR
using PhyloPOMP.GuidedSEIR.Demes: Expos, Infec
import PartiallyObservedMarkovProcesses as POMP

@testset verbose=true "SEIR model with guided proposals" begin

    seed!(2121916527)

    g1 = parse_newick(seir_trees[1], time = 50.0)
    @test g1 isa Genealogy{PhyloPOMP.Unstructured.DemeSet}

    g2 = guide(
        g1,
        fsmarkov(Expos=>0.1,Infec=>1,(Expos,Infec)=>1),
        seir_convert!
    )

    p = seir(g2,E0=0,I0=0)
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=100))==-Inf

    p = seir(g2)
    @test p isa POMP.PompObject

    @info h2("simulate test")
    sm = simulate(p, nsim = 3)
    @time sm = simulate(p, nsim = 3)
    @test sm isa Matrix{<:POMP.PompObject}

    @info h2("pfilter test")
    pf = pfilter(p, Np = 100)
    @time pf = pfilter(p, Np = 100)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(logLik(pf))

    @info h2("pfilter benchmark")
    @btime pfilter($p, Np = 1000)

    ll = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
    llest,llse = logmeanexp(ll,se=true)
    @info "logLik = $(round(llest,digits=2)) ± $(round(llse,sigdigits=3))"

end

end
