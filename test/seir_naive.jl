module NaiveSEIRTest

import ..Main: h1, h2

@info h1("SEIR model with naïve proposals")

using Test
using BenchmarkTools
using Random: seed!
using PhyloPOMP
using PhyloPOMP.NaiveSEIR
import PartiallyObservedMarkovProcesses as POMP

@testset verbose=true "SEIR model with naïve proposals" begin

    seed!(2121916527)

    g = parse_newick(seir_trees[1], time = 50.0)
    @test g isa Genealogy{PhyloPOMP.Unstructured}

    p = seir(g,E0=0,I0=0)
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=100))==-Inf

    p = seir(g,χ=0.01)
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
