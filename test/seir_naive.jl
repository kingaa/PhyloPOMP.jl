module NaiveSEIRTest

import ..Main: h1, h2

@info h1("SEIR model with naïve proposals")

using Test
using BenchmarkTools
using Statistics: std
using Random: seed!
using PhyloPOMP
import PartiallyObservedMarkovProcesses as POMP
using PhyloPOMP.NaiveSEIR: seir

@testset verbose=true "SEIR model with naïve proposals" begin

    seed!(2121916527)

    g = parse_newick(readlines("seir1.nwk"), time = 50.0)
    @test g isa Genealogy{PhyloPOMP.Unstructured.T}

    p = seir(g)
    @test p isa POMP.PompObject

    @info h2("simulate test")
    sm = simulate(p, nsim = 3)
    @time sm = simulate(p, nsim = 3)
    @test sm isa Matrix{<:POMP.PompObject}

    @info h2("pfilter test")
    pf = pfilter(p, Np = 100)
    @time pf = pfilter(p, Np = 100)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(pf.logLik)

    @info h2("pfilter benchmark")
    @btime pfilter($p, Np = 1000)

    logmeanexp(x) = begin
        xmax = maximum(x)
        xmax + log(sum(exp.(x .- xmax))) - log(length(x))
    end

    x = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
    @info "logLik = $(round(logmeanexp(x),digits=2)) ± $(round(std(x),sigdigits=3))"

end

end
