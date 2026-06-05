module GuidedSEIRTest

import ..Main: h1, h2

@info h1("SEIR model with guided proposals")

using Test
using BenchmarkTools
using Statistics: std
using Random: seed!
using PhyloPOMP
import PartiallyObservedMarkovProcesses as POMP
using PhyloPOMP: Sample, Node, Root
using PhyloPOMP.GuidedSEIR: seir
using PhyloPOMP.GuidedSEIR.SEIR: Expos, Infec

@testset verbose=true "SEIR model with guided proposals" begin

    seed!(2121916527)

    g1 = parse_newick(readlines("seir1.nwk"), time = 50.0)
    @test g1 isa Genealogy{PhyloPOMP.Unstructured.DemeSet}

    ## This function should return true if the guide probabilities
    ## will be fixed at this node and false otherwise. If the former,
    ## it should fill the vector `v` with an appropriate probability
    ## vector.
    seir_convert!(
        v; deme, type, time,
    ) = begin
        if type==Sample || type==Node
            demekron!(v,Infec)
            true
        else
            false
        end
    end

    g2 = guide(
        g1,
        fsmarkov(Expos=>0.1,Infec=>1,(Expos,Infec)=>1),
        seir_convert!
    )

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
