module GuidedSEIRTest

import ..Main: h1, h2

@info h1("SEIR model with guided proposals")

using PhyloPOMP
import PartiallyObservedMarkovProcesses as POMP
using Test
using BenchmarkTools
using Statistics: std

using PhyloPOMP.GuidedSEIR: seir, seir_convert
using PhyloPOMP.GuidedSEIR.SEIR: Expos, Infec, T as DemeType

@testset verbose=true "SEIR model with guided proposals" begin

    g1 = parse_newick(readlines("seir1.nwk"), time = 50.0)
    @test g1 isa Genealogy{PhyloPOMP.Unstructured.T}

    g2 = geneal_convert(g1,seir_convert,DemeType)
    @test g2 isa Genealogy{DemeType}

    g3 = guide(g2,fsmarkov(Expos=>0.1,Infec=>1,(Expos,Infec)=>1))

    p = seir(g3)
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
