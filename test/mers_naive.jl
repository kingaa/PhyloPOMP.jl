module NaiveMERSTest

import ..Main: h1, h2

@info h1("MERS model with naïve proposals")

using Test
using BenchmarkTools
using Random: seed!
using PhyloPOMP
using PhyloPOMP.NaiveMERS
import PartiallyObservedMarkovProcesses as POMP

@testset verbose=true "MERS model with naïve proposals" begin

    seed!(2121916527)

    p = NaiveMERS.filter_pomp(I_c0=0,I_h0=0)
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=100))==-Inf

    p = NaiveMERS.filter_pomp()
    @test p isa POMP.PompObject

    @info h2("simulate test")
    sm = simulate(p, nsim = 3)
    @time sm = simulate(p, nsim = 3)
    @test sm isa Matrix{<:POMP.PompObject}

    @info h2("pfilter test")
    pf = pfilter(p, Np = 100)
    @time pf = pfilter(p, Np = 100)
    @test pf isa POMP.PfilterdPompObject

    @info h2("pfilter benchmark")
    @btime pfilter($p, Np = 1000)

    @time ll = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
    llest,llse = logmeanexp(ll,se=true)
    @info "logLik = $(round(llest,digits=2)) ± $(round(llse,sigdigits=3))"

end

end
