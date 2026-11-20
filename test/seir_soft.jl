module SoftSEIRTest

import ..Main: h1, h2

@info h1("SEIR model with soft proposals")

using Test
using BenchmarkTools
using Random: seed!
using PhyloPOMP
using PhyloPOMP: SoftSEIR
import PartiallyObservedMarkovProcesses as POMP

heavy = occursin(r"y|yes|t|true", get(ENV,"RUN_HEAVY_TESTS","no"))

@testset verbose=true "SEIR model with soft proposals" begin

    seed!(2123986527)

    g = @guide SoftSEIR SoftSEIR.geneal[1] fsmarkov(Expos=>0.1, Infec=>1, (Expos,Infec)=>1)
    @test g isa Guide

    p = filter_pomp(
        SoftSEIR, g,
        params=(
            β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
            pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.0, R0 = 0.08,
        )
    )
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=100))==-Inf

    p = filter_pomp(
        SoftSEIR, g,
        params=(
            β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.0, χ = 0.0,
            pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
        )
    )
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=100))==-Inf

    p = filter_pomp(
        SoftSEIR, g,
        params=(
            β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.01,
            pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
        )
    )
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

    if heavy
        @time ll = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
        llest,llse = logmeanexp(ll,se=true)
        @info "logLik = $(round(llest,digits=2)) ± $(round(llse,sigdigits=3))"
    end

end

end
