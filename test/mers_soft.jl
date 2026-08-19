module SoftMERSTest

import ..Main: h1, h2

@info h1("MERS model with soft proposals")

using Test
using BenchmarkTools
using Random: seed!
using PhyloPOMP
using PhyloPOMP.SoftMERS
using PhyloPOMP.SoftMERS.Demes: Camel, Human
import PartiallyObservedMarkovProcesses as POMP

heavy = occursin(r"y|yes|t|true", get(ENV,"RUN_HEAVY_TESTS","yes"))

## A small binary tree with both camel and human tips, requiring at least
## one cross-deme transmission to explain -- small enough for a modest
## particle count to reliably find compatible colorings (the full empirical
## `mers_tree`, used only for the -Inf/constructor checks below, is far too
## large for that).
const small_tree = "(([&&PhyloPOMP deme=camel]1:1.0,[&&PhyloPOMP deme=camel]2:1.0):0.5,[&&PhyloPOMP deme=human]3:1.5);"

@testset verbose=true "MERS model with soft proposals" begin

    seed!(2121916527)

    m0 = fsmarkov(Camel=>0.3, Human=>0.7, (Camel,Human)=>1)
    p = SoftMERS.filter_pomp(SoftMERS.mers_tree, m0, I_c0=0, I_h0=0)
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=50))==-Inf

    g = parse_newick(small_tree, demes=SoftMERS.Demes)
    m = fsmarkov(Camel=>0.3, Human=>0.7, (Camel,Human)=>1)
    common = (
        β_hc=1.0, β_ch=1.0, χ_c=1.0, χ_h=1.0,
        S_c0=1.0, S_h0=1.0, I_c0=0.5, I_h0=0.5, N_c=100.0, N_h=100.0,
    )
    p = SoftMERS.filter_pomp(g, m; common...)
    @test p isa POMP.PompObject

    @info h2("simulate test")
    sm = simulate(p, nsim = 3)
    @time sm = simulate(p, nsim = 3)
    @test sm isa Matrix{<:POMP.PompObject}

    @info h2("pfilter test")
    pf = pfilter(p, Np = 5000)
    @time pf = pfilter(p, Np = 5000)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(logLik(pf))

    @info h2("soft vs guided: mathematically distinct kernels")
    gG = parse_newick(small_tree, demes=PhyloPOMP.GuidedMERS.Demes)
    CamelG, HumanG = PhyloPOMP.GuidedMERS.Camel, PhyloPOMP.GuidedMERS.Human
    mG = fsmarkov(CamelG=>0.01, HumanG=>0.99, (CamelG,HumanG)=>0.01)
    mS = fsmarkov(Camel=>0.01, Human=>0.99, (Camel,Human)=>0.01)
    seed!(7)
    pS2 = SoftMERS.filter_pomp(g, mS; common...)
    llS = [logLik(pfilter(pS2,Np=300)) for _ ∈ 1:5]
    seed!(7)
    pG2 = PhyloPOMP.GuidedMERS.filter_pomp(gG, mG; common...)
    llG = [logLik(pfilter(pG2,Np=300)) for _ ∈ 1:5]
    @test !all(isapprox.(llS, llG; atol=0.5))

    if heavy
        @info h2("pfilter benchmark")
        @btime pfilter($p, Np = 1000)

        @time ll = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
        llest,llse = logmeanexp(ll,se=true)
        @info "logLik = $(round(llest,digits=2)) ± $(round(llse,sigdigits=3))"
    end

end

end
