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

## A small tree with camel and human tips that needs at least one cross-deme
## transmission. The full `mers_tree` is too large for a modest particle count,
## so it is used only for the -Inf and constructor checks.
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

    @info h2("default parameters give a finite logLik")
    ## Shared defaults (mers_funs.jl). N_c = N_h = 10000 makes a 3-tip tree hard
    ## to filter, so assert only that the estimate over 5 replicates is finite.
    seed!(1)
    pdef = SoftMERS.filter_pomp(g, m)
    lldef,_ = logmeanexp([logLik(pfilter(pdef,Np=2000)) for _ ∈ 1:5], se=true)
    @info "defaults: logLik = $(round(lldef,digits=2))"
    @test isfinite(lldef)

    @info h2("soft, hard and guided agree on the likelihood")
    ## Soft, hard and guided are proposals for one likelihood, so their pfilter
    ## estimates agree within Monte Carlo error (z < 3). 15 replicates at
    ## Np=2000 keep the heavy-tailed estimator stable.
    kernels = (
        ("soft",   PhyloPOMP.SoftMERS),
        ("hard",   PhyloPOMP.HardMERS),
        ("guided", PhyloPOMP.GuidedMERS),
    )
    agree = map(kernels) do (nm, M)
        gK = parse_newick(small_tree, demes=M.Demes)
        mK = fsmarkov(M.Camel=>0.3, M.Human=>0.7, (M.Camel,M.Human)=>1)
        seed!(101)
        pK = M.filter_pomp(gK, mK; common...)
        llest, llse = logmeanexp(
            [logLik(pfilter(pK,Np=2000)) for _ ∈ 1:15], se=true,
        )
        @info "$nm: logLik = $(round(llest,digits=2)) ± $(round(llse,sigdigits=3))"
        (nm, llest, llse)
    end
    for r ∈ agree
        @test isfinite(r[2]) && isfinite(r[3])
    end
    for i ∈ 1:length(agree), j ∈ (i+1):length(agree)
        a, b = agree[i], agree[j]
        z = abs(a[2]-b[2])/sqrt(a[3]^2+b[3]^2)
        @info "$(a[1]) vs $(b[1]): Δ = $(round(a[2]-b[2],digits=2)), z = $(round(z,digits=2))"
        @test z < 3
    end

    if heavy
        @info h2("pfilter benchmark")
        @btime pfilter($p, Np = 1000)

        @time ll = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
        llest,llse = logmeanexp(ll,se=true)
        @info "logLik = $(round(llest,digits=2)) ± $(round(llse,sigdigits=3))"
    end

end

end
