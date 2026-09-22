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

    @info h2("default parameters give a finite logLik")
    ## SoftMERS now shares GuidedMERS's defaults (mers_funs.jl). The
    ## pre-refactor defaults had χ_h = 0 (and β_hc = β_ch = 0), which
    ## forced -Inf at the first human tip of any mixed-deme tree; the
    ## shared defaults do not.
    ##
    ## Those defaults use N_c = N_h = 10000, which makes a 3-tip tree a
    ## severe filtering problem: a single Np=2000 run survives to the last
    ## tip only about half the time (measured over 13 seeds). The estimate
    ## over 5 replicates is -Inf only if every one of them collapses, and
    ## that is the property asserted here --- that the new defaults are not
    ## *structurally* degenerate, as the old ones were.
    seed!(1)
    pdef = SoftMERS.filter_pomp(g, m)
    lldef,_ = logmeanexp([logLik(pfilter(pdef,Np=2000)) for _ ∈ 1:5], se=true)
    @info "defaults: logLik = $(round(lldef,digits=2))"
    @test isfinite(lldef)

    @info h2("soft, hard and guided agree on the likelihood")
    ## SoftMERS, HardMERS and GuidedMERS are three importance-sampling
    ## proposals for the SAME genealogy likelihood, so their pfilter
    ## *estimates* must agree (up to Monte Carlo error), even though their
    ## per-replicate values never do.  (The test that used to live here
    ## compared per-replicate logLiks under a shared seed and asserted they
    ## differ -- true of any two distinct kernels, and so no evidence of
    ## anything.)
    ##
    ## The logmeanexp estimator is heavy-tailed on this fixture, so the
    ## replicate count was chosen for margin rather than minimality: over a
    ## 12-seed sweep at Np=2000, R=15 gave a worst-case pairwise
    ## |Δ| / sqrt(se₁²+se₂²) of 2.27, with no seed reaching 3.
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
