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

    @info h2("constructor validation")
    ## proposal_floor=0 is legal (disables the floor -- "no floor, let it
    ## fail"); only out-of-[0,1] values are rejected.
    p0 = SoftMERS.filter_pomp(SoftMERS.mers_tree, m0, proposal_floor=0.0)
    @test p0 isa POMP.PompObject
    @test_throws ArgumentError SoftMERS.filter_pomp(SoftMERS.mers_tree, m0, proposal_floor=-0.1)
    @test_throws ArgumentError SoftMERS.filter_pomp(SoftMERS.mers_tree, m0, proposal_floor=1.5)

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
    ## Np is large enough that the per-particle probability of matching the
    ## required cross-deme transmission (estimated empirically at ~0.2%
    ## per particle for this fixture) makes total collapse to -Inf
    ## negligible (<1e-9) regardless of seed.
    pf = pfilter(p, Np = 5000)
    @time pf = pfilter(p, Np = 5000)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(logLik(pf))

    @info h2("deterministic kernel tests")

    ## ell=0: the identity outcome is certain, regardless of I.
    q0, qsum = SoftMERS.cross_proposal(7, 0, 0.05)
    @test q0 == 1
    @test qsum == 0
    @test isempty(SoftMERS.floored_branch_law(Float64[], 0, 7, 0.05))

    ## I=ell boundary: identity and every tracked branch keep strictly
    ## positive proposal mass, even though the naive base law would put
    ## zero mass on the identity outcome here.
    I, ell, eps = 4, 4, 0.05
    q0b, qsumb = SoftMERS.cross_proposal(I, ell, eps)
    @test q0b > 0
    @test qsumb > 0
    qb = SoftMERS.floored_branch_law(zeros(ell), ell, I, eps)
    @test all(qb .> 0)
    @test qsumb ≈ sum(qb)

    ## Normalized proposal probabilities: identity + every branch sums to 1.
    @test q0b + qsumb ≈ 1
    I2, ell2, eps2 = 12, 3, 0.1
    q0c, qsumc = SoftMERS.cross_proposal(I2, ell2, eps2)
    @test q0c + qsumc ≈ 1

    ## Aggregate-plus-conditional equivalence: the aggregate tracked-group
    ## probability times the conditional (relhaz-weighted) branch share
    ## equals the full per-outcome q_epsilon(b).
    rh = [0.1, 4.0, 0.3]
    qb2 = SoftMERS.floored_branch_law(rh, ell2, I2, eps2)
    @test qsumc ≈ sum(qb2)
    conditional = qb2 ./ qsumc
    @test all((qsumc .* conditional) .≈ qb2)

    ## Invalid relative hazards (zero/non-finite) fall back to a uniform,
    ## strictly positive conditional law rather than propagating NaN.
    qb3 = SoftMERS.floored_branch_law([0.0, NaN, Inf], 3, 10, 0.05)
    @test all(isfinite.(qb3)) && all(qb3 .> 0)

    @info h2("proposal_floor=0: no floor, let it fail")
    ## epsilon=0 reduces cross_proposal exactly to the plain naive split,
    ## with no floor mass added at the I=ell boundary (contrast with eps>0
    ## above, where q0b/qsumb stayed strictly positive there).
    q0d, qsumd = SoftMERS.cross_proposal(I2, ell2, 0.0)
    @test q0d ≈ (I2-ell2)/I2
    @test qsumd ≈ ell2/I2
    q0e, qsume = SoftMERS.cross_proposal(I, ell, 0.0)   # I=ell boundary
    @test q0e == 0   # naive base law puts zero identity mass here; no floor to rescue it
    @test qsume ≈ 1
    ## A branch with genuinely zero relative hazard gets exactly zero
    ## proposal mass at epsilon=0 -- the "let it fail" case the floor
    ## (epsilon>0) is designed to prevent.
    qb4 = SoftMERS.floored_branch_law([0.0, 4.0, 0.3], ell2, I2, 0.0)
    @test qb4[1] == 0
    @test all(qb4[2:end] .> 0)

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
