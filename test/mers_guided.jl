module GuidedMERSTest

import ..Main: h1, h2

@info h1("MERS model with guided proposals")

using Test
using BenchmarkTools
using Random: seed!
using PhyloPOMP
using PhyloPOMP.GuidedMERS
using PhyloPOMP.GuidedMERS.Demes: Camel, Human
import PartiallyObservedMarkovProcesses as POMP

heavy = occursin(r"y|yes|t|true", get(ENV,"RUN_HEAVY_TESTS","yes"))

## See test/mers_soft.jl for why a small reduced fixture (rather than the
## full empirical `mers_tree`) is used for the simulate/pfilter checks.
const small_tree = "(([&&PhyloPOMP deme=camel]1:1.0,[&&PhyloPOMP deme=camel]2:1.0):0.5,[&&PhyloPOMP deme=human]3:1.5);"

@testset verbose=true "MERS model with guided proposals" begin

    seed!(2121916527)

    m0 = fsmarkov(Camel=>0.3, Human=>0.7, (Camel,Human)=>1)
    p = GuidedMERS.filter_pomp(GuidedMERS.mers_tree, m0, I_c0=0, I_h0=0)
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=50))==-Inf

    @info h2("constructor validation")
    ## proposal_floor=0 is legal (disables the floor -- "no floor, let it
    ## fail"); only out-of-[0,1] values are rejected.
    p0 = GuidedMERS.filter_pomp(GuidedMERS.mers_tree, m0, proposal_floor=0.0)
    @test p0 isa POMP.PompObject
    @test_throws ArgumentError GuidedMERS.filter_pomp(GuidedMERS.mers_tree, m0, proposal_floor=-0.1)
    @test_throws ArgumentError GuidedMERS.filter_pomp(GuidedMERS.mers_tree, m0, proposal_floor=1.5)

    g = parse_newick(small_tree, demes=GuidedMERS.Demes)
    m = fsmarkov(Camel=>0.3, Human=>0.7, (Camel,Human)=>1)
    common = (
        β_hc=1.0, β_ch=1.0, χ_c=1.0, χ_h=1.0,
        S_c0=1.0, S_h0=1.0, I_c0=0.5, I_h0=0.5, N_c=100.0, N_h=100.0,
    )
    p = GuidedMERS.filter_pomp(g, m; common...)
    @test p isa POMP.PompObject

    @info h2("simulate test")
    sm = simulate(p, nsim = 3)
    @time sm = simulate(p, nsim = 3)
    @test sm isa Matrix{<:POMP.PompObject}

    @info h2("pfilter test")
    ## Np is large enough that the per-particle probability of matching the
    ## required cross-deme transmission makes total collapse to -Inf
    ## negligible regardless of seed (see test/mers_soft.jl).
    pf = pfilter(p, Np = 5000)
    @time pf = pfilter(p, Np = 5000)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(logLik(pf))

    @info h2("deterministic kernel tests")

    ## ell=0: no branches to jointly weight, so the identity outcome is
    ## certain regardless of I.
    joint0 = GuidedMERS.floored_shares([9.0], 0.05)
    @test joint0 == [1.0]

    ## I=ell boundary: the untracked-host weight is exactly 0 (naive
    ## would assign zero identity mass), but the joint floor still keeps
    ## every outcome strictly positive.
    I, ell, eps = 4, 4, 0.05
    joint = GuidedMERS.floored_shares(vcat(Float64(I-ell), zeros(ell)), eps)
    @test length(joint) == ell+1
    @test all(joint .> 0)
    @test sum(joint) ≈ 1

    ## Normalized proposal probabilities: identity + every branch sums to 1.
    I2, ell2, eps2 = 10, 3, 0.05
    rh = [0.01, 0.01, 20.0]   # one branch's relative hazard dominates
    joint2 = GuidedMERS.floored_shares(vcat(Float64(I2-ell2), rh), eps2)
    @test sum(joint2) ≈ 1
    @test all(joint2 .> 0)

    ## Guided jointly normalizes identity against the branch hazards --
    ## unlike Soft (`cross_proposal`+`floored_branch_law`), which fixes the
    ## identity/tracked split to the naive ratio (I-ell)/I regardless of
    ## the relative hazards. With one hazard dominating, Guided's identity
    ## share must fall well below the fixed naive ratio.
    naive_identity_share = (I2-ell2)/I2
    @test joint2[1] < naive_identity_share

    ## Invalid relative hazards (zero/non-finite) fall back to a uniform,
    ## strictly positive joint law rather than propagating NaN.
    joint3 = GuidedMERS.floored_shares(vcat(3.0, [0.0, NaN, Inf]), 0.05)
    @test all(isfinite.(joint3)) && all(joint3 .> 0)

    @info h2("proposal_floor=0: no floor, let it fail")
    ## epsilon=0 reduces floored_shares exactly to plain normalized shares --
    ## a weight of exactly 0 stays exactly 0, unlike the eps>0 cases above.
    w = vcat(Float64(I2-ell2), rh)
    joint4 = GuidedMERS.floored_shares(w, 0.0)
    @test joint4 ≈ w ./ sum(w)
    ## I=ell boundary (identity weight exactly 0) with nonzero branch
    ## hazards: at eps=0 the identity share is genuinely 0, not rescued by
    ## a floor -- contrast with the `joint`/`eps` case above, which was
    ## exactly this boundary with eps=0.05 and stayed strictly positive.
    joint5 = GuidedMERS.floored_shares(vcat(0.0, rh), 0.0)
    @test joint5[1] == 0
    @test all(joint5[2:end] .> 0)

    ## floored_shares_feasible (the singular_part! root/fork helper): at
    ## eps=0 an infeasible entry stays exactly 0 (as always), and a feasible
    ## entry with vanishing guide weight now genuinely gets 0 too, instead
    ## of the strictly-positive floor value it would get at eps>0.
    feas = GuidedMERS.floored_shares_feasible([0.0, 3.0, 0.0], [true,true,false], 0.0)
    @test feas[1] == 0    # feasible, but zero guide weight -- no floor to rescue it
    @test feas[2] > 0
    @test feas[3] == 0    # infeasible: always exactly 0, at any eps

    if heavy
        @info h2("pfilter benchmark")
        @btime pfilter($p, Np = 1000)

        @time ll = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
        llest,llse = logmeanexp(ll,se=true)
        @info "logLik = $(round(llest,digits=2)) ± $(round(llse,sigdigits=3))"
    end

end

end
