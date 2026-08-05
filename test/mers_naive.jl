module NaiveMERSTest

import ..Main: h1, h2

@info h1("MERS model with naïve proposals")

using Test
using BenchmarkTools
using Random: seed!
using PhyloPOMP
using PhyloPOMP.NaiveMERS
import PartiallyObservedMarkovProcesses as POMP

heavy = occursin(r"y|yes|t|true", get(ENV,"RUN_HEAVY_TESTS","yes"))

## Unlike Soft/Guided/Hard, NaiveMERS.filter_pomp takes no genealogy argument
## and no `proposal_floor` -- it always filters the full empirical `mers_tree`
## (548 nodes) with R phylopomp-style parameter names. There is no way to
## substitute a small fixture here, so the tests below work directly against
## the real tree.
@testset verbose=true "MERS model with naïve proposals" begin

    seed!(2121916527)

    p = NaiveMERS.filter_pomp(I_c0=0,I_h0=0)
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=50))==-Inf

    @info h2("constructor validation")
    ## filter_pomp accepts both the Ic0/Ih0 spelling (matching R phylopomp)
    ## and the underscored I_c0/I_h0 alias; passing both with conflicting
    ## values is rejected rather than silently picking one.
    @test_throws ArgumentError NaiveMERS.filter_pomp(Ic0=1.0, I_c0=2.0)
    @test_throws ArgumentError NaiveMERS.filter_pomp(Ih0=1.0, I_h0=2.0)

    p = NaiveMERS.filter_pomp()
    @test p isa POMP.PompObject

    @info h2("simulate test")
    sm = simulate(p, nsim = 3)
    @time sm = simulate(p, nsim = 3)
    @test sm isa Matrix{<:POMP.PompObject}

    @info h2("pfilter test")
    ## The naive kernel has no guide, so on the full ~274-tip tree the
    ## per-particle probability of matching every observed tip's deme is
    ## astronomically small -- this is precisely why the Soft/Guided/Hard
    ## kernels exist. logLik here is expected to be -Inf regardless of Np;
    ## we only check that pfilter runs and returns the right type.
    ## Separately, under current defaults (chi_h = 0.0) the first human
    ## Sample node deterministically forces ll = log(chi_h*Ih) = -Inf, so
    ## -Inf fires immediately (node 11 of 548) rather than requiring the
    ## particle filter to fail to find a compatible coloring.
    pf = pfilter(p, Np = 200)
    @time pf = pfilter(p, Np = 200)
    @test pf isa POMP.PfilterdPompObject

end

end
