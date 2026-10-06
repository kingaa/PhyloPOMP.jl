"""
`simulate` on `@mgp MTBD`, then `GuidedMTBD.filter_pomp`.
With `r1 < 1` some samples keep their lineage. The MTBD filters need `r1 = r2 = 1`.
"""
module MTBDSimulateTest

## `h1`/`h2` come from runtests.jl; fall back when run standalone.
const h1 = isdefined(Main, :h1) ? Main.h1 : identity
const h2 = isdefined(Main, :h2) ? Main.h2 : identity
isdefined(Main, :SimulateChecks) || Base.include(Main, joinpath(@__DIR__, "simulate_checks.jl"))

@info h1("MTBD forward simulator")

using Test
using Random: MersenneTwister, seed!
using PhyloPOMP
using PhyloPOMP: MTBDDemes
using Main.SimulateChecks
using PhyloPOMP.GuidedMTBD
import PartiallyObservedMarkovProcesses as POMP

## Rates as `GuidedMTBD.filter_pomp` takes them; `simulate` also needs r1, r2.
const rates = (lambda11=1.2, lambda12=0.3, lambda21=0.2, lambda22=0.9, m12=0.2, m21=0.1,
               mu1=0.5, mu2=0.5, psi1=0.3, psi2=0.3)

@testset verbose=true "MTBD forward simulator" begin

    @info h2("invariants: one type-1 founder, and a forest with both types")
    θ = merge(rates, (r1 = 0.7, r2 = 1.0))
    cov = run_invariants(PhyloPOMP.MTBD, θ, [((I1=1, I2=0), [1, 0]), ((I1=2, I2=1), [2, 1])];
                         ntree = 100, tmax = 5.0, demeset = MTBDDemes,
                         samplemap = [MTBDDemes.I1, MTBDDemes.I2])
    @test cov.nonempty > 50
    @test cov.inline > 0        # r1 < 1: some type-1 samples keep their lineage

    @testset "simulate -> GuidedMTBD (r1 = r2 = 1): check passes and logLik finite" begin
        @info h2("simulate -> GuidedMTBD (r1 = r2 = 1): check passes and logLik finite")
        θ1 = merge(rates, (r1 = 1.0, r2 = 1.0))
        sim(s) = simulate(PhyloPOMP.MTBD, θ1; x0 = (I1=1, I2=0), graft = [1, 0], tmax = 5.0,
                          rng = MersenneTwister(s), demeset = GuidedMTBD.Demes,
                          samplemap = [GuidedMTBD.Camel, GuidedMTBD.Human])
        ## First seed giving 3 to 15 samples from both types.
        g = nothing
        for s ∈ 1:500
            gs = sim(s)
            d = [gs[i].deme for i ∈ samples(gs)]
            if 3 ≤ nsample(gs) ≤ 15 && GuidedMTBD.Camel ∈ d && GuidedMTBD.Human ∈ d
                g = gs; break
            end
        end
        @test !isnothing(g)
        @test GuidedMTBD.check(g) === nothing
        m = fsmarkov(GuidedMTBD.Camel=>0.5, GuidedMTBD.Human=>0.5,
                     (GuidedMTBD.Camel,GuidedMTBD.Human)=>3.0)
        p = GuidedMTBD.filter_pomp(g, m; rates...)
        @test p isa POMP.PompObject
        seed!(20261005)
        @test isfinite(logLik(pfilter(p, Np=1000)))
    end

end

end
