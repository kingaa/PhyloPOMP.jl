## Contributed by Simon Frost (@sdwfrost) in Issue #8 (https://github.com/kingaa/PhyloPOMP.jl/issues/8).
## He reports:
## The exact log likelihood is −3.1014.
## I computed it by solving the filter equation directly on the finite joint space of (S,E,I,R) and
## lineage colourings (uniformization), cross-checked to 1e-15 against an independent matrix-exponential
## implementation.

module ExactTest

import ..Main: h1, h2

@info h1("SEIR filters vs known likelihood")

using Test
using PhyloPOMP, PhyloPOMP.NaiveSEIR, PhyloPOMP.SoftSEIR, PhyloPOMP.GuidedSEIR, PhyloPOMP.HardSEIR
using Random: seed!
seed!(641922521)

@testset verbose=true "SEIR filters vs known likelihood" begin

    ## a single infective carries the root lineage,
    ## so transmissions occur with I == ellI
    g = parse_newick("(([&&PhyloPOMP type=sample]:0.9)[&&PhyloPOMP type=sample]:0.6);", time=1.8)
    par = (
        β=4.0, σ=0.7, γ=0.5, ω=0.2, ψ=0.4, χ=0.0,
        pop=6, S0=4/6, E0=0.0, I0=1/6, R0=1/6,
    )
    m(M) = fsmarkov(
        M.Demes.Expos=>0.1,
        M.Demes.Infec=>1,
        (M.Demes.Expos,M.Demes.Infec)=>1
    )
    est(p) = logmeanexp([logLik(pfilter(p, Np=500)) for _ in 1:200], se=true)

    ll,se = est(NaiveSEIR.filter_pomp(g; par...))
    @info "naïve filter estimate = $(round(ll,digits=2)) ± $(round(se,sigdigits=2))"
    @test isapprox(ll,-3.1014; atol=3*se)

    ll,se = est(SoftSEIR.filter_pomp(g, m(SoftSEIR); par...))
    @info "soft filter estimate = $(round(ll,digits=2)) ± $(round(se,sigdigits=2))"
    @test isapprox(ll,-3.1014; atol=3*se)

    ll, se = est(GuidedSEIR.filter_pomp(g, m(GuidedSEIR); par...))
    @info "semisoft filter estimate = $(round(ll,digits=2)) ± $(round(se,sigdigits=2))"
    @test isapprox(ll,-3.1014; atol=3*se)

    ll,se = est(HardSEIR.filter_pomp(g, m(HardSEIR); par...))
    @info "hard filter estimate = $(round(ll,digits=2)) ± $(round(se,sigdigits=2))"
    @test isapprox(ll,-3.1014; atol=3*se)

end

end
