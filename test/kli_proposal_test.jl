"""
Tests for ProposalStrategy, boost and driver: driver = alpha*pi,
boost = Phi/pi. boost(Phi(cross), pi) is compared with seir_naive.jl's
tracked-parent infection branch in exact Rational arithmetic.
"""
module KliProposalTest

import ..Main: h1, h2

@info h1("Proposal backend: ProposalStrategy / boost / driver")

using Test
using PhyloPOMP
using PhyloPOMP: ProposalStrategy, NaiveProposal, SoftProposal, GuidedProposal,
    HardProposal, boost, driver, reduced_transitions, ReducedTransition

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

@testset verbose=true "Proposal backend" begin

    @info h2("ProposalStrategy taxonomy: four disjoint singleton tags")
    @testset "ProposalStrategy taxonomy" begin
        @test NaiveProposal <: ProposalStrategy
        @test SoftProposal <: ProposalStrategy
        @test GuidedProposal <: ProposalStrategy
        @test HardProposal <: ProposalStrategy
        # Distinct, instantiable singleton types.
        tags = ProposalStrategy[NaiveProposal(), SoftProposal(), GuidedProposal(), HardProposal()]
        @test length(unique(typeof, tags)) == 4
    end

    @info h2("driver(alpha, pi) = alpha*pi")
    @testset "driver" begin
        @test driver(3 // 4, 2 // 5) == (3 // 4) * (2 // 5)
        @test driver(0.5, 0.2) ≈ 0.1
    end

    @info h2("boost(Phi, pi) = Phi/pi, plus the ReducedTransition overload")
    @testset "boost" begin
        @test boost(1 // 10, 1 // 5) == 1 // 2
        @test boost(3 // 20, 1 // 4) == 3 // 5

        infection = find_event(PhyloPOMP.SEIR, :infection)
        rts = reduced_transitions(infection, [2, 2], [6, 5])
        cross = rts[findfirst(rt -> rt.kind == :cross, rts)]
        @test cross.Φ == 1 // 10
        @test boost(cross, 1 // 5) == boost(cross.Φ, 1 // 5)
        @test boost(cross, 1 // 5) == 1 // 2
    end

    @info h2("Cross-check against seir_naive.jl's tracked-parent (CrossDeme) infection branch")
    @testset "seir_naive.jl infection k==2 cross-check" begin
        # n_E=6, ell_E=2, n_I=5, post-event ell_I=2; Phi(cross)=1/10.
        # Pre-event tracked I count is 3, since the swap removes one.
        I  = 5 // 1          # n_I, unchanged by infection
        a  = 3 // 1          # ellI PRE-event (seir_naive.jl's `ellI` before swap!)
        Epost = 6 // 1        # E AFTER `E += 1` (== n_E == 6)
        ellI_post = a - 1     # == 2, matches the reduced_transitions ell_I above

        pi2 = a / I
        @test pi2 == 3 // 5

        # Naive correction factor = (1/pi2) * a * (1 - ellI_post/I) * (1/E), without the decay term.
        naive_factor = (1 / pi2) * a * (1 - ellI_post / I) * (1 / Epost)
        @test naive_factor == 1 // 2

        # pi_u(cross) = pi2 * (1/a) = 1/n_I; the ell factors cancel.
        pi_u_cross = 1 // I
        @test pi_u_cross == 1 // 5

        infection = find_event(PhyloPOMP.SEIR, :infection)
        rts = reduced_transitions(infection, [2, 2], [6, 5])
        cross = rts[findfirst(rt -> rt.kind == :cross, rts)]
        @test cross.Φ == 1 // 10

        B_u = boost(cross, pi_u_cross)
        @test B_u == 1 // 2
        @test B_u == naive_factor
    end

end

end # module KliProposalTest
