"""
M07 Part B acceptance-gate tests: `ProposalStrategy`/`boost`/`driver`
(`src/examples/mgp_proposal.jl`), the generic driver/boost composition
(KLI's `beta_u = alpha_u*pi_u`, `B_u = Phi_u/pi_u`) plus the
naive/soft/guided/hard `ProposalStrategy` vocabulary, per M00's decision to
KEEP the existing taxonomy rather than invent a new one.

The single required numerical cross-check: `seir_naive.jl`'s `infection`
event, `k==2` branch (tracked-parent, CrossDeme I->E), at the SAME concrete
`(ell,n)` instance M04/M05 already hand-verified
(`reduced_transitions(infection, [2,2], [6,5])` gives
`Phi(cross) = 1//10`). `seir_naive.jl:150-156`'s actual combined
log-likelihood correction for this branch (everything EXCEPT the
`-decay*step` term, which is unrelated to the boost) is reconstructed here
algebraically and compared, in EXACT `Rational{Int}` arithmetic, against
`boost(Phi_u, pi_u)` from this milestone's generic function.
"""
module KliProposalTest

import ..Main: h1, h2

@info h1("Proposal backend: ProposalStrategy / boost / driver (M07 Part B)")

using Test
using PhyloPOMP
using PhyloPOMP: ProposalStrategy, NaiveProposal, SoftProposal, GuidedProposal,
    HardProposal, boost, driver, reduced_transitions, ReducedTransition

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

@testset verbose=true "Proposal backend (M07 Part B)" begin

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

    @info h2("driver(alpha, pi) = alpha*pi -- trivial, but exact")
    @testset "driver" begin
        @test driver(3 // 4, 2 // 5) == (3 // 4) * (2 // 5)
        @test driver(0.5, 0.2) ≈ 0.1
    end

    @info h2("boost(Phi, pi) = Phi/pi -- trivial, but exact, plus the " *
             "ReducedTransition convenience overload")
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

    @info h2("Numerical cross-check against seir_naive.jl's k==2 " *
             "(tracked-parent, CrossDeme) infection branch, seir_naive.jl:150-156")
    @testset "seir_naive.jl infection k==2 cross-check" begin
        # ------------------------------------------------------------------
        # Setup: the SAME concrete instance M04/M05 already verified,
        # n_E=6, ell_E=2, n_I=5, ell_I(post-event)=2 -- reduced_transitions(
        # infection, [2,2], [6,5]) gives Phi(cross) = 1/10 (M04's table,
        # re-derived above).
        #
        # seir_naive.jl's `event_rates!` (seir_naive.jl:111,115-116) computes,
        # BEFORE this jump:
        #   pi[2] = ellI/I          -- P(tracked-parent branch proposed)
        # where `ellI` there is the PRE-event tracked-I count and `I` is the
        # (event-invariant, since infection's Delta doesn't touch I) total I
        # occupancy. Since the reduced Phi_u(cross) above is defined at the
        # POST-event ell_I=2 (per mers_filter_suite.tex Step C's convention),
        # and this branch's `swap!(cols,Infec,Expos,b)` REMOVES one tracked
        # I lineage, PRE-event ellI = POST-event ell_I + 1 = 3.
        # ------------------------------------------------------------------
        I  = 5 // 1          # n_I, unchanged by infection
        a  = 3 // 1          # ellI PRE-event (seir_naive.jl's `ellI` before swap!)
        Epost = 6 // 1        # E AFTER `E += 1` (== n_E == 6)
        ellI_post = a - 1     # == 2, matches the reduced_transitions ell_I above

        pi2 = a / I           # seir_naive.jl:116, pi[2] = ellI/I
        @test pi2 == 3 // 5

        # seir_naive.jl:145 (shared line, all k): ll -= decay*step + log(pi[k])
        # seir_naive.jl:150-156 (k==2 branch):
        #   ll += log(ellI)                      # ellI == a (still pre-swap)
        #   ellE, ellI = swap!(...)               # ellI := a-1 == ellI_post
        #   ll += log(1-ellI/I) - log(E)          # ellI == ellI_post, E == Epost
        # Combined multiplicative factor (exponentiating the ll delta,
        # EXCLUDING the -decay*step term, which has nothing to do with the
        # boost):
        naive_factor = (1 / pi2) * a * (1 - ellI_post / I) * (1 / Epost)
        @test naive_factor == 1 // 2

        # This generic milestone's boost: Phi_u(cross) / pi_u(cross), where
        # pi_u(cross) is the FULLY-marginalized proposal probability of the
        # reduced ("cross", aggregate-count-level) outcome -- i.e. pi[2]
        # (which branch) TIMES q=1/a (which specific tracked lineage is
        # chosen within that branch, seir_naive.jl:152's `rand(cols[Infec])`)
        # = (a/I)*(1/a) = 1/I. This is what seir_naive.jl's naive_factor
        # above reduces to as well (the `a`/`ellI` factors cancel), and it
        # is the general (parameter-independent) closed form: pi_u(cross) =
        # 1/n_I always, for this branch.
        pi_u_cross = 1 // I
        @test pi_u_cross == 1 // 5

        infection = find_event(PhyloPOMP.SEIR, :infection)
        rts = reduced_transitions(infection, [2, 2], [6, 5])
        cross = rts[findfirst(rt -> rt.kind == :cross, rts)]
        @test cross.Φ == 1 // 10

        B_u = boost(cross, pi_u_cross)
        @test B_u == 1 // 2
        @test B_u == naive_factor   # <-- the required cross-check
    end

end

end # module KliProposalTest
