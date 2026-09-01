"""
M04 acceptance-gate tests: `reduce_event_indicator` / `reduced_transitions`
(`src/examples/mgp_reduce.jl`), the explicit reduction over the event
indicator `m` that groups M03's `full_transitions` output by REDUCED
color-only outcome and sums `phi` (exact `Rational{Int}`) within each group
to get `Phi_u`, per `mers_filter_suite.tex`:423-432.

Every concrete `(ell, n)` instance and every `phi` value below is REUSED,
unchanged, from `test/kli_full_transitions_test.jl` (M03's own hand-verified
term-by-term numbers) -- this file only adds the collapsing/summation layer
on top, per the milestone's instruction to "compute the reduced Phi_u table
by hand (sum the exact ... phi_u values you already have from M03's
term-by-term match)".
"""
module KliReduceTest

import ..Main: h1, h2

@info h1("Explicit reduction over the event indicator m: reduce_event_indicator (M04)")

using Test
using PhyloPOMP
using PhyloPOMP: full_transitions, reduce_event_indicator, reduced_transitions,
    ReducedTransition, IdentityTransition, InlineSameDemeTransition,
    CrossDemeTransition, ForkTransition, KLITransition

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

by_kind(rts) = Dict(rt.kind => rt for rt in rts)

@testset verbose=true "Explicit reduction over m (M04)" begin

    @info h2("TCC (r=(2,0)): 3 full -> 2 reduced (noop collapses" *
             " Identity+InlineSameDeme, fork stays alone)")
    @testset "TCC" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        # Same concrete instance as kli_full_transitions_test.jl "TCC":
        # phi(0)=3/10 (Identity), phi(1)=3/10 (InlineSameDeme, deme C),
        # phi(2)=1/10 (Fork).
        ts = full_transitions(tcc, [2, 0], [5, 3])
        rts = reduce_event_indicator(ts)
        @test length(rts) == 2

        g = by_kind(rts)
        @test g[:noop].Φ == 3//10 + 3//10  # == 3/5
        @test g[:noop].Φ == 3 // 5
        @test g[:fork].Φ == 1 // 10        # alone, not collapsed with anything

        # Provenance: the noop group's contributing transitions are EXACTLY
        # the Identity and InlineSameDeme instances, by type.
        noop_types = Set(typeof(t) for t in g[:noop].transitions)
        @test noop_types == Set([IdentityTransition, InlineSameDemeTransition])
        @test length(g[:noop].transitions) == 2
        fork_types = Set(typeof(t) for t in g[:fork].transitions)
        @test fork_types == Set([ForkTransition])
        @test length(g[:fork].transitions) == 1

        # Total Phi is conserved (partition, no drop/double-count).
        @test sum(rt.Φ for rt in rts) == sum(t.phi for t in ts)

        # convenience wrapper agrees (compare by (kind, Φ) pairs -- struct
        # `==` is identity-based by default since ReducedTransition carries
        # Vector fields, so a Set-of-tuples comparison is used instead of
        # `==` on the whole structs/vectors).
        rts2 = reduced_transitions(tcc, [2, 0], [5, 3])
        @test Set((rt.kind, rt.Φ) for rt in rts2) == Set((rt.kind, rt.Φ) for rt in rts)
    end

    @info h2("THH (r=(0,2)): mirror of TCC")
    @testset "THH" begin
        thh = find_event(PhyloPOMP.MERS, :transmission_hh)
        # Same instance as kli_full_transitions_test.jl "THH":
        # phi(0)=2/5, phi(1)=4/15 (InlineSameDeme, deme H), phi(2)=1/15 (Fork).
        ts = full_transitions(thh, [0, 2], [3, 6])
        rts = reduce_event_indicator(ts)
        @test length(rts) == 2

        g = by_kind(rts)
        @test g[:noop].Φ == 2//5 + 4//15
        @test g[:noop].Φ == 2 // 3
        @test g[:fork].Φ == 1 // 15

        noop_types = Set(typeof(t) for t in g[:noop].transitions)
        @test noop_types == Set([IdentityTransition, InlineSameDemeTransition])

        @test sum(rt.Φ for rt in rts) == sum(t.phi for t in ts)
    end

    @info h2("THC (r=(1,1)): 4 full -> 3 reduced -- THE critical collapse" *
             " test: (1,0) InlineSameDeme collapses with (0,0) Identity," *
             " NOT with the (0,1) CrossDeme case")
    @testset "THC" begin
        thc = find_event(PhyloPOMP.MERS, :transmission_hc)
        # Same instance as kli_full_transitions_test.jl "THC":
        # phi(0,0)=9/20 (Identity), phi(0,1)=3/20 (CrossDeme C->H),
        # phi(1,0)=3/20 (InlineSameDeme, deme C), phi(1,1)=1/20 (Fork).
        ts = full_transitions(thc, [2, 1], [5, 4])
        @test length(ts) == 4
        rts = reduce_event_indicator(ts)
        @test length(rts) == 3  # 4 full transitions -> 3 reduced

        g = by_kind(rts)
        @test length(g) == 3  # :noop, :cross, :fork -- one of each
        @test g[:noop].Φ == 9//20 + 3//20
        @test g[:noop].Φ == 3 // 5
        @test g[:cross].Φ == 3 // 20   # alone -- NOT collapsed into noop
        @test g[:fork].Φ == 1 // 20    # alone

        # The critical distinction: deme-match (not saturation index) drives
        # collapsing. (1,0)=InlineSameDeme must be IN the noop group; the
        # (0,1)=CrossDeme transition must NOT be in it.
        noop_types = Set(typeof(t) for t in g[:noop].transitions)
        @test noop_types == Set([IdentityTransition, InlineSameDemeTransition])
        noop_ss = Set(t.s for t in g[:noop].transitions)
        @test noop_ss == Set([[0, 0], [1, 0]])
        @test only(g[:cross].transitions).s == [0, 1]
        @test only(g[:fork].transitions).s == [1, 1]

        # keys distinguish the noop group from the cross group even though
        # both "involve" deme C in some sense.
        @test g[:noop].key != g[:cross].key

        @test sum(rt.Φ for rt in rts) == sum(t.phi for t in ts)
    end

    @info h2("TCH (r=(1,1)): mirror of THC, roles swapped")
    @testset "TCH" begin
        tch = find_event(PhyloPOMP.MERS, :transmission_ch)
        # Same instance as kli_full_transitions_test.jl "TCH":
        # phi(0,0)=9/20 (Identity), phi(1,0)=3/20 (CrossDeme H->C),
        # phi(0,1)=3/20 (InlineSameDeme, deme H), phi(1,1)=1/20 (Fork).
        ts = full_transitions(tch, [2, 1], [5, 4])
        rts = reduce_event_indicator(ts)
        @test length(rts) == 3

        g = by_kind(rts)
        @test g[:noop].Φ == 9//20 + 3//20
        @test g[:noop].Φ == 3 // 5
        @test g[:cross].Φ == 3 // 20
        @test g[:fork].Φ == 1 // 20

        noop_ss = Set(t.s for t in g[:noop].transitions)
        @test noop_ss == Set([[0, 0], [0, 1]])
        @test only(g[:cross].transitions).s == [1, 0]
        @test only(g[:fork].transitions).s == [1, 1]

        @test sum(rt.Φ for rt in rts) == sum(t.phi for t in ts)
    end

    @info h2("SEIR infection (r=(1,1), parent deme I): structurally" *
             " identical shape to THC/TCH, 4 full -> 3 reduced")
    @testset "SEIR infection" begin
        infection = find_event(PhyloPOMP.SEIR, :infection)
        # Same instance as kli_full_transitions_test.jl "SEIR infection":
        # phi(0,0)=2/5 (Identity), phi(1,0)=1/10 (CrossDeme I->E, new E child),
        # phi(0,1)=2/15 (InlineSameDeme, deme I, continuing parent),
        # phi(1,1)=1/30 (Fork).
        ts = full_transitions(infection, [2, 2], [6, 5])
        rts = reduce_event_indicator(ts)
        @test length(rts) == 3

        g = by_kind(rts)
        @test g[:noop].Φ == 2//5 + 2//15
        @test g[:noop].Φ == 8 // 15
        @test g[:cross].Φ == 1 // 10
        @test g[:fork].Φ == 1 // 30

        noop_ss = Set(t.s for t in g[:noop].transitions)
        @test noop_ss == Set([[0, 0], [0, 1]])
        @test only(g[:cross].transitions).s == [1, 0]
        @test only(g[:fork].transitions).s == [1, 1]

        @test sum(rt.Φ for rt in rts) == sum(t.phi for t in ts)
    end

    @info h2("SEIR progression (r=(0,1), MIGRATION): exactly 2 full and 2" *
             " reduced -- nothing collapses (no InlineSameDeme/Fork" *
             " possible, single slot always at target deme)")
    @testset "SEIR progression" begin
        progression = find_event(PhyloPOMP.SEIR, :progression)
        # Same instance as kli_full_transitions_test.jl "SEIR progression":
        # phi(0,0)=3/5 (Identity), phi(0,1)=1/5 (CrossDeme E->I).
        ts = full_transitions(progression, [3, 2], [9, 5])
        @test length(ts) == 2

        # Structural confirmation: only Identity/CrossDeme are possible --
        # no InlineSameDemeTransition or ForkTransition, since progression
        # has a single production slot (r=(0,1)) that lives at the TARGET
        # deme (I), never the ancestral deme (E) -- so sum(s) is 0 or 1
        # only, and the one occupied slot's deme (I) can never equal
        # event.from (E).
        @test all(t isa Union{IdentityTransition,CrossDemeTransition} for t in ts)
        @test !any(t isa InlineSameDemeTransition for t in ts)
        @test !any(t isa ForkTransition for t in ts)

        rts = reduce_event_indicator(ts)
        @test length(rts) == 2  # UNCHANGED from the 2 full transitions

        g = by_kind(rts)
        @test g[:noop].Φ == 3 // 5
        @test length(g[:noop].transitions) == 1  # singleton group, no collapse
        @test only(g[:noop].transitions) isa IdentityTransition
        @test g[:cross].Φ == 1 // 5
        @test length(g[:cross].transitions) == 1
        @test only(g[:cross].transitions) isa CrossDemeTransition

        @test sum(rt.Φ for rt in rts) == sum(t.phi for t in ts)

        rts2 = reduced_transitions(progression, [3, 2], [9, 5])
        @test Set((rt.kind, rt.Φ) for rt in rts2) == Set((rt.kind, rt.Φ) for rt in rts)
    end

    @info h2("generic invariants")
    @testset "invariants" begin
        # Every ReducedTransition's Φ equals the sum of its provenance
        # transitions' phi (definition, not a coincidence) -- spot-checked
        # across every case above via one representative call.
        thc = find_event(PhyloPOMP.MERS, :transmission_hc)
        ts = full_transitions(thc, [2, 1], [5, 4])
        rts = reduce_event_indicator(ts)
        for rt in rts
            @test rt.Φ == sum(t.phi for t in rt.transitions)
        end

        # Every input transition appears in exactly one output group's
        # provenance list (partition property).
        all_provenance = vcat((rt.transitions for rt in rts)...)
        @test length(all_provenance) == length(ts)
        @test Set(objectid.(all_provenance)) == Set(objectid.(ts))
    end

end

end # module
