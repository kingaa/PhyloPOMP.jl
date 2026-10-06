"""
Tests for full_transitions: each saturation is classified Identity,
InlineSameDeme, CrossDeme or Fork with an exact Rational phi. MERS cases
follow the worked table in mers_filter_suite.tex; SEIR cases use the same
binomial ratio.
"""
module KliFullTransitionsTest

import ..Main: h1, h2

@info h1("Full KLI compatibility lowering: full_transitions")

using Test
using PhyloPOMP
using PhyloPOMP: production_slots, enumerate_saturations, kli_binomial_ratio,
    full_transitions, IdentityTransition, InlineSameDemeTransition,
    CrossDemeTransition, ForkTransition, ChopTransition, KLITransition

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

# classification-signature helper: (kind symbol, s, phi) for order-independent
# comparison against a hand-built expected set.
kind(t::IdentityTransition) = :identity
kind(t::InlineSameDemeTransition) = :inline
kind(t::CrossDemeTransition) = :cross
kind(t::ForkTransition) = :fork
sig(t::KLITransition) = (kind(t), t.s, t.phi)

@testset verbose=true "Full KLI compatibility lowering" begin

    @info h2("MERS transmission_cc (TCC, r=(2,0)): 3 saturations, term by term")
    @testset "TCC" begin
        # TCC, I_C=5, ell_C=2. The H deme is irrelevant (r_H=0).
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        @test production_slots(tcc) == [2, 0]
        @test tcc.from == 1  # camel deme (I_c), the fork's source

        ts = full_transitions(tcc, [2, 0], [5, 3])
        @test length(ts) == 3

        expected = Set([
            (:identity, [0, 0], 3 // 10),  # y'=y
            (:inline,   [1, 0], 3 // 10),  # sigma^b_{CC}y = y, "identity"
            (:fork,     [2, 0], 1 // 10),  # kappa^{bb'}_{CCC}y
        ])
        @test Set(sig(t) for t in ts) == expected

        # Inline slot is in deme C; both fork slots are in C.
        inline_t = only(filter(t -> t isa InlineSameDemeTransition, ts))
        @test inline_t.deme == 1
        fork_t = only(filter(t -> t isa ForkTransition, ts))
        @test fork_t.ancestral_deme == 1
        @test fork_t.slot_demes == [1, 1]
    end

    @info h2("MERS transmission_hh (THH, r=(0,2)): mirror of TCC")
    @testset "THH" begin
        # THH needs ell_H >= 2 so the fork saturation is enumerated.
        thh = find_event(PhyloPOMP.MERS, :transmission_hh)
        @test production_slots(thh) == [0, 2]
        @test thh.from == 2  # human deme (I_h)

        # I_H=6, ell_H=2: phi = 2/5, 4/15, 1/15; the ell-weighted sum is 1.
        ts = full_transitions(thh, [0, 2], [3, 6])
        @test length(ts) == 3
        expected = Set([
            (:identity, [0, 0], 2 // 5),
            (:inline,   [0, 1], 4 // 15),
            (:fork,     [0, 2], 1 // 15),
        ])
        @test Set(sig(t) for t in ts) == expected

        inline_t = only(filter(t -> t isa InlineSameDemeTransition, ts))
        @test inline_t.deme == 2
        fork_t = only(filter(t -> t isa ForkTransition, ts))
        @test fork_t.ancestral_deme == 2
        @test fork_t.slot_demes == [2, 2]
    end

    @info h2("MERS transmission_hc (THC, r=(1,1)): all 4 saturations, the" *
             " (0,1)-CrossDeme vs (1,0)-InlineSameDeme distinction checked" *
             " explicitly")
    @testset "THC" begin
        # C slot is the continuing parent; H slot is the new child, whose ancestral deme is C.
        thc = find_event(PhyloPOMP.MERS, :transmission_hc)
        @test production_slots(thc) == [1, 1]
        @test thc.from == 1  # camel deme (I_c) is the parent

        ts = full_transitions(thc, [2, 1], [5, 4])
        @test length(ts) == 4
        expected = Set([
            (:identity, [0, 0], 9 // 20),  # y'=y
            (:cross,    [0, 1], 3 // 20),  # sigma^b_{CH}y, b in col_C -- NEW
                                            # human child traces to camel
                                            # parent: ancestral C != post-event H
            (:inline,   [1, 0], 3 // 20),  # sigma^b_{CC}y = y -- continuing
                                            # parent: ancestral C == post-event C
            (:fork,     [1, 1], 1 // 20),  # kappa^{bb'}_{CHC}y
        ])
        @test Set(sig(t) for t in ts) == expected

        # (0,1) must be CrossDeme (C to H) and (1,0) InlineSameDeme, never the reverse.
        cross_t = only(filter(t -> t isa CrossDemeTransition, ts))
        @test cross_t.s == [0, 1]
        @test cross_t.ancestral_deme == 1  # C
        @test cross_t.target_deme == 2     # H
        inline_t = only(filter(t -> t isa InlineSameDemeTransition, ts))
        @test inline_t.s == [1, 0]
        @test inline_t.deme == 1  # C

        fork_t = only(filter(t -> t isa ForkTransition, ts))
        @test fork_t.ancestral_deme == 1        # C
        @test Set(fork_t.slot_demes) == Set([1, 2])  # {C, H}, tex's kappa^{bb'}_{CHC}
    end

    @info h2("MERS transmission_ch (TCH, r=(1,1)): mirrored roles")
    @testset "TCH" begin
        # Mirror of THC with C and H swapped.
        tch = find_event(PhyloPOMP.MERS, :transmission_ch)
        @test production_slots(tch) == [1, 1]
        @test tch.from == 2  # human deme (I_h) is the parent

        ts = full_transitions(tch, [2, 1], [5, 4])
        @test length(ts) == 4
        expected = Set([
            (:identity, [0, 0], 9 // 20),
            (:cross,    [1, 0], 3 // 20),  # sigma^b_{HC}y, b in col_H -- NEW
                                            # camel child from human parent
            (:inline,   [0, 1], 3 // 20),  # sigma^b_{HH}y = y
            (:fork,     [1, 1], 1 // 20),  # kappa^{bb'}_{HCH}y
        ])
        @test Set(sig(t) for t in ts) == expected

        cross_t = only(filter(t -> t isa CrossDemeTransition, ts))
        @test cross_t.s == [1, 0]
        @test cross_t.ancestral_deme == 2  # H
        @test cross_t.target_deme == 1     # C
        inline_t = only(filter(t -> t isa InlineSameDemeTransition, ts))
        @test inline_t.s == [0, 1]
        @test inline_t.deme == 2  # H

        fork_t = only(filter(t -> t isa ForkTransition, ts))
        @test fork_t.ancestral_deme == 2  # H
        @test Set(fork_t.slot_demes) == Set([1, 2])
    end

    @info h2("SEIR infection (r=(1,1) at parent deme I): same shape as THC")
    @testset "SEIR infection" begin
        # infection: one slot is the new E child (ancestral deme I), the other the continuing parent in I.
        infection = find_event(PhyloPOMP.SEIR, :infection)
        @test production_slots(infection) == [1, 1]
        @test infection.from == 2  # I

        # n_E=6, ell_E=2, n_I=5, ell_I=2; phi = product of per-deme ratios:
        # phi(s_E,s_I) = C(n_E-ell_E,1-s_E)/C(n_E,1) * C(n_I-ell_I,1-s_I)/C(n_I,1):
        #   (0,0): (4/6)*(3/5) = 2/5
        #   (1,0): (1/6)*(3/5) = 1/10
        #   (0,1): (4/6)*(1/5) = 2/15
        #   (1,1): (1/6)*(1/5) = 1/30
        ts = full_transitions(infection, [2, 2], [6, 5])
        @test length(ts) == 4
        expected = Set([
            (:identity, [0, 0], 2 // 5),
            (:cross,    [1, 0], 1 // 10),  # new E child tracked: ancestral I != E
            (:inline,   [0, 1], 2 // 15),  # continuing parent tracked: I == I
            (:fork,     [1, 1], 1 // 30),
        ])
        @test Set(sig(t) for t in ts) == expected

        cross_t = only(filter(t -> t isa CrossDemeTransition, ts))
        @test cross_t.ancestral_deme == 2  # I
        @test cross_t.target_deme == 1     # E
        inline_t = only(filter(t -> t isa InlineSameDemeTransition, ts))
        @test inline_t.deme == 2  # I
    end

    @info h2("SEIR progression (r=(0,1), MIGRATION)")
    @testset "SEIR progression" begin
        # progression: single slot at the target deme I, so s_I=1 is CrossDeme and s_I=0 is Identity.
        progression = find_event(PhyloPOMP.SEIR, :progression)
        @test production_slots(progression) == [0, 1]
        @test progression.from == 1  # E

        # n_I=5, ell_I=2: phi(s_I=0)=3/5, phi(s_I=1)=1/5, with no ell dependence in the second.
        ts = full_transitions(progression, [3, 2], [9, 5])
        @test length(ts) == 2
        expected = Set([
            (:identity, [0, 0], 3 // 5),
            (:cross,    [0, 1], 1 // 5),
        ])
        @test Set(sig(t) for t in ts) == expected

        cross_t = only(filter(t -> t isa CrossDemeTransition, ts))
        @test cross_t.ancestral_deme == 1  # E
        @test cross_t.target_deme == 2     # I

        # phi(s_I=1) = 1/n_I, with no ell dependence (r=(0,1) gives C(n-ell,0)/C(n,1)).
    end

    @info h2("boundary: ell=0 forces only the Identity transition to be" *
             " feasible")
    @testset "boundary ell=0" begin
        progression = find_event(PhyloPOMP.SEIR, :progression)
        # ell_I=0: min(r_I=1, ell_I=0)=0, so enumerate_saturations forces
        # s_I=0 -- the CrossDeme outcome cannot even be ENUMERATED (no
        # tracked I-lineage exists to move), let alone realized.
        ts = full_transitions(progression, [0, 0], [9, 5])
        @test length(ts) == 1
        @test ts[1] isa IdentityTransition
        @test ts[1].s == [0, 0]
        @test ts[1].phi == 1 // 1  # C(5,1)/C(5,1) = 1, no tracked lineages to miss
    end

    @info h2("boundary: ell=n forces phi=0 on the untracked-outcome" *
             " saturation (still enumerated, but assigned zero weight)")
    @testset "boundary ell=n" begin
        progression = find_event(PhyloPOMP.SEIR, :progression)
        # ell_I = n_I = 5 (deme I fully tracked): s_I=0 is still COMBINATORIALLY
        # enumerated (min(1,5)=1, so {0,1}), but phi(s_I=0) = C(0,1)/C(5,1) = 0
        # -- forced infeasible, since there is no untracked I-individual left
        # for an untracked lineage to have progressed. Only s_I=1 (CrossDeme)
        # carries positive weight.
        ts = full_transitions(progression, [3, 5], [9, 5])
        @test length(ts) == 2
        by_kind = Dict(kind(t) => t for t in ts)
        @test by_kind[:identity].phi == 0 // 1
        @test by_kind[:cross].phi == 1 // 5

        # richer illustration on MERS THC: fully-tracked camel deme (ell_C=
        # n_C=5) forces BOTH s_C=0 outcomes (Identity and the H-only CrossDeme)
        # to phi=0 -- there is no untracked camel individual left to be the
        # (untracked) parent -- while both s_C=1 outcomes (InlineSameDeme,
        # Fork) retain positive weight.
        thc = find_event(PhyloPOMP.MERS, :transmission_hc)
        ts2 = full_transitions(thc, [5, 1], [5, 4])
        @test length(ts2) == 4
        by_kind2 = Dict(kind(t) => t for t in ts2)
        @test by_kind2[:identity].phi == 0 // 1
        @test by_kind2[:cross].phi == 0 // 1
        @test by_kind2[:inline].phi == 3 // 20
        @test by_kind2[:fork].phi == 1 // 20
    end

    @info h2("scope: full_transitions rejects DEATH/SAMPLE/NEUTRAL events" *
             " (r=(0,0) always; decay/rate-driven, not saturation-driven)")
    @testset "out-of-scope event types" begin
        # Events are rejected by type (DEATH/SAMPLE/NEUTRAL), not by r.
        # SEIR sampling has r=(0,1) and is still rejected.
        for (model, name, expected_r) in (
            (PhyloPOMP.SEIR, :recovery, [0, 0]), (PhyloPOMP.SEIR, :waning, [0, 0]),
            (PhyloPOMP.SEIR, :sampling, [0, 1]), (PhyloPOMP.MERS, :removal_c, [0, 0]),
            (PhyloPOMP.MERS, :sampling_c, [0, 0]), (PhyloPOMP.MERS, :birth_c, [0, 0]),
            (PhyloPOMP.MERS, :death_c, [0, 0]),
        )
            ev = find_event(model, name)
            @test production_slots(ev) == expected_r
            @test_throws ArgumentError full_transitions(ev, zeros(Int, length(model.demes)),
                                                          ones(Int, length(model.demes)))
        end

        # ChopTransition exists (IR vocabulary is complete) but is never
        # constructed by full_transitions -- it is a documented stub type.
        recovery = find_event(PhyloPOMP.SEIR, :recovery)
        @test ChopTransition(recovery) isa KLITransition
        @test ChopTransition(recovery) isa ChopTransition
    end

end

end # module
