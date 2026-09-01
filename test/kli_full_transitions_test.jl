"""
M03 acceptance-gate tests: `full_transitions` (`src/examples/mgp_transitions.jl`),
the generic classifier that turns a `production_slots`/`enumerate_saturations`/
`kli_binomial_ratio` triple (M02, `src/examples/mgp_phi.jl`) into a classified
`KLITransition` (`IdentityTransition` / `InlineSameDemeTransition` /
`CrossDemeTransition` / `ForkTransition`) per `mers_filter_suite.tex`'s Step D
(lines 459-467).

Every MERS case here is checked TERM BY TERM against
`src/examples/mers_filter_suite.tex`'s worked derivations (TCC/THH/THC/TCH,
lines 468-509, Complete Reference Table lines 541-577) -- both the exact
`Rational{Int}` phi_u value AND the operator classification. The SEIR cases
(`infection`, `progression`) have no existing tex derivation, so they are
hand-derived from the same general binomial-ratio formula (already vetted by
M02) and cross-checked against the actual, hand-coded, tested
`src/examples/seir_naive.jl` filter as a secondary oracle -- see the comments
at each SEIR testset for the exact file:line correspondence.
"""
module KliFullTransitionsTest

import ..Main: h1, h2

@info h1("Full KLI compatibility lowering: full_transitions (M03)")

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

@testset verbose=true "Full KLI compatibility lowering (M03)" begin

    @info h2("MERS transmission_cc (TCC, r=(2,0)): 3 saturations, term by term")
    @testset "TCC" begin
        # mers_filter_suite.tex:472-480, table rows 554-556.
        # I_C=5, ell_C=2 (H irrelevant, r_H=0) -- same concrete instance M02
        # already verified the phi values against (kli_phi_test.jl:90-101).
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        @test production_slots(tcc) == [2, 0]
        @test tcc.from == 1  # camel deme (I_c), the fork's source

        ts = full_transitions(tcc, [2, 0], [5, 3])
        @test length(ts) == 3

        expected = Set([
            (:identity, [0, 0], 3 // 10),  # y'=y
            (:inline,   [1, 0], 3 // 10),  # sigma^b_{CC}y = y (tex:479, "identity")
            (:fork,     [2, 0], 1 // 10),  # kappa^{bb'}_{CCC}y (tex:479-480)
        ])
        @test Set(sig(t) for t in ts) == expected

        # operator structural detail: the s_C=1 inline transition's deme is C
        # (=tcc.from); the s_C=2 fork's slot_demes are [C,C] (both slots land
        # in the camel deme) with ancestral deme C, matching tex:480's
        # "kappa^{bb'}_{CCC}" (all three subscripts camel).
        inline_t = only(filter(t -> t isa InlineSameDemeTransition, ts))
        @test inline_t.deme == 1
        fork_t = only(filter(t -> t isa ForkTransition, ts))
        @test fork_t.ancestral_deme == 1
        @test fork_t.slot_demes == [1, 1]
    end

    @info h2("MERS transmission_hh (THH, r=(0,2)): mirror of TCC")
    @testset "THH" begin
        # mers_filter_suite.tex:482-485, table rows 557-559. Note: unlike
        # M02's kli_phi_test.jl (which exercises the raw kli_binomial_ratio
        # formula directly at ell_H=1, where s_H=2 is combinatorially
        # unreachable -- min(r_H=2,ell_H=1)=1), full_transitions goes
        # through enumerate_saturations, so ell_H must be >= 2 here for the
        # s_H=2 (fork) saturation to actually be enumerated.
        thh = find_event(PhyloPOMP.MERS, :transmission_hh)
        @test production_slots(thh) == [0, 2]
        @test thh.from == 2  # human deme (I_h)

        # I_H=6, ell_H=2: phi(0)=C(4,2)/C(6,2)=6/15=2/5; phi(1)=C(4,1)/C(6,2)
        # =4/15; phi(2)=C(4,0)/C(6,2)=1/15 (Chu-Vandermonde check: C(2,0)*2/5
        # + C(2,1)*4/15 + C(2,2)*1/15 = 6/15+8/15+1/15 = 1).
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
        # mers_filter_suite.tex:487-498, table rows 560-563. C-slot is the
        # continuing camel parent (ancestral C); H-slot is the new human
        # child, ANCESTRAL DEME C (not H) -- tex:488-489.
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

        # THE distinction the milestone spec calls out most explicitly: get
        # this backwards and the whole compiler mis-derives every spillover
        # event. s=(0,1) (H-slot filled) MUST be CrossDeme (C->H); s=(1,0)
        # (C-slot filled) MUST be InlineSameDeme (C==C), never the reverse.
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
        # mers_filter_suite.tex:500-509, table rows 564-567. Exact mirror of
        # THC with C<->H swapped; the C-slot is now the child (ancestral H).
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

    @info h2("SEIR infection (r=(1,1) at parent deme I): structurally" *
             " identical to THC, hand-derived + cross-checked against" *
             " seir_naive.jl")
    @testset "SEIR infection" begin
        # SEIR demes=(E,I). infection: move=fork(I => E, I), so from=I (index
        # 2), production r=[1,1] -- one slot is the new E child (ancestral
        # deme I, the parent's deme -- tex Step A, mers_filter_suite.tex:
        # 438-441, "a new child in a different deme d' has ancestral deme d
        # (the parent's deme), not d'"), the other is the continuing-parent
        # slot which stays in I (ancestral deme I == post-event deme I).
        infection = find_event(PhyloPOMP.SEIR, :infection)
        @test production_slots(infection) == [1, 1]
        @test infection.from == 2  # I

        # n_E=6, ell_E=2, n_I=5, ell_I=2 (post-event quantities, per the tex's
        # stated convention -- mers_filter_suite.tex:454, "I_d is the
        # post-event occupancy and ell_d ... is the post-event pruned count").
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

        # Secondary oracle: seir_naive.jl's regular_part! (the tested, working
        # hand-coded filter) k==2 branch is exactly this CrossDeme case (a
        # tracked I-parent's lineage identity is reassigned to the new E
        # child): `swap!(cols,Infec,Expos,b)` at seir_naive.jl:153. Its
        # importance-weight bookkeeping there (`ll+=log(ellI)` at line 151 --
        # the -log(q) proposal-selection correction for choosing among ellI
        # tracked I-lineages -- combined with the `-log(pi[2])` already
        # subtracted at line 145, pi[2]=ellI/I) algebraically collapses,
        # together with the remaining `ll += log(1-ellI/I)-log(E)` at line
        # 156 (ellI, I, E all POST-swap/post-increment there), to exactly
        # phi_u(s=(1,0)) = [1/n_E]*[(n_I-ell_I)/n_I] -- the same product this
        # test asserts above (matching term by term once pi_u's contribution
        # and the -log(q) lineage-selection term are algebraically removed).
        # k==1 (untracked parent, seir_naive.jl:146-149) is this test's
        # Identity case with pi[1]=1-ellI/I already carrying the I-side
        # factor and `log(1-ellE/E)` (line 149) carrying the E-side factor.
        # seir_naive.jl's regular_part! never visits the InlineSameDeme or
        # Fork saturations during a REGULAR step -- see "Mathematical
        # decisions" in handoffs/M03_full_kli_lowering.md for why that is a
        # proposal-design choice (Q_u/regular-vs-singular gating, explicitly
        # out of scope for M03), not a disagreement with the phi_u values
        # asserted here.
        @test true  # documents the cross-check above; no new numeric assertion
    end

    @info h2("SEIR progression (r=(0,1), MIGRATION): resolves the migration-" *
             "semantics question against seir_naive.jl")
    @testset "SEIR progression" begin
        # move=swap(E => I): from=E (index 1), production r=[0,1] (single
        # slot, at the TARGET deme I, not the source E -- see the milestone's
        # "Migration semantics" resolution in the handoff). Hypothesis
        # confirmed: s_I=1 -> CrossDeme (swap!(y,E,I,b)); s_I=0 -> Identity;
        # there is no independent E-slot (nothing is "produced" at E).
        progression = find_event(PhyloPOMP.SEIR, :progression)
        @test production_slots(progression) == [0, 1]
        @test progression.from == 1  # E

        # n_I=5, ell_I=2 -- SAME concrete instance already verified in
        # kli_phi_test.jl:69-74 (phi(0,0)=3/5, phi(0,1)=1/5).
        # phi(s_I) = C(n_I-ell_I,1-s_I)/C(n_I,1):
        #   s_I=0: C(3,1)/C(5,1) = 3/5   (I-ell_I)/I, ell_I unchanged (untracked)
        #   s_I=1: C(3,0)/C(5,1) = 1/5   NOT ell_I/I=2/5 -- see below.
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

        # Cross-checked against seir_naive.jl:157-167 (regular_part!, k==3
        # untracked / k==4 tracked): k==3 (untracked E individual progresses)
        # adds `log(1-ellI/I)` with I POST-increment, ellI unchanged --
        # exactly (n_I-ell_I)/n_I, matching this test's Identity phi=3/5
        # formula pattern. k==4 (tracked E lineage progresses): `ll +=
        # log(ellE)` (the -log(q) proposal-selection correction, q=1/ellE)
        # then `swap!(cols,Expos,Infec,b)` at line 164 (E->I, exactly this
        # test's CrossDeme swap!(y,E,I,b)) then `ll -= log(I)` (I post-
        # increment) -- i.e. phi_u(s_I=1) = 1/n_I, EXACTLY matching this
        # test's 1/5 value, once the -log(pi[4])=-log(ellE/E) selection term
        # and the +log(ellE) proposal-correction term are algebraically
        # removed (they cancel each other's ellE dependence, leaving only
        # 1/I = 1/n_I).
        #
        # NOTE: the milestone spec's own hypothesized closed form for the
        # CrossDeme case was "phi_u = ell/I". That hypothesis is WRONG --
        # confirmed both by the general M02 binomial-ratio formula (r=(0,1)
        # gives phi(s=1)=C(n-ell,0)/C(n,1)=1/n, with NO ell-dependence) and
        # independently by seir_naive.jl's actual `ll -= log(I)` term (no
        # `ellE`/`ellI` factor survives in the final I-side contribution).
        # The correct formula, asserted above, is phi_u(s=1) = 1/I
        # (post-event I), not ell/I.
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
        # mers_filter_suite.tex:511-539 (RC/RH/SC/SH) -- explicitly deferred,
        # see ChopTransition's docstring in mgp_transitions.jl. Rejection is
        # by event TYPE (DEATH/SAMPLE/NEUTRAL), not by "r happens to be all
        # zero": MERS's sample_remove-based sampling_c/h and SEIR's
        # recovery/waning genuinely have r=(0,0), but SEIR's `sampling`
        # (built from `move=sample(I)`, not `sample_remove`) has r=(0,1) --
        # a real production slot for the sampled leaf -- and STILL must be
        # rejected, confirming the scope boundary is type-based.
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
