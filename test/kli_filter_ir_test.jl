"""
M05 acceptance-gate tests: `classify_filter_terms` / `filter_spec`
(`src/examples/mgp_filter_ir.jl`), the Filter IR classification pass that
turns M04's `Vector{ReducedTransition}` into `RegularFlow`/`SingularFlow`/
`OutflowImbalanceTerm` buckets, per `mers_filter_suite.tex`'s driver/boost/decay
sections (cited in full in `mgp_filter_ir.jl`'s header comment).

Every concrete `(ell, n)` instance and every `Φ` value below is REUSED,
unchanged, from `test/kli_reduce_test.jl` (M04's own hand-verified numbers,
which are themselves M03's `phi` numbers, summed) -- this file only adds the
regular/singular/decay classification layer on top, and re-derives the
numbers via the actual function call chain rather than trusting the handoff
tables.
"""
module KliFilterIrTest

import ..Main: h1, h2

@info h1("Filter IR: classify_filter_terms / filter_spec (M05)")

using Test
using PhyloPOMP
using PhyloPOMP: full_transitions, reduce_event_indicator, reduced_transitions,
    ReducedTransition, KLITransition,
    FilterTerm, RegularFlow, SingularFlow, OutflowImbalanceTerm, FilterSpec,
    classify_filter_terms, filter_spec,
    Event, EventType, BIRTH, MIGRATION, SAMPLE, kli_decay

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

@testset verbose=true "Filter IR classification (M05)" begin

    @info h2("SEIR infection (regular BIRTH): noop+cross -> RegularFlow," *
             " fork -> decay bucket")
    @testset "SEIR infection" begin
        infection = find_event(PhyloPOMP.SEIR, :infection)
        @test infection.regular == true

        # Same instance as kli_reduce_test.jl "SEIR infection":
        # Phi(noop)=8/15, Phi(cross)=1/10, Phi(fork)=1/30.
        rts = reduced_transitions(infection, [2, 2], [6, 5])
        @test length(rts) == 3

        spec = classify_filter_terms(infection, rts)
        @test spec isa FilterSpec

        # Exactly 2 RegularFlow terms (noop, cross).
        @test length(spec.regular) == 2
        reg_by_kind = Dict(rf.reduced.kind => rf for rf in spec.regular)
        @test Set(keys(reg_by_kind)) == Set([:noop, :cross])
        @test reg_by_kind[:noop].Φ == 8 // 15
        @test reg_by_kind[:cross].Φ == 1 // 10

        # Exactly 1 non-regular term (fork), landing in the decay bucket.
        @test isempty(spec.singular)
        @test length(spec.outflow_imbalance) == 1
        @test spec.outflow_imbalance[1].reduced.kind == :fork
        @test spec.outflow_imbalance[1].Φ == 1 // 30
        @test spec.outflow_imbalance[1].reason == :fork_unobserved_at_regular_time

        # convenience wrapper agrees
        spec2 = filter_spec(infection, [2, 2], [6, 5])
        @test length(spec2.regular) == 2
        @test length(spec2.outflow_imbalance) == 1
        @test isempty(spec2.singular)
    end

    @info h2("SEIR progression (regular MIGRATION): both reduced" *
             " transitions -> RegularFlow (no fork/decay case exists)")
    @testset "SEIR progression" begin
        progression = find_event(PhyloPOMP.SEIR, :progression)
        @test progression.regular == true

        # Same instance as kli_reduce_test.jl "SEIR progression":
        # Phi(noop)=3/5, Phi(cross)=1/5.
        rts = reduced_transitions(progression, [3, 2], [9, 5])
        @test length(rts) == 2
        @test Set(rt.kind for rt in rts) == Set([:noop, :cross])  # no :fork

        spec = classify_filter_terms(progression, rts)
        @test length(spec.regular) == 2
        @test isempty(spec.singular)
        @test isempty(spec.outflow_imbalance)  # nothing to defer -- matches
        # seir_naive.jl:157-167's regular_part! k==3/k==4 branches, which
        # handle exactly the untracked-E (Identity/noop) and tracked-E
        # (CrossDeme) cases and nothing else.

        reg_by_kind = Dict(rf.reduced.kind => rf for rf in spec.regular)
        @test reg_by_kind[:noop].Φ == 3 // 5
        @test reg_by_kind[:cross].Φ == 1 // 5
    end

    @info h2("MERS transmission_cc (TCC, regular BIRTH): 1 RegularFlow" *
             " (noop), 1 decay-bucket term (fork) -- no cross case exists")
    @testset "MERS transmission_cc" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        @test tcc.regular == true

        # Same instance as kli_reduce_test.jl "TCC": Phi(noop)=3/5,
        # Phi(fork)=1/10.
        rts = reduced_transitions(tcc, [2, 0], [5, 3])
        @test length(rts) == 2

        spec = classify_filter_terms(tcc, rts)
        @test length(spec.regular) == 1
        @test spec.regular[1].reduced.kind == :noop
        @test spec.regular[1].Φ == 3 // 5

        @test isempty(spec.singular)
        @test length(spec.outflow_imbalance) == 1
        @test spec.outflow_imbalance[1].reduced.kind == :fork
        @test spec.outflow_imbalance[1].Φ == 1 // 10
        @test spec.outflow_imbalance[1].reason == :fork_unobserved_at_regular_time

        # Matches mers_naive.jl's regular_part! k==1 branch (mers_naive.jl:
        # ~180-183): "Sc -= 1; Ic += 1; ll += log(1-(ellc*(ellc-1)/Ic/(Ic-1)))"
        # -- no `fork!`/`swap!` call at all, i.e. only the noop outcome is
        # ever realized at regular time. Note `Φ(noop)+Φ(fork) = 7/10 ≠ 1`
        # here (this instance's `ℓ_C=2, I_C=5` is not the `ellc*(ellc-1)/
        # (Ic*(Ic-1))` bracket's own normalization point) -- `Σ_s φ_u(s)`
        # is not in general a probability distribution over `s` (it is a
        # sum of binomial-ratio compatibility factors, not a categorical
        # kernel), so no such identity is asserted.
    end

    @info h2("MERS transmission_hc (THC, regular BIRTH): noop+cross ->" *
             " RegularFlow, fork -> decay bucket (cross-checked against" *
             " kli_reduce_test.jl's THC numbers)")
    @testset "MERS transmission_hc" begin
        thc = find_event(PhyloPOMP.MERS, :transmission_hc)
        @test thc.regular == true

        rts = reduced_transitions(thc, [2, 1], [5, 4])
        @test length(rts) == 3
        spec = classify_filter_terms(thc, rts)
        @test length(spec.regular) == 2
        @test length(spec.outflow_imbalance) == 1
        @test isempty(spec.singular)

        reg_by_kind = Dict(rf.reduced.kind => rf for rf in spec.regular)
        @test reg_by_kind[:noop].Φ == 3 // 5
        @test reg_by_kind[:cross].Φ == 3 // 20
        @test spec.outflow_imbalance[1].Φ == 1 // 20

        # Matches mers_naive.jl's regular_part! k==3 (noop, no cols
        # mutation) and k==4 (cross, `swap!(cols, Camel, Human, b)`) --
        # `fork!` never appears in regular_part! for this event.
    end

    @info h2("Provenance chain: FilterTerm -> ReducedTransition ->" *
             " KLITransition -> source Event, unbroken")
    @testset "provenance" begin
        infection = find_event(PhyloPOMP.SEIR, :infection)
        ts  = full_transitions(infection, [2, 2], [6, 5])
        rts = reduce_event_indicator(ts)
        spec = classify_filter_terms(infection, rts)

        # A RegularFlow term traces back to real KLITransitions of `event`.
        rf = only(rf for rf in spec.regular if rf.reduced.kind == :noop)
        @test rf.event === infection
        @test rf.reduced isa ReducedTransition
        @test !isempty(rf.reduced.transitions)
        @test all(t isa KLITransition for t in rf.reduced.transitions)
        @test all(t.event === infection for t in rf.reduced.transitions)
        @test all(t in ts for t in rf.reduced.transitions)  # traces to M03's own output
        @test rf.Φ == rf.reduced.Φ == sum(t.phi for t in rf.reduced.transitions)

        # A decay-bucket term traces back the same way.
        dt = only(spec.outflow_imbalance)
        @test dt.event === infection
        @test dt.reduced isa ReducedTransition
        @test dt.reduced.kind == :fork
        @test all(t.event === infection for t in dt.reduced.transitions)
        @test all(t in ts for t in dt.reduced.transitions)
        @test dt.Φ == dt.reduced.Φ == sum(t.phi for t in dt.reduced.transitions)

        # Mismatched-provenance guard: classifying under the wrong event
        # throws rather than silently mislabeling.
        progression = find_event(PhyloPOMP.SEIR, :progression)
        @test_throws ArgumentError classify_filter_terms(progression, rts)
    end

    @info h2("Singular BIRTH/MIGRATION branch: no such event currently" *
             " exists in SEIR/MERS (both models' only regular==false" *
             " events are SAMPLE-type, out of full_transitions's scope" *
             " per M03) -- validated here with a synthetic Event")
    @testset "synthetic singular BIRTH event" begin
        seir_infection = find_event(PhyloPOMP.SEIR, :infection)

        # Confirm the documented gap: every regular==false event in both
        # shipped models is SAMPLE-type.
        for model in (PhyloPOMP.SEIR, PhyloPOMP.MERS)
            for ev in model.events
                if !ev.regular
                    @test ev.type == SAMPLE
                end
            end
        end

        # Build a synthetic BIRTH event with regular=false, identical
        # structural fields to SEIR's `infection` otherwise, purely to
        # exercise the `!event.regular` branch of `classify_filter_terms`.
        synthetic = Event(:synthetic_singular_infection, seir_infection.Δ,
                           seir_infection.hazard, seir_infection.r, BIRTH,
                           seir_infection.from, seir_infection.into,
                           false, true)
        rts = reduced_transitions(synthetic, [2, 2], [6, 5])
        @test length(rts) == 3  # same (ell, n) as "SEIR infection" above

        spec = classify_filter_terms(synthetic, rts)
        @test isempty(spec.regular)
        @test isempty(spec.outflow_imbalance)
        @test length(spec.singular) == 3
        @test Set(sf.reduced.kind for sf in spec.singular) == Set([:noop, :cross, :fork])
        # Every reduced transition, including :fork, is available only as
        # SingularFlow -- consistent with mgp_filter.jl's regular_step!
        # already zeroing a non-regular event's contribution to the
        # continuous driver at every time other than its fixed singular
        # time (`alpha = [ev.regular ? kli_hazard(...) : 0.0 ...]`,
        # mgp_filter.jl:166-167), so there is no separate imbalance bucket
        # for a singular event's OWN outcomes.
    end

    @info h2("M06 Part 0 correction: OutflowImbalanceTerm's `mechanism` is" *
             " never `:lambda` -- this bucket is structurally distinct from" *
             " KLI's actual decay term (kli_decay, mgp_filter.jl:87-88," *
             " unverified stub, untouched by this milestone)")
    @testset "OutflowImbalanceTerm mechanism is not lambda (M06 correction)" begin
        infection = find_event(PhyloPOMP.SEIR, :infection)
        tcc       = find_event(PhyloPOMP.MERS, :transmission_cc)

        spec_seir = filter_spec(infection, [2, 2], [6, 5])
        spec_mers = filter_spec(tcc, [2, 0], [5, 3])

        # Every OutflowImbalanceTerm this pipeline can currently construct
        # carries mechanism = :inflow_outflow_imbalance -- never :lambda.
        # This is the type-level assertion that stands in for "no current
        # code path conflates this bucket with an actual λ computation":
        # since KLI's λ is not implemented anywhere yet (kli_decay is still
        # an `error(...)` stub), the only thing checkable today is that this
        # bucket's own tag never claims to *be* λ.
        for spec in (spec_seir, spec_mers)
            @test !isempty(spec.outflow_imbalance)
            for term in spec.outflow_imbalance
                @test term isa OutflowImbalanceTerm
                @test term.mechanism == :inflow_outflow_imbalance
                @test term.mechanism != :lambda
            end
        end

        # `FilterSpec` has no field literally named `decay` -- a future
        # milestone cannot accidentally `sum(spec.decay)` into a `lambda`
        # total, because that field does not exist (it must find
        # `outflow_imbalance` instead, and read why it is named that way).
        @test :decay ∉ fieldnames(FilterSpec)
        @test :outflow_imbalance in fieldnames(FilterSpec)

        # kli_decay (mgp_filter.jl's actual, still-unimplemented λ slot,
        # Eq. 47/B2) is a completely separate code path from
        # classify_filter_terms/OutflowImbalanceTerm -- it takes different
        # arguments (alpha, pi, cols, x, model) and is still an unverified
        # `error(...)` stub, confirming M05/M06 built no bridge between the
        # two.
        @test_throws ErrorException kli_decay(Float64[], Float64[], nothing,
                                               nothing, PhyloPOMP.SEIR)
    end

end

end # module
