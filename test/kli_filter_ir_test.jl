"""
Tests for classify_filter_terms and filter_spec: reduced transitions go into
RegularFlow, SingularFlow and OutflowImbalanceTerm buckets. (ell, n) instances
and Phi values match kli_reduce_test.jl.
"""
module KliFilterIrTest

import ..Main: h1, h2

@info h1("Filter IR: classify_filter_terms / filter_spec")

using Test
using PhyloPOMP
using PhyloPOMP: full_transitions, reduce_event_indicator, reduced_transitions,
    ReducedTransition, KLITransition,
    FilterTerm, RegularFlow, SingularFlow, OutflowImbalanceTerm, FilterSpec,
    classify_filter_terms, filter_spec,
    Event, EventType, BIRTH, MIGRATION, SAMPLE, kli_decay

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

@testset verbose=true "Filter IR classification" begin

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
        @test isempty(spec.outflow_imbalance)  # no fork term, so no decay-bucket term

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

        # Phi(noop)+Phi(fork) = 7/10, not 1: the sum of phi over s is not a probability.
    end

    @info h2("MERS transmission_hc (THC): noop+cross to RegularFlow, fork to the decay bucket")
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
    end

    @info h2("Trace-back chain: FilterTerm -> ReducedTransition ->" *
             " KLITransition -> source Event, unbroken")
    @testset "trace back to source event" begin
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
        @test all(t in ts for t in rf.reduced.transitions)  # traces to full_transitions' output
        @test rf.Φ == rf.reduced.Φ == sum(t.phi for t in rf.reduced.transitions)

        # A decay-bucket term traces back the same way.
        dt = only(spec.outflow_imbalance)
        @test dt.event === infection
        @test dt.reduced isa ReducedTransition
        @test dt.reduced.kind == :fork
        @test all(t.event === infection for t in dt.reduced.transitions)
        @test all(t in ts for t in dt.reduced.transitions)
        @test dt.Φ == dt.reduced.Φ == sum(t.phi for t in dt.reduced.transitions)

        # Wrong-event guard: classifying under the wrong event
        # throws rather than silently mislabeling.
        progression = find_event(PhyloPOMP.SEIR, :progression)
        @test_throws ArgumentError classify_filter_terms(progression, rts)
    end

    @info h2("Singular BIRTH/MIGRATION branch, via a synthetic Event (no shipped model has one)")
    @testset "synthetic singular BIRTH event" begin
        seir_infection = find_event(PhyloPOMP.SEIR, :infection)

        # Every regular==false event in both
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
        # A non-regular event has no imbalance bucket; all its outcomes are singular.
    end

    @info h2("OutflowImbalanceTerm.mechanism is never :lambda")
    @testset "OutflowImbalanceTerm mechanism is not lambda" begin
        infection = find_event(PhyloPOMP.SEIR, :infection)
        tcc       = find_event(PhyloPOMP.MERS, :transmission_cc)

        spec_seir = filter_spec(infection, [2, 2], [6, 5])
        spec_mers = filter_spec(tcc, [2, 0], [5, 3])

        # The mechanism tag never claims to be KLI's lambda; lambda is not implemented.
        for spec in (spec_seir, spec_mers)
            @test !isempty(spec.outflow_imbalance)
            for term in spec.outflow_imbalance
                @test term isa OutflowImbalanceTerm
                @test term.mechanism == :inflow_outflow_imbalance
                @test term.mechanism != :lambda
            end
        end

        # FilterSpec has no field named decay.
        @test :decay ∉ fieldnames(FilterSpec)
        @test :outflow_imbalance in fieldnames(FilterSpec)

        # kli_decay is a separate code path (mgp_filter.jl); at a state it equals compiled_decay.
        cols = PhyloPOMP.Coloring(PhyloPOMP.NaiveSEIR.Demes)
        push!(cols[PhyloPOMP.NaiveSEIR.Infec], 1)
        x = (S = 30, E = 4, I = 6, R = 5)
        θ = (β = 3.0, σ = 1.0, γ = 1.0, ω = 0.5, ψ = 0.2, χ = 0.1, N = 100.0)
        slots = PhyloPOMP.kli_slots(PhyloPOMP.SEIR)
        al, pv = zeros(length(slots)), zeros(length(slots))
        PhyloPOMP.kli_rates!(al, pv, slots, cols, x, θ, PhyloPOMP.SEIR)
        @test kli_decay(al, pv, slots, cols, x, θ, PhyloPOMP.SEIR) ≈
              PhyloPOMP.compiled_decay(PhyloPOMP.SEIR, x, θ, [0, 1], [4, 6])
    end

end

end # module
