"""
M06 acceptance-gate tests: `explain` (`src/examples/mgp_explain.jl`, Part 1)
and `mgpaudit`/`@mgpaudit` (`src/examples/mgp_mgpaudit.jl`, Part 2).

`explain` walks the full provenance chain `FilterTerm -> ReducedTransition ->
KLITransition -> Event` built by M01-M05 (and M06 Part 0's
`DecayTerm -> OutflowImbalanceTerm` rename); `mgpaudit`/`@mgpaudit` walks
`audit_model` (M01) through `classify_filter_terms` (M05/M06) for every
event of a model and reports either a full derivation (BIRTH/MIGRATION) or
an explicit out-of-scope note (DEATH/SAMPLE/NEUTRAL).

Every concrete `(ℓ, n)` instance used below for `explain`'s own tests is
REUSED, unchanged, from `test/kli_reduce_test.jl` / `test/kli_filter_ir_test.jl`
(SEIR `infection`: `[2, 2]`/`[6, 5]`; MERS `transmission_cc`: `[2, 0]`/`[5, 3]`)
-- per the milestone's explicit instruction not to invent new numbers for
this check. `mgpaudit`'s own DEFAULT `(ℓ, n)` state (`default_audit_state`)
is a separate, independently documented choice (see that function's
docstring in `mgp_mgpaudit.jl`) -- not required to match these test values.
"""
module MgpAuditTest

import ..Main: h1, h2

@info h1("Provenance / explain and mgpaudit / @mgpaudit (M06)")

using Test
using PhyloPOMP
using PhyloPOMP: filter_spec, RegularFlow, SingularFlow, OutflowImbalanceTerm,
    FilterSpec, ReducedTransition, KLITransition,
    explain, TransitionExplanation, ReducedExplanation, TermExplanation,
    mgpaudit, @mgpaudit, default_audit_state,
    MarkAudit, MarkAuditInScope, MarkAuditOutOfScope, MgpAuditReport,
    Event, EventType, BIRTH, MIGRATION, DEATH, SAMPLE, NEUTRAL

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

const SEIR_EVENT_NAMES = [:infection, :progression, :recovery, :waning, :sampling]
const MERS_EVENT_NAMES = [:transmission_cc, :transmission_hh, :transmission_hc,
                           :transmission_ch, :removal_c, :removal_h,
                           :sampling_c, :sampling_h, :birth_c, :birth_h,
                           :death_c, :death_h]

@testset verbose=true "Provenance / explain and mgpaudit (M06)" begin

    @info h2("explain(::KLITransition) / explain(::ReducedTransition) --" *
             " bottom-of-chain provenance, SEIR infection instance")
    @testset "explain: KLITransition / ReducedTransition" begin
        infection = find_event(PhyloPOMP.SEIR, :infection)
        # Same instance as kli_reduce_test.jl / kli_filter_ir_test.jl
        # "SEIR infection": Φ(noop)=8/15, Φ(cross)=1/10, Φ(fork)=1/30.
        spec = filter_spec(infection, [2, 2], [6, 5])

        rf = only(rf for rf in spec.regular if rf.reduced.kind == :noop)
        te = explain(rf.reduced)
        @test te isa ReducedExplanation
        @test te.event_name == :infection
        @test te.event_type == BIRTH
        @test te.kind == :noop
        @test te.Φ == 8 // 15
        @test length(te.members) == 2  # Identity + InlineSameDeme collapsed
        @test Set(m.kind for m in te.members) == Set([:Identity, :InlineSameDeme])
        for m in te.members
            @test m.event_name == :infection
            @test m.event_type == BIRTH
        end

        # Every individual KLITransition explains correctly too.
        for t in rf.reduced.transitions
            et = explain(t)
            @test et isa TransitionExplanation
            @test et.event_name == :infection
            @test et.event_type == BIRTH
            @test et.r == [1, 1]
            @test et.phi == t.phi
            @test et.s == t.s
        end
    end

    @info h2("explain(::RegularFlow) correctly identifies event name and" *
             " transition kind -- SEIR infection")
    @testset "explain: RegularFlow (SEIR infection)" begin
        infection = find_event(PhyloPOMP.SEIR, :infection)
        spec = filter_spec(infection, [2, 2], [6, 5])
        rf = only(rf for rf in spec.regular if rf.reduced.kind == :cross)

        te = explain(rf)
        @test te isa TermExplanation
        @test te.bucket == :regular_flow
        @test te.event_name == :infection
        @test te.event_type == BIRTH
        @test te.r == [1, 1]
        @test te.Φ == 1 // 10
        @test te.reduced.kind == :cross
        @test te.mechanism === nothing
        @test te.reason === nothing

        # Human-readable output actually contains the identifying info.
        str = sprint(show, MIME("text/plain"), te)
        @test occursin("infection", str)
        @test occursin("BIRTH", str)
        @test occursin("regular_flow", str) || occursin("RegularFlow", str)
    end

    @info h2("explain(::OutflowImbalanceTerm) correctly identifies event" *
             " name, transition kind, and the M06 Part-0 mechanism tag --" *
             " MERS transmission_cc")
    @testset "explain: OutflowImbalanceTerm (MERS transmission_cc)" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        # Same instance as kli_reduce_test.jl "TCC" / kli_filter_ir_test.jl
        # "MERS transmission_cc": Φ(noop)=3/5, Φ(fork)=1/10.
        spec = filter_spec(tcc, [2, 0], [5, 3])
        oi = only(spec.outflow_imbalance)

        te = explain(oi)
        @test te isa TermExplanation
        @test te.bucket == :outflow_imbalance
        @test te.event_name == :transmission_cc
        @test te.event_type == BIRTH
        @test te.r == [2, 0]
        @test te.Φ == 1 // 10
        @test te.reduced.kind == :fork
        @test te.mechanism == :inflow_outflow_imbalance
        @test te.reason == :fork_unobserved_at_regular_time
        # Structural distinctness from lambda re-asserted here too (not just
        # kli_filter_ir_test.jl's dedicated testset): explain() never labels
        # this mechanism :lambda.
        @test te.mechanism != :lambda

        str = sprint(show, MIME("text/plain"), te)
        @test occursin("transmission_cc", str)
        @test occursin("BIRTH", str)
        @test occursin("inflow_outflow_imbalance", str)
        @test occursin("NOT KLI's λ", str) || occursin("NOT λ", str)
    end

    @info h2("mgpaudit(SEIR) / @mgpaudit SEIR run without error and cover" *
             " every one of SEIR's 5 events by name")
    @testset "mgpaudit(SEIR) completeness" begin
        report = mgpaudit(PhyloPOMP.SEIR)
        @test report isa MgpAuditReport
        @test report.model_name == :SEIR
        @test length(report.marks) == length(SEIR_EVENT_NAMES)

        seen = Symbol[m.event_audit.name for m in report.marks]
        @test Set(seen) == Set(SEIR_EVENT_NAMES)

        # BIRTH/MIGRATION -> in scope; DEATH/SAMPLE/NEUTRAL -> explicit
        # out-of-scope note (never silently omitted).
        for m in report.marks
            if m.event_audit.type in (BIRTH, MIGRATION)
                @test m isa MarkAuditInScope
                @test !isempty(m.reduced)
                @test !isempty(m.terms)
            else
                @test m isa MarkAuditOutOfScope
                @test occursin("out of scope", m.note)
            end
        end

        # infection (BIRTH) must show a RegularFlow and an
        # OutflowImbalanceTerm at this report's default state.
        infection_mark = only(m for m in report.marks if m.event_audit.name == :infection)
        @test any(t.bucket == :regular_flow for t in infection_mark.terms)
        @test any(t.bucket == :outflow_imbalance for t in infection_mark.terms)

        # progression (MIGRATION) structurally never has an
        # OutflowImbalanceTerm (M04's finding, carried through M05/M06).
        progression_mark = only(m for m in report.marks if m.event_audit.name == :progression)
        @test all(t.bucket != :outflow_imbalance for t in progression_mark.terms)

        # The printed report also mentions every event name -- a
        # human-readable completeness check, not just the struct's.
        str = sprint(show, MIME("text/plain"), report)
        for name in SEIR_EVENT_NAMES
            @test occursin(String(name), str)
        end

        # @mgpaudit is a thin macro forwarding to mgpaudit -- same result.
        report_macro = @mgpaudit SEIR
        @test report_macro isa MgpAuditReport
        @test Set(m.event_audit.name for m in report_macro.marks) == Set(SEIR_EVENT_NAMES)
    end

    @info h2("mgpaudit(MERS) / @mgpaudit MERS run without error and cover" *
             " every one of MERS's 12 events by name")
    @testset "mgpaudit(MERS) completeness" begin
        report = mgpaudit(PhyloPOMP.MERS)
        @test report isa MgpAuditReport
        @test report.model_name == :MERS
        @test length(report.marks) == length(MERS_EVENT_NAMES)

        seen = Symbol[m.event_audit.name for m in report.marks]
        @test Set(seen) == Set(MERS_EVENT_NAMES)

        for m in report.marks
            if m.event_audit.type in (BIRTH, MIGRATION)
                @test m isa MarkAuditInScope
                @test !isempty(m.reduced)
                @test !isempty(m.terms)
            else
                @test m isa MarkAuditOutOfScope
                @test occursin("out of scope", m.note)
            end
        end

        # All four transmission marks are BIRTH and, at the report's default
        # state, each shows a RegularFlow and an OutflowImbalanceTerm.
        for name in (:transmission_cc, :transmission_hh, :transmission_hc, :transmission_ch)
            mark = only(m for m in report.marks if m.event_audit.name == name)
            @test mark isa MarkAuditInScope
            @test any(t.bucket == :regular_flow for t in mark.terms)
            @test any(t.bucket == :outflow_imbalance for t in mark.terms)
        end

        # removal_c/h are DEATH, sampling_c/h are SAMPLE, birth_c/h and
        # death_c/h (demography) are NEUTRAL -- all explicitly out of scope.
        for name in (:removal_c, :removal_h)
            mark = only(m for m in report.marks if m.event_audit.name == name)
            @test mark isa MarkAuditOutOfScope
            @test mark.event_audit.type == DEATH
        end
        for name in (:sampling_c, :sampling_h)
            mark = only(m for m in report.marks if m.event_audit.name == name)
            @test mark isa MarkAuditOutOfScope
            @test mark.event_audit.type == SAMPLE
        end
        for name in (:birth_c, :birth_h, :death_c, :death_h)
            mark = only(m for m in report.marks if m.event_audit.name == name)
            @test mark isa MarkAuditOutOfScope
            @test mark.event_audit.type == NEUTRAL
        end

        str = sprint(show, MIME("text/plain"), report)
        for name in MERS_EVENT_NAMES
            @test occursin(String(name), str)
        end

        report_macro = @mgpaudit MERS
        @test report_macro isa MgpAuditReport
        @test Set(m.event_audit.name for m in report_macro.marks) == Set(MERS_EVENT_NAMES)
    end

    @info h2("@mgpaudit forwards explicit ℓ/n keyword arguments (thin" *
             " macro, no derivation logic of its own)")
    @testset "@mgpaudit with explicit state" begin
        # Reuse the exact validated MERS THC instance from
        # kli_reduce_test.jl / kli_filter_ir_test.jl.
        report = @mgpaudit MERS ℓ=[2, 1] n=[5, 4]
        @test report.ℓ == [2, 1]
        @test report.n == [5, 4]
        thc = only(m for m in report.marks if m.event_audit.name == :transmission_hc)
        @test thc isa MarkAuditInScope
        reg_by_kind = Dict(rt.kind => rt for rt in thc.reduced)
        @test reg_by_kind[:noop].Φ == 3 // 5
        @test reg_by_kind[:cross].Φ == 3 // 20
        @test reg_by_kind[:fork].Φ == 1 // 20
    end

    @info h2("mgpaudit validates ℓ/n shape and the ℓ<=n invariant")
    @testset "mgpaudit argument validation" begin
        @test_throws ArgumentError mgpaudit(PhyloPOMP.SEIR; ℓ = [1, 1, 1])
        @test_throws ArgumentError mgpaudit(PhyloPOMP.SEIR; ℓ = [5, 5], n = [1, 1])
    end

    @info h2("default_audit_state matches this milestone's documented" *
             " per-model choices")
    @testset "default_audit_state" begin
        @test default_audit_state(PhyloPOMP.SEIR) == ([2, 2], [6, 5])
        @test default_audit_state(PhyloPOMP.MERS) == ([2, 2], [5, 5])
    end

end

end # module
