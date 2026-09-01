"""
Structural completeness check of the Population IR (`Event`/`MGPModel`)
for both models currently expressible via `@mgp`: `SEIR` and `MERS`.

Confirms that for every event the compiler can recover/compute Δ_u
(state jump), α_u (hazard), r_u (production vector), and W_u (from/into
deme wiring) -- the four static per-mark quantities the Population IR is
supposed to carry per `docs/compiler/compiler_roadmap.md` item 2. τ_u
(event time) is correctly *not* checked here: it is a runtime quantity,
not part of the static `Event` description (see `Event`'s docstring in
mgp.jl).

Also exercises `audit_model`/`validate_model` (`src/examples/mgp_audit.jl`,
new in M01) against both models: `audit_model(SEIR)` is cross-checked
against the hand-written `SEIR_REFERENCE` table, and both models are
checked to `validate_model` clean (empty issue list).
"""
module PopulationIRTest

import ..Main: h1, h2

@info h1("Population IR structural audit (SEIR, MERS)")

using Test
using PhyloPOMP

@testset verbose=true "Population IR structural audit" begin

    @info h2("every event exposes Δ, α, r, W (SEIR + MERS)")
    for model in (PhyloPOMP.SEIR, PhyloPOMP.MERS)
        ndemes = length(model.demes)
        @test !isempty(model.events)
        for ev in model.events
            # Δ_u: state jump, every referenced compartment declared.
            @test ev.Δ isa Vector{Pair{Symbol,Int}}
            @test all(p -> p.first in model.compartments, ev.Δ)

            # α_u: hazard, a callable closure (x, θ) -> Float64.
            @test ev.hazard isa Function

            # r_u: production vector, one entry per lineage-carrying deme.
            @test ev.r isa Vector{Int}
            @test length(ev.r) == ndemes
            @test all(>=(0), ev.r)

            # W_u: source/target deme wiring, valid indices into model.demes
            # (or the documented "none" sentinels: from=0, into=[]).
            @test ev.from == 0 || (1 <= ev.from <= ndemes)
            @test all(i -> 1 <= i <= ndemes, ev.into)

            # τ_u is deliberately NOT stored on Event (runtime quantity) --
            # confirm no field pretends to be it.
            @test !hasproperty(ev, :τ) && !hasproperty(ev, :t) && !hasproperty(ev, :time)

            # Semantic type / regular-singular / observed metadata present.
            @test ev.type isa PhyloPOMP.EventType
            @test ev.regular isa Bool
            @test ev.observed isa Bool
        end
    end

    @info h2("hazards evaluate to finite, non-negative reals at a real state")
    x_seir = (S = 90, E = 4, I = 5, R = 1)
    θ_seir = (β = 4.0, σ = 0.8, γ = 0.5, ω = 0.2, ψ = 0.03, χ = 0.0, N = 100.0)
    for ev in PhyloPOMP.SEIR.events
        v = ev.hazard(x_seir, θ_seir)
        @test v isa Real
        @test isfinite(v) && v >= 0
    end

    x_mers = (S_c = 19, I_c = 1, S_h = 20, I_h = 0)
    θ_mers = (β_cc = 3.0, β_ch = 0.5, β_hc = 0.5, β_hh = 3.0,
              γ_c = 1.0, γ_h = 1.0, χ_c = 0.3, χ_h = 0.3,
              B_c = 0.0, B_h = 0.0, N_c = 20.0, N_h = 20.0)
    for ev in PhyloPOMP.MERS.events
        v = ev.hazard(x_mers, θ_mers)
        @test v isa Real
        @test isfinite(v) && v >= 0
    end

    @info h2("audit_model(SEIR) matches SEIR_REFERENCE")
    audit = PhyloPOMP.audit_model(PhyloPOMP.SEIR)
    ref = PhyloPOMP.SEIR_REFERENCE
    @test audit.name == ref.name
    @test audit.compartments == ref.compartments
    @test audit.demes == ref.demes
    @test length(audit.events) == length(ref.events)
    for (a, ev) in zip(audit.events, ref.events)
        @test a.name == ev.name
        @test a.type == ev.type
        @test a.regular == ev.regular
        @test a.observed == ev.observed
        @test a.delta == ev.Δ
        @test a.has_hazard == (ev.hazard isa Function)
        @test a.r == ev.r
        @test a.from == ev.from
        @test a.into == ev.into
        @test a.from == 0 || a.from_deme == ref.demes[ev.from]
        @test a.into_demes == [ref.demes[i] for i in ev.into]
    end

    @info h2("audit_model(MERS) structural invariants")
    audit_mers = PhyloPOMP.audit_model(PhyloPOMP.MERS)
    @test audit_mers.name == :MERS
    @test length(audit_mers.events) == length(PhyloPOMP.MERS.events)
    for a in audit_mers.events
        @test a.from == 0 || 1 <= a.from <= length(PhyloPOMP.MERS.demes)
        @test all(i -> 1 <= i <= length(PhyloPOMP.MERS.demes),
                  [findfirst(==(d), PhyloPOMP.MERS.demes) for d in a.into_demes])
        @test length(a.r) == length(PhyloPOMP.MERS.demes)
    end

    @info h2("validate_model reports no issues for SEIR or MERS")
    @test PhyloPOMP.validate_model(PhyloPOMP.SEIR) == String[]
    @test PhyloPOMP.validate_model(PhyloPOMP.MERS) == String[]
    @test PhyloPOMP.validate_model(PhyloPOMP.SEIR_REFERENCE) == String[]

    @info h2("validate_model catches a deliberately broken model")
    broken = PhyloPOMP.MGPModel(:Broken, [:A, :B], [:A],
        [PhyloPOMP.Event(:bad, [:A=>-1, :ZZZ=>+1], (x,θ)->1.0,
                          [1, 2], PhyloPOMP.BIRTH, 5, [7], true, false)])
    issues = PhyloPOMP.validate_model(broken)
    @test !isempty(issues)
    @test any(occursin("r (production vector) has length", i) for i in issues)
    @test any(occursin("from index", i) for i in issues)
    @test any(occursin("into index", i) for i in issues)
    @test any(occursin("unknown compartment", i) for i in issues)

    @info h2("audit_model / EventAudit / ModelAudit show output is non-empty")
    io = IOBuffer()
    show(io, MIME("text/plain"), audit)
    s = String(take!(io))
    @test occursin("ModelAudit(SEIR)", s)
    @test occursin("infection", s)
    @test occursin("sampling", s)

end

end
