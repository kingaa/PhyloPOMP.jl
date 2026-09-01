# mgp_audit.jl
# =============================================================================
# Structural audit / semantic validation of an MGPModel (Population IR).
#
#   Read-only inspection only: no likelihood math, no genealogy/coloring
#   state, no saturation/φ/Q_u computation (those are runtime, genealogy-
#   dependent quantities -- explicitly out of scope until M02+, see the
#   `Event` docstring in mgp.jl). This file answers one question: "can the
#   compiler mechanically inspect/audit the static Population IR that
#   already exists in `Event`/`MGPModel`?" -- and, as a stretch, "can it
#   also catch structural mistakes in that IR without throwing?"
#
#   `audit_model` is an ordinary function, not a macro, per this project's
#   rule that macros are syntax, not architecture. The `@mgpaudit` macro
#   this docstring originally forward-referenced is now built (M06 Part 2)
#   -- see `src/examples/mgp_mgpaudit.jl`, a SEPARATE file (not appended
#   here) because `mgpaudit`'s per-mark derivation walks M02-M06 IR
#   (`ReducedTransition`/`FilterSpec`/`explain`/...) that is only defined by
#   files included LATER than this one in `Examples.jl` -- Julia struct
#   field types must already exist at struct-definition time, so a struct
#   carrying e.g. `Vector{ReducedExplanation}` cannot itself live in a file
#   included before `mgp_explain.jl`. `mgp_audit.jl` stays exactly where M01
#   put it (its `audit_model`/`validate_model` genuinely have no such
#   dependency); `mgp_mgpaudit.jl` is included near the end of the chain,
#   after `mgp_explain.jl`, where its dependencies are already satisfied.
#
# STATUS: new in M01. Cross-checked against `SEIR_REFERENCE` (mgp.jl) and
#   against the `MERS` model (mgp_mers.jl) by `test/population_ir_test.jl`.
#
# Primary source: King, Lin & Ionides, "Exact phylodynamic likelihood via
#   structured Markov genealogy processes" (StructuredMGPs.pdf), for the
#   terminology only (Δ_u, α_u, r_u, W_u) -- no equation from the paper is
#   implemented in this file.
# =============================================================================

export audit_model, validate_model, EventAudit, ModelAudit

"""
    EventAudit

One row of a model audit: everything statically knowable about a single
`Event` without evaluating its hazard closure or touching genealogy state.

Fields mirror the Population IR terminology used throughout the compiler
project: `delta` = Δ_u (state jump), `has_hazard` = whether α_u is present
as a callable, `r` = r_u (production vector), `from`/`into` = W_u (deme
wiring), resolved here to both raw indices and the corresponding deme
*names* for readability.
"""
struct EventAudit
    name       :: Symbol
    type       :: EventType
    regular    :: Bool
    observed   :: Bool
    delta      :: Vector{Pair{Symbol,Int}}
    has_hazard :: Bool
    r          :: Vector{Int}
    from       :: Int
    from_deme  :: Union{Symbol,Nothing}
    into       :: Vector{Int}
    into_demes :: Vector{Symbol}
end

"""
    ModelAudit

Structural audit report for an `MGPModel`: its compartments, its
lineage-carrying demes, and one `EventAudit` per event. Returned by
`audit_model`. A plain, printable, structured value -- consume it
programmatically (e.g. from a future `@mgpaudit` macro) or just `show` it
at the REPL.
"""
struct ModelAudit
    name         :: Symbol
    compartments :: Vector{Symbol}
    demes        :: Vector{Symbol}
    events       :: Vector{EventAudit}
end

"""
    audit_model(model::MGPModel) -> ModelAudit

Inspect `model`'s Population IR and return a structured `ModelAudit`
exposing, for every event: its semantic type, regular/singular and
observed flags, state transition Δ_u, whether a hazard α_u is present,
production vector r_u, and source/target deme wiring W_u (`from`/`into`,
resolved to deme names as well as raw indices).

Purely inspective: never calls a hazard closure, never touches
genealogy/coloring state, never computes saturation/φ_u/Φ_u/Q_u (those are
runtime, genealogy-dependent quantities -- out of scope for this pass, see
`Event`'s docstring in mgp.jl). Does not throw on structural problems --
see `validate_model` for a checked, issue-reporting pass.

`show(io, MIME("text/plain"), audit)` renders a human-readable table.
"""
function audit_model(model::MGPModel)
    demes = model.demes
    events = EventAudit[]
    for ev in model.events
        from_deme = (1 <= ev.from <= length(demes)) ? demes[ev.from] : nothing
        into_demes = Symbol[(1 <= i <= length(demes)) ? demes[i] : Symbol("<invalid:$i>")
                             for i in ev.into]
        push!(events, EventAudit(
            ev.name, ev.type, ev.regular, ev.observed,
            ev.Δ, ev.hazard isa Function, ev.r,
            ev.from, from_deme, ev.into, into_demes,
        ))
    end
    ModelAudit(model.name, model.compartments, demes, events)
end

function Base.show(io::IO, ::MIME"text/plain", ev::EventAudit)
    tag = ev.regular ? "regular" : "singular"
    obs = ev.observed ? "observed" : "background"
    println(io, "  - ", ev.name, "  [", ev.type, ", ", tag, ", ", obs, "]")
    delta_str = isempty(ev.delta) ? "(none)" : join(("$(p.first)=$(p.second >= 0 ? "+" : "")$(p.second)" for p in ev.delta), ", ")
    println(io, "      Δ (state jump):     ", delta_str)
    println(io, "      α (hazard):         ", ev.has_hazard ? "present (closure)" : "MISSING")
    println(io, "      r (production):     ", ev.r)
    from_str = ev.from_deme === nothing ? "(none)" : "$(ev.from_deme) [deme #$(ev.from)]"
    println(io, "      W.from (source):    ", from_str)
    into_str = isempty(ev.into) ? "(none)" :
        join(("$(d) [deme #$(i)]" for (d, i) in zip(ev.into_demes, ev.into)), ", ")
    println(io, "      W.into (targets):   ", into_str)
end

Base.show(io::IO, ev::EventAudit) = show(io, MIME("text/plain"), ev)

function Base.show(io::IO, ::MIME"text/plain", a::ModelAudit)
    println(io, "ModelAudit(", a.name, ")")
    println(io, "  compartments: ", join(a.compartments, ", "))
    println(io, "  demes:        ", join(a.demes, ", "))
    println(io, "  events (", length(a.events), "):")
    for ev in a.events
        show(io, MIME("text/plain"), ev)
    end
end

Base.show(io::IO, a::ModelAudit) = show(io, MIME("text/plain"), a)

# =============================================================================
# validate_model -- separately-callable semantic validation.
#
#   M00 found the equivalent checks folded into @mgp's macroexpansion-time
#   error() calls (mgp_macro.jl:161,170,176,177-178,121-122) with no way to
#   run them against an already-constructed MGPModel (e.g. a hand-written
#   one like SEIR_REFERENCE). This closes that gap additively, without
#   touching the macro or its existing (working, tested) error() checks.
# =============================================================================

"""
    validate_model(model::MGPModel) -> Vector{String}

Run structural/semantic checks against `model` and return a list of
human-readable issue strings; an empty vector means the model passed every
check. Never throws.

Covers: duplicate compartment/deme/event names, every deme being a
declared compartment, every event's `from`/`into` deme indices being valid
indices into `model.demes`, every event's `r` vector length matching
`length(model.demes)`, every event's hazard being a callable `Function`,
every Δ entry referencing a declared compartment, and a few per-`EventType`
wiring sanity checks (e.g. a `BIRTH`/`MIGRATION`/`DEATH`/`SAMPLE` event
must have a source deme).

This is purely static-IR validation -- it does not (and cannot) check
anything dynamic/genealogy-dependent (saturation, φ_u, Φ_u, Q_u); that
machinery does not exist yet (see `docs/compiler/compiler_roadmap.md`
items 4-7).
"""
function validate_model(model::MGPModel)
    issues = String[]
    ndemes = length(model.demes)
    comp_set = Set(model.compartments)

    allunique(model.compartments) ||
        push!(issues, "duplicate compartment names: $(model.compartments)")
    allunique(model.demes) ||
        push!(issues, "duplicate deme names: $(model.demes)")
    for d in model.demes
        d in comp_set || push!(issues, "deme `$d` is not a declared compartment")
    end
    allunique(ev.name for ev in model.events) ||
        push!(issues, "duplicate event names among model.events")

    for ev in model.events
        tag = "event `$(ev.name)`"

        ev.hazard isa Function ||
            push!(issues, "$tag: hazard (α_u) is not a callable Function")

        length(ev.r) == ndemes ||
            push!(issues, "$tag: r (production vector) has length $(length(ev.r)), " *
                           "expected $(ndemes) (= number of demes)")

        (ev.from == 0 || (1 <= ev.from <= ndemes)) ||
            push!(issues, "$tag: from index $(ev.from) out of range 1:$(ndemes)")

        for i in ev.into
            (1 <= i <= ndemes) ||
                push!(issues, "$tag: into index $(i) out of range 1:$(ndemes)")
        end

        for p in ev.Δ
            p.first in comp_set ||
                push!(issues, "$tag: Δ references unknown compartment `$(p.first)`")
        end

        if ev.type in (BIRTH, MIGRATION, DEATH, SAMPLE) && ev.from == 0
            push!(issues, "$tag: $(ev.type) event has from=0 (no source deme)")
        end
        if ev.type == MIGRATION && isempty(ev.into)
            push!(issues, "$tag: MIGRATION event has no into deme")
        end
    end

    issues
end
