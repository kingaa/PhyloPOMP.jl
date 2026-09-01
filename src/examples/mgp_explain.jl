# mgp_explain.jl
# =============================================================================
# M06 Part 1: Provenance / `explain`.
#
#   `explain(x)` walks the full derivation chain M01-M05 built --
#
#       FilterTerm -> ReducedTransition -> KLITransition -> Event
#
#   -- and returns a structured, printable explanation of `x`, covering (per
#   the milestone spec): the mark/event name and its semantic `EventType`;
#   the production vector `r_u`; the saturation `s` of each contributing full
#   transition; the full KLI transition kind (`Identity`/`InlineSameDeme`/
#   `CrossDeme`/`Fork`) and its `φ_u`; which full transitions collapsed into
#   a given `ReducedTransition` and its `Φ_u`; and, for a `FilterTerm`, the
#   `RegularFlow`/`SingularFlow`/`OutflowImbalanceTerm` classification
#   (renamed from `DecayTerm`, see `mgp_filter_ir.jl`'s header and
#   `handoffs/M06_audit_provenance.md`'s "Correction to M05" section).
#
#   `explain` is defined for every `FilterTerm` concrete subtype
#   (`RegularFlow`, `SingularFlow`, `OutflowImbalanceTerm`), and, since the
#   provenance chain is available one level down too, also for a bare
#   `ReducedTransition` (M04) and a bare `KLITransition` (M03) -- so a future
#   milestone can call `explain` on ANY node of the chain, not just the top.
#
#   Follows M01's `audit_model`/`ModelAudit` convention (`mgp_audit.jl`):
#   an ordinary function returns a structured, immutable value; a
#   `Base.show(io, MIME"text/plain", ...)` method renders it for humans.
#   `explain` is read-only inspection over already-computed IR -- it
#   performs no new KLI math, gates nothing, and computes no decay/lambda
#   weight (still M07's job).
#
# Primary sources: this repo's own M01-M05 IR (`mgp_audit.jl`, `mgp_phi.jl`,
#   `mgp_transitions.jl`, `mgp_reduce.jl`, `mgp_filter_ir.jl`), consumed
#   read-only; `mers_filter_suite.tex` line citations already established by
#   those files (repeated in the docstrings below only where load-bearing).
# =============================================================================

export explain, TransitionExplanation, ReducedExplanation, TermExplanation

# -----------------------------------------------------------------------------
# Per-KLITransition-subtype "kind" tag and human-readable detail string.
# Mirrors the vocabulary already established by mgp_transitions.jl (M03) and
# mgp_reduce.jl (M04) -- not a new taxonomy, just a presentation layer over
# the existing dispatch.
# -----------------------------------------------------------------------------
kli_kind_symbol(::IdentityTransition)        = :Identity
kli_kind_symbol(::InlineSameDemeTransition)  = :InlineSameDeme
kli_kind_symbol(::CrossDemeTransition)       = :CrossDeme
kli_kind_symbol(::ForkTransition)            = :Fork

kli_detail_string(t::IdentityTransition) =
    "no production slot occupied (y'=y)"
kli_detail_string(t::InlineSameDemeTransition) =
    "1 slot occupied, ancestral deme == post-event deme (deme #$(t.deme)); " *
    "reduced no-op (y'=y), but a distinct full transition (m changes)"
kli_detail_string(t::CrossDemeTransition) =
    "1 slot occupied, moves deme #$(t.ancestral_deme) -> deme #$(t.target_deme)"
kli_detail_string(t::ForkTransition) =
    "$(length(t.slot_demes)) slots occupied from ancestral deme #$(t.ancestral_deme), " *
    "landing in demes $(t.slot_demes) (branch point)"

# -----------------------------------------------------------------------------
# TransitionExplanation -- explain(t::KLITransition)
# -----------------------------------------------------------------------------

"""
    TransitionExplanation

Structured explanation of a single `KLITransition` (M03): the source mark
(`event_name`/`event_type`/`r`), the saturation `s` that produced it, its
exact compatibility factor `phi` (`φ_u(s)`, `Rational{Int}`), and which of
the four KLI-meaningful full-transition kinds it is (`kind` ∈
`(:Identity, :InlineSameDeme, :CrossDeme, :Fork)`), plus a human-readable
`detail` string.

Returned by `explain(t::KLITransition)`.
"""
struct TransitionExplanation
    event_name :: Symbol
    event_type :: EventType
    r          :: Vector{Int}
    kind       :: Symbol
    s          :: Vector{Int}
    phi        :: Rational{Int}
    detail     :: String
end

"""
    explain(t::KLITransition) -> TransitionExplanation

Explain a single full (uncollapsed) KLI-compatible coloring transition: which
mark (`event`) produced it, its saturation `s`, its exact compatibility
factor `φ_u(s)` (`t.phi`), and which of `IdentityTransition`/
`InlineSameDemeTransition`/`CrossDemeTransition`/`ForkTransition` it is.
Bottom of the M06 provenance chain -- `t.event` is the ultimate source
`Event`, already directly reachable, so `explain` needs no further descent.
"""
function explain(t::KLITransition)
    TransitionExplanation(t.event.name, t.event.type, production_slots(t.event),
                           kli_kind_symbol(t), t.s, t.phi, kli_detail_string(t))
end

function Base.show(io::IO, ::MIME"text/plain", e::TransitionExplanation)
    println(io, "TransitionExplanation for mark `", e.event_name, "` (",
                 e.event_type, ", r=", e.r, ")")
    println(io, "  full-transition kind: ", e.kind)
    println(io, "  saturation s:         ", e.s)
    println(io, "  φ_u(s):               ", e.phi)
    println(io, "  detail:               ", e.detail)
end
Base.show(io::IO, e::TransitionExplanation) = show(io, MIME("text/plain"), e)

# -----------------------------------------------------------------------------
# ReducedExplanation -- explain(rt::ReducedTransition)
# -----------------------------------------------------------------------------

"""
    ReducedExplanation

Structured explanation of a `ReducedTransition` (M04): the source mark, the
reduced `key`/`kind` (`:noop`/`:cross`/`:fork`), the summed `Φ` (`Φ_u`), and
the `TransitionExplanation` of every full `KLITransition` that collapsed
into it (M03's `phi_u` values that were summed to produce `Φ`).

Returned by `explain(rt::ReducedTransition)`.
"""
struct ReducedExplanation
    event_name :: Symbol
    event_type :: EventType
    key        :: Tuple
    kind       :: Symbol
    Φ          :: Rational{Int}
    members    :: Vector{TransitionExplanation}
end

"""
    explain(rt::ReducedTransition) -> ReducedExplanation

Explain a `ReducedTransition`: which reduced outcome it represents
(`key`/`kind`), its summed compatibility factor `Φ_u` (`rt.Φ`), and which
full `KLITransition`s (each explained via `explain(::KLITransition)`)
collapsed into it, per M04's `reduce_event_indicator` grouping
(`mers_filter_suite.tex`'s `Phi_u = sum phi_u` reduction, lines 423-432).
`rt.transitions` is assumed non-empty (M04's `reduce_event_indicator` never
constructs an empty group; throws `ArgumentError` if violated, since an empty
group has no `event` to report).
"""
function explain(rt::ReducedTransition)
    isempty(rt.transitions) &&
        throw(ArgumentError("explain(::ReducedTransition): rt.transitions is empty -- " *
                             "cannot determine the source event (this should not happen; " *
                             "reduce_event_indicator never constructs an empty group)."))
    ev = rt.transitions[1].event
    members = [explain(t) for t in rt.transitions]
    ReducedExplanation(ev.name, ev.type, rt.key, rt.kind, rt.Φ, members)
end

function Base.show(io::IO, ::MIME"text/plain", e::ReducedExplanation)
    println(io, "ReducedExplanation for mark `", e.event_name, "` (", e.event_type, ")")
    println(io, "  reduced kind: ", e.kind, "   key: ", e.key)
    println(io, "  Φ_u (summed): ", e.Φ)
    println(io, "  collapsed from ", length(e.members), " full transition(s):")
    for m in e.members
        println(io, "    - ", m.kind, "  s=", m.s, "  φ_u=", m.phi, "  (", m.detail, ")")
    end
end
Base.show(io::IO, e::ReducedExplanation) = show(io, MIME("text/plain"), e)

# -----------------------------------------------------------------------------
# TermExplanation -- explain(term::FilterTerm)
# -----------------------------------------------------------------------------

"""
    TermExplanation

Structured explanation of a `FilterTerm` (M05/M06's `RegularFlow`,
`SingularFlow`, or `OutflowImbalanceTerm`): the classification `bucket`, the
source mark, and the full `ReducedExplanation` (which, in turn, carries every
contributing `TransitionExplanation`) it was built from.

`mechanism`/`reason` are only populated (non-`nothing`) for
`OutflowImbalanceTerm` -- see that type's docstring in `mgp_filter_ir.jl` for
why this bucket is NOT KLI's `lambda` decay term, and
`handoffs/M06_audit_provenance.md`'s "Correction to M05" section for the
rename this milestone made (`DecayTerm` -> `OutflowImbalanceTerm`) precisely
so this distinction cannot be silently lost.

Returned by `explain(term::FilterTerm)`.
"""
struct TermExplanation
    bucket     :: Symbol   # :regular_flow, :singular_flow, :outflow_imbalance
    event_name :: Symbol
    event_type :: EventType
    r          :: Vector{Int}
    Φ          :: Rational{Int}
    reduced    :: ReducedExplanation
    mechanism  :: Union{Symbol,Nothing}
    reason     :: Union{Symbol,Nothing}
end

"""
    explain(term::RegularFlow) -> TermExplanation
    explain(term::SingularFlow) -> TermExplanation
    explain(term::OutflowImbalanceTerm) -> TermExplanation

Explain a classified `FilterTerm`: which bucket it landed in
(`:regular_flow`/`:singular_flow`/`:outflow_imbalance`), the source mark
(`term.event`), its `Φ`, and the full provenance of the `ReducedTransition`
(and, transitively, every collapsed `KLITransition`) it was built from
(`explain(term.reduced)`). For `OutflowImbalanceTerm`, also reports
`mechanism` (always `:inflow_outflow_imbalance` in the current milestone --
see `OutflowImbalanceTerm`'s docstring for why this is NOT `:lambda`) and
`reason` (`:fork_unobserved_at_regular_time`).

This is the top of the M06 provenance chain the milestone spec asks for:
`FilterTerm -> ReducedTransition -> KLITransition -> Event`, fully walked and
printable via `Base.show`.
"""
explain(term::RegularFlow) =
    TermExplanation(:regular_flow, term.event.name, term.event.type,
                     production_slots(term.event), term.Φ, explain(term.reduced),
                     nothing, nothing)

explain(term::SingularFlow) =
    TermExplanation(:singular_flow, term.event.name, term.event.type,
                     production_slots(term.event), term.Φ, explain(term.reduced),
                     nothing, nothing)

explain(term::OutflowImbalanceTerm) =
    TermExplanation(:outflow_imbalance, term.event.name, term.event.type,
                     production_slots(term.event), term.Φ, explain(term.reduced),
                     term.mechanism, term.reason)

const BUCKET_LABEL = Dict(
    :regular_flow      => "RegularFlow (continuous-time regular-event filter, KLI Eq. 45)",
    :singular_flow     => "SingularFlow (fixed observed data-event time only, KLI \"β^ev\")",
    :outflow_imbalance => "OutflowImbalanceTerm (fork/branch-point mass resolved via the " *
                           "regular inflow/outflow balance -- NOT KLI's λ; see mgp_filter_ir.jl)",
)

function Base.show(io::IO, ::MIME"text/plain", e::TermExplanation)
    println(io, "TermExplanation for mark `", e.event_name, "` (", e.event_type, ", r=", e.r, ")")
    println(io, "  bucket: ", get(BUCKET_LABEL, e.bucket, String(e.bucket)))
    println(io, "  Φ:      ", e.Φ)
    if e.bucket == :outflow_imbalance
        println(io, "  mechanism: ", e.mechanism,
                     e.mechanism == :inflow_outflow_imbalance ?
                         " (NOT KLI's λ -- see mers_filter_suite.tex lines 760-762, 841-844)" :
                         "")
        println(io, "  reason:    ", e.reason)
    end
    println(io, "  --- derivation (ReducedTransition -> KLITransition -> Event) ---")
    r = e.reduced
    println(io, "  reduced kind: ", r.kind, "   key: ", r.key, "   Φ_u: ", r.Φ)
    for m in r.members
        println(io, "    - ", m.kind, "  s=", m.s, "  φ_u=", m.phi, "  (", m.detail, ")")
    end
end
Base.show(io::IO, e::TermExplanation) = show(io, MIME("text/plain"), e)
