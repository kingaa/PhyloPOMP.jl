# explain(x) returns a printable record of where a FilterTerm,
# ReducedTransition or KLITransition came from, following the chain
# FilterTerm -> ReducedTransition -> KLITransition -> Event. Read-only.

export explain, TransitionExplanation, ReducedExplanation, TermExplanation

# Kind tag and detail string per KLITransition subtype.
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

"""
    TransitionExplanation

Structured explanation of a single `KLITransition`: the source mark
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
`t.event` is the source `Event`.
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

"""
    ReducedExplanation

Structured explanation of a `ReducedTransition`: the source mark, the
reduced `key`/`kind` (`:noop`/`:cross`/`:fork`), the summed `Φ` (`Φ_u`), and
the `TransitionExplanation` of every full `KLITransition` that collapsed
into it.

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
collapsed into it by `reduce_event_indicator` (`Φ_u = Σ φ_u`).
Throws `ArgumentError` if `rt.transitions` is empty, since there is no `event`
to report.
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

"""
    TermExplanation

Structured explanation of a `FilterTerm`: the `bucket`, the source mark, `Φ`
and the `ReducedExplanation`.
`mechanism` and `reason` are `nothing` except for `OutflowImbalanceTerm`.

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
(`term.event`), its `Φ`, and the explanation of its `ReducedTransition`
(`explain(term.reduced)`). For `OutflowImbalanceTerm`, also reports
`mechanism` (`:inflow_outflow_imbalance`) and `reason`
(`:fork_unobserved_at_regular_time`).
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
    :regular_flow      => "RegularFlow (continuous-time regular-event filter)",
    :singular_flow     => "SingularFlow (fixed observed data-event time only)",
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
                         " (NOT KLI's λ)" :
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
