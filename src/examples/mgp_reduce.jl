# mgp_reduce.jl
# =============================================================================
# M04: Explicit Reduction over the event-indicator m.
#
#   Consumes M03's `full_transitions(event, ell, n; Q) -> Vector{KLITransition}`
#   (src/examples/mgp_transitions.jl) and groups its elements by the outcome
#   they produce on the REDUCED, color-only state (the `Coloring{D,N}` type in
#   src/coloring.jl -- no event-indicator `m` component), summing the exact
#   `Rational{Int}` phi_u compatibility factors within each group to obtain
#   Phi_u, per `mers_filter_suite.tex`'s "Coloring convention" paragraph and
#   Eq. (lines 417-432):
#
#     "For a transition of the reduced color-only process, let Phi_u(y-,y+)
#      denote the sum of the KLI compatibility factors phi_u over all
#      full-coloring transitions that collapse to the same reduced
#      transition... When several collapse under the reduction, Phi_u is
#      their sum, evaluated in closed form by the Chu-Vandermonde identity."
#
#   This is a purely MECHANICAL grouping pass over whatever `full_transitions`
#   returns -- it does NOT know about Q_u regular-vs-singular gating (M03
#   found and documented that e.g. Fork is unreachable at regular events in
#   seir_naive.jl's proposal, but that gating is later-milestone territory,
#   not this file's concern) and does NOT implement Filter IR classification
#   (regular/singular/decay -- that is M05).
#
#   Collapsing rule (derived from `swap!`'s actual semantics in
#   src/coloring.jl, and mers_filter_suite.tex:419-421's "on this reduced
#   state a same-deme swap is the identity, sigma^b_{CC}y=y ... even though
#   the corresponding operator is not an identity on KLI's full coloring"):
#     - IdentityTransition (no slot occupied, y'=y) and InlineSameDemeTransition
#       (one slot occupied, ancestral deme == post-event deme, so
#       swap!(y,d,d,b) is a no-op on the BitSet-based Coloring: delete!
#       followed by push! of the same lineage id into the same deme's set)
#       BOTH leave the reduced coloring unchanged. They collapse into ONE
#       reduced "no-op in deme d" transition, where d = event.from (the
#       ancestral deme shared by both -- Step A of the tex, lines 438-441,
#       "ancestral deme" is ALWAYS event.from for every occupied slot of a
#       given event, so this is well-defined and consistent regardless of
#       which slot was filled).
#     - CrossDemeTransition (one slot occupied, ancestral != target deme)
#       genuinely moves a lineage between demes on the reduced state -- a
#       DIFFERENT reduced outcome from "no change" for every distinct
#       (ancestral_deme, target_deme) pair. Never collapsed with the noop
#       group or with a CrossDeme of a different deme pair.
#     - ForkTransition (two-or-more slots occupied) ADDS a lineage to the
#       reduced coloring (a branch point) -- again a genuinely different
#       reduced outcome from a "move". Grouped by (ancestral_deme,
#       sorted(slot_demes)) so that, generically, two DIFFERENT Fork
#       saturations of the SAME event that happen to produce the identical
#       reduced outcome (same ancestral deme, same multiset of post-event
#       slot demes) WOULD correctly collapse -- this does not occur for any
#       current SEIR/MERS event (each has at most one Fork-classified
#       saturation: TCC/THH have s summing to r_u=2 in a single deme only
#       when that deme's r_d=2, and THC/TCH's only 2-slot saturation is
#       (1,1)), but the grouping key is generic, not a special case of
#       "there's only ever one Fork".
#
# Primary sources: mers_filter_suite.tex, "Coloring convention" paragraph
#   (lines 417-421) and the Phi_u/boost definition (lines 423-432); Step A
#   (lines 438-441, ancestral deme = event.from uniformly); `src/coloring.jl`
#   (`swap!`'s actual same-deme no-op behavior, read directly rather than
#   assumed); M03's `mgp_transitions.jl` (consumed, not modified).
# =============================================================================

export ReducedTransition, reduce_event_indicator, reduced_transitions

"""
    ReducedTransition

The result of collapsing one or more `KLITransition`s (M03) that produce the
identical outcome on the REDUCED, color-only state (no event-indicator `m`)
into a single reduced transition, per `mers_filter_suite.tex`'s
`Phi_u = sum phi_u` reduction (lines 423-432).

Fields:
- `key`         : an opaque, generic grouping key (a `Tuple`) -- two
                   `KLITransition`s collapse into the same `ReducedTransition`
                   iff their `key`s are `==`. Structure (for inspection):
                     - `(:noop, d)` for the shared "no reduced-coloring
                       change, ancestral/event deme `d`" class (collapses
                       `IdentityTransition` and `InlineSameDemeTransition`).
                     - `(:cross, d_from, d_to)` for a single-lineage move.
                     - `(:fork, d_anc, slot_demes_tuple)` for a branch point,
                       `slot_demes_tuple` a sorted `Tuple` (order-independent,
                       multiplicity-preserving) of post-event slot demes.
- `kind`        : `:noop`, `:cross`, or `:fork` -- the human-readable tag
                   redundant with `key[1]`, provided for convenient dispatch
                   without unpacking `key`.
- `Φ`           : `Rational{Int}`, the exact sum of `phi` over every
                   contributing `KLITransition` (KLI's `Phi_u`).
- `transitions` : `Vector{KLITransition}`, the full-coloring transitions that
                   collapsed into this reduced transition, in the order they
                   appeared in the input -- PROVENANCE, kept (not discarded)
                   because later audit/explain milestones need to answer
                   "which full transitions collapsed into this Phi_u".
"""
struct ReducedTransition
    key         :: Tuple
    kind        :: Symbol
    Φ           :: Rational{Int}
    transitions :: Vector{KLITransition}
end

"""
    reduced_key(t::KLITransition) -> Tuple

The generic grouping key for `t`: two `KLITransition`s reduce to the same
`ReducedTransition` iff `reduced_key` returns `==` values for both. See
`ReducedTransition`'s docstring for the key structure per subtype, and this
file's header comment for the mathematical justification (in particular, why
`IdentityTransition` and `InlineSameDemeTransition` share a key -- both leave
the `Coloring` unchanged, per `swap!`'s actual same-deme no-op behavior in
`src/coloring.jl`).
"""
reduced_key(t::IdentityTransition) = (:noop, t.event.from)
reduced_key(t::InlineSameDemeTransition) = (:noop, t.deme)
reduced_key(t::CrossDemeTransition) = (:cross, t.ancestral_deme, t.target_deme)
reduced_key(t::ForkTransition) =
    (:fork, t.ancestral_deme, Tuple(sort(t.slot_demes)))

reduced_kind(t::IdentityTransition) = :noop
reduced_kind(t::InlineSameDemeTransition) = :noop
reduced_kind(t::CrossDemeTransition) = :cross
reduced_kind(t::ForkTransition) = :fork

"""
    reduce_event_indicator(transitions::Vector{KLITransition}) -> Vector{ReducedTransition}

Group `transitions` (M03's `full_transitions` output, for one event/state
pair) by `reduced_key`, and sum `phi` (exact `Rational{Int}` arithmetic)
within each group to obtain `Phi_u` for that reduced transition, per
`mers_filter_suite.tex`:423-432. This is a purely mechanical grouping pass:
it does not classify Filter-IR terms (regular/singular/decay -- M05), and it
does not gate transitions by `Q_u` reachability (M03's documented, deferred
concern) -- it operates on whatever `full_transitions` happened to return.

Returns one `ReducedTransition` per distinct reduced outcome, in the order
those outcomes were first encountered in `transitions`. Every input
`KLITransition` contributes to exactly one output group (no transition is
dropped or double-counted); `sum(rt.Φ for rt in result) ==
sum(t.phi for t in transitions)` always holds (both being a sum over the
full input transitions, just partitioned differently).

Example (MERS THC, `r=(1,1)`, 4 full transitions -> 3 reduced transitions):
`IdentityTransition` (s=(0,0)) and `InlineSameDemeTransition` (s=(1,0), the
continuing camel parent) both collapse into `(:noop, 1)` (deme C); the
`CrossDemeTransition` (s=(0,1), camel-to-human) and `ForkTransition`
(s=(1,1), the spillover branch point) each remain their own singleton
reduced transition.
"""
function reduce_event_indicator(transitions::AbstractVector{<:KLITransition})
    order = Tuple[]
    groups = Dict{Tuple,Vector{KLITransition}}()
    for t in transitions
        k = reduced_key(t)
        if !haskey(groups, k)
            groups[k] = KLITransition[]
            push!(order, k)
        end
        push!(groups[k], t)
    end
    result = Vector{ReducedTransition}(undef, length(order))
    for (i, k) in enumerate(order)
        members = groups[k]
        Φ = sum(t.phi for t in members)
        result[i] = ReducedTransition(k, reduced_kind(members[1]), Φ, members)
    end
    return result
end

"""
    reduced_transitions(event::Event, ell, n; Q = 1) -> Vector{ReducedTransition}

Convenience composition of M03's `full_transitions` and this file's
`reduce_event_indicator`: `full_transitions(event, ell, n; Q=Q) |>
reduce_event_indicator`. Goes directly from an event/state pair to the
grouped, Phi_u-summed reduced transitions in one call.
"""
reduced_transitions(event::Event, ℓ::AbstractVector{<:Integer},
                     n::AbstractVector{<:Integer}; Q::Real = 1) =
    reduce_event_indicator(full_transitions(event, ℓ, n; Q = Q))
