# Group full_transitions output by reduced (color-only) outcome and sum phi into Φ_u.
# IdentityTransition and InlineSameDemeTransition both leave the coloring unchanged.
# They share the key (:noop, d), d = event.from.

export ReducedTransition, reduce_event_indicator, reduced_transitions

"""
    ReducedTransition

Result of grouping KLITransitions that give the same outcome on the color-only state.

Fields:
- `key`: grouping key, a `Tuple`. Shapes: `(:noop, d)`, `(:cross, d_from, d_to)`,
  `(:fork, d_anc, slot_demes)` with `slot_demes` a sorted `Tuple`.
- `kind`: `:noop`, `:cross` or `:fork` (same as `key[1]`).
- `Φ`: exact `Rational{Int}` sum of `phi` over the group.
- `transitions`: the grouped KLITransitions, in input order.
"""
struct ReducedTransition
    key         :: Tuple
    kind        :: Symbol
    Φ           :: Rational{Int}
    transitions :: Vector{KLITransition}
end

"""
    reduced_key(t::KLITransition) -> Tuple

Grouping key for `t`; see `ReducedTransition` for the shapes.
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

Group `transitions` (`full_transitions` output for one event/state) by `reduced_key`.
`Φ` is the exact sum of `phi` per group.

Returns one `ReducedTransition` per distinct outcome, in first-seen order.
`sum(rt.Φ for rt in result) == sum(t.phi for t in transitions)`.

Example: MERS THC gives 4 transitions.
Identity and InlineSameDeme merge into `(:noop, 1)`, giving 3 reduced transitions.
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

`full_transitions(event, ℓ, n; Q)` followed by `reduce_event_indicator`.
"""
reduced_transitions(event::Event, ℓ::AbstractVector{<:Integer},
                     n::AbstractVector{<:Integer}; Q::Real = 1) =
    reduce_event_indicator(full_transitions(event, ℓ, n; Q = Q))
