# Filter IR: classifies reduced transitions into regular, singular and outflow-imbalance terms.
#
# Fork outcomes of a regular event are never regular flow: the filter realizes forks
# only at observed branch points. They are classed as OutflowImbalanceTerm.

export FilterTerm, RegularFlow, SingularFlow, OutflowImbalanceTerm, FilterSpec,
       classify_filter_terms, filter_spec

"""
    FilterTerm

Abstract supertype of `RegularFlow`, `SingularFlow` and `OutflowImbalanceTerm`.
Each carries `event`, `reduced` and `Φ`.
"""
abstract type FilterTerm end

"""
    RegularFlow(event, reduced, Φ)

A `:noop` or `:cross` outcome of a regular BIRTH/MIGRATION event.
`:fork` is never regular flow.
"""
struct RegularFlow <: FilterTerm
    event   :: Event
    reduced :: ReducedTransition
    Φ       :: Rational{Int}
end

"""
    SingularFlow(event, reduced, Φ)

An outcome of a singular event (`event.regular == false`).
Realized only at the observed data-event time.
"""
struct SingularFlow <: FilterTerm
    event   :: Event
    reduced :: ReducedTransition
    Φ       :: Rational{Int}
end

"""
    OutflowImbalanceTerm(event, reduced, Φ, mechanism, reason)

Fork mass of a regular event that cannot be regular flow.
`Φ` is copied from the `ReducedTransition`.
`mechanism` is always `:inflow_outflow_imbalance`.
`reason` is `:fork_unobserved_at_regular_time`.
"""
struct OutflowImbalanceTerm <: FilterTerm
    event     :: Event
    reduced   :: ReducedTransition
    Φ         :: Rational{Int}
    mechanism :: Symbol
    reason    :: Symbol
end

"""
    FilterSpec(event, regular, singular, outflow_imbalance)

One event's classified terms: `regular`, `singular` and `outflow_imbalance` vectors.
"""
struct FilterSpec
    event             :: Event
    regular           :: Vector{RegularFlow}
    singular          :: Vector{SingularFlow}
    outflow_imbalance :: Vector{OutflowImbalanceTerm}
end

"""
    classify_filter_terms(event::Event, reduced::AbstractVector{ReducedTransition}) -> FilterSpec

Classify each reduced transition.
If `event.regular`, `:noop`/`:cross` go to `RegularFlow` and `:fork` to `OutflowImbalanceTerm`.
Otherwise all go to `SingularFlow`.
Throws `ArgumentError` if a transition's event is not `event`.
"""
function classify_filter_terms(event::Event, reduced::AbstractVector{ReducedTransition})
    for rt in reduced, t in rt.transitions
        t.event === event ||
            throw(ArgumentError("classify_filter_terms: a ReducedTransition " *
                                 "comes from event `$(t.event.name)`, " *
                                 "not the `event` argument `$(event.name)`."))
    end

    regular           = RegularFlow[]
    singular          = SingularFlow[]
    outflow_imbalance = OutflowImbalanceTerm[]

    if event.regular
        for rt in reduced
            if rt.kind == :noop || rt.kind == :cross
                push!(regular, RegularFlow(event, rt, rt.Φ))
            elseif rt.kind == :fork
                push!(outflow_imbalance,
                      OutflowImbalanceTerm(event, rt, rt.Φ, :inflow_outflow_imbalance,
                                            :fork_unobserved_at_regular_time))
            else
                throw(ArgumentError("classify_filter_terms: unrecognized " *
                                     "ReducedTransition.kind `$(rt.kind)`"))
            end
        end
    else
        for rt in reduced
            push!(singular, SingularFlow(event, rt, rt.Φ))
        end
    end

    FilterSpec(event, regular, singular, outflow_imbalance)
end

"""
    filter_spec(event::Event, ℓ, n; Q = 1) -> FilterSpec

Compose `reduced_transitions` and `classify_filter_terms` for `(event, ℓ, n)`.
"""
filter_spec(event::Event, ℓ::AbstractVector{<:Integer}, n::AbstractVector{<:Integer};
            Q::Real = 1) =
    classify_filter_terms(event, reduced_transitions(event, ℓ, n; Q = Q))
