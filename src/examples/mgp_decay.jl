# Decay λ(t,x,y).
# DEATH: α·1{n[d] ≤ ℓ[d]}, d = event.from.
# SAMPLE: the full α.
# BIRTH, MIGRATION, NEUTRAL: 0.

export DecayContribution, decay_contribution, total_decay

"""
    DecayContribution

Contribution of one DEATH or SAMPLE event to `λ`.

Fields
- `event`  : the source `Event`.
- `alpha`  : the hazard `event.hazard(x, θ)`.
- `gate`   : `1//1` or `0//1`.
- `value`  : `alpha * gate`.
- `reason` : `:sampling_hazard_full` or `:sub_threshold_removal`.
"""
struct DecayContribution
    event  :: Event
    alpha  :: Real
    gate   :: Rational{Int}
    value  :: Real
    reason :: Symbol
end

"""
    decay_contribution(event::Event, x, θ, ℓ::AbstractVector{<:Integer},
                        n::AbstractVector{<:Integer}) -> DecayContribution

The contribution of one event to `λ(t,x,y)`.

- SAMPLE: `gate = 1//1`.
- DEATH: `gate = 1//1` if `n[d] <= ℓ[d]`, else `0//1`, where `d = event.from`.

Exact if `x` and `θ` are `Rational`.
Throws `ArgumentError` for other event types, for `event.from == 0`, for `d` out of
range, and for `ℓ[d] > n[d]`.
"""
function decay_contribution(event::Event, x, θ,
                             ℓ::AbstractVector{<:Integer},
                             n::AbstractVector{<:Integer})
    event.type in (DEATH, SAMPLE) ||
        throw(ArgumentError("decay_contribution: event `$(event.name)` has " *
                             "type `$(event.type)`, not DEATH or SAMPLE -- " *
                             "lambda only accumulates DEATH/SAMPLE mass " *
                             "(BIRTH/MIGRATION belong to the Filter IR, " *
                             "NEUTRAL never contributes to lambda)."))

    α = event.hazard(x, θ)

    if event.type == SAMPLE
        gate   = one(Rational{Int})
        reason = :sampling_hazard_full
    else # DEATH
        d = event.from
        d >= 1 ||
            throw(ArgumentError("decay_contribution: DEATH event " *
                                 "`$(event.name)` has no source deme " *
                                 "(`event.from == 0`)."))
        (d <= length(ℓ) && d <= length(n)) ||
            throw(ArgumentError("decay_contribution: deme index $d " *
                                 "(event.from) exceeds length(ℓ)=$(length(ℓ)) " *
                                 "or length(n)=$(length(n))."))
        ℓd, nd = ℓ[d], n[d]
        ℓd <= nd ||
            throw(ArgumentError("decay_contribution: model invariant " *
                                 "violated -- ℓ[$d]=$ℓd > n[$d]=$nd."))
        gate   = nd <= ℓd ? one(Rational{Int}) : zero(Rational{Int})
        reason = :sub_threshold_removal
    end

    DecayContribution(event, α, gate, α * gate, reason)
end

"""
    total_decay(model::MGPModel, x, θ, ℓ::AbstractVector{<:Integer},
                n::AbstractVector{<:Integer}) -> Real

`λ(t,x,y)`: sum of `decay_contribution(...).value` over the DEATH and SAMPLE
events of `model`. Other event types are skipped.
"""
function total_decay(model::MGPModel, x, θ,
                      ℓ::AbstractVector{<:Integer},
                      n::AbstractVector{<:Integer})
    total = zero(Rational{Int})
    for event in model.events
        event.type in (DEATH, SAMPLE) || continue
        total += decay_contribution(event, x, θ, ℓ, n).value
    end
    total
end
