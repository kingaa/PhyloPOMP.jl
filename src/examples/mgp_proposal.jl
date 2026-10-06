# Proposal-strategy tags, and the driver and boost compositions.

export ProposalStrategy, NaiveProposal, SoftProposal, GuidedProposal, HardProposal,
       boost, driver

"""
    ProposalStrategy

Tag for a hand-designed proposal family: `NaiveProposal`, `SoftProposal`, `GuidedProposal`, `HardProposal`.
`π_u` is not computed generically.
"""
abstract type ProposalStrategy end

"Naive family: lineage-count-proportional selection."
struct NaiveProposal <: ProposalStrategy end

"Soft family: `SoftSEIR`, `SoftMERS`."
struct SoftProposal <: ProposalStrategy end

"Guided family: π from a reverse-time sweep over the deme process."
struct GuidedProposal <: ProposalStrategy end

"Hard family: `HardSEIR`, `HardMERS`."
struct HardProposal <: ProposalStrategy end

"""
    driver(α::Real, π::Real) -> Real

Driver rate `β_u = α_u · π_u`.
"""
driver(α::Real, π::Real) = α * π

"""
    boost(Φ::Rational{Int}, π::Real) -> Real

Boost `B_u = Φ_u / π_u`.
`π` is the total probability of proposing the reduced outcome, marginalized over within-outcome choices.
"""
boost(Φ::Rational{Int}, π::Real) = Φ / π

"""
    boost(rt::ReducedTransition, π::Real) -> Real

Same as `boost(rt.Φ, π)`.
"""
boost(rt::ReducedTransition, π::Real) = boost(rt.Φ, π)
