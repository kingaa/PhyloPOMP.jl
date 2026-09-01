# mgp_proposal.jl
# =============================================================================
# M07 (Part B): the proposal-strategy vocabulary (naive/soft/guided/hard,
# per M00's already-made decision to KEEP this taxonomy rather than invent a
# new one -- `handoffs/M00_reconnaissance.md`'s "Known failures" section
# flagged this as an open design fork; the master project resolved it in
# favor of formalizing the existing names) plus the two generic compositions
# every proposal kernel, whichever strategy designed it, must satisfy:
#
#   driver(alpha, pi) = alpha * pi     -- KLI's beta_u^reg = alpha_u^reg*pi_u
#                                          (Eq. 45; `mers_filter_suite.tex`
#                                          line 849's notebox: "the per-mark
#                                          driver is KLI's beta^reg_u =
#                                          alpha^reg_u*pi_u").
#   boost(Phi, pi)   = Phi / pi        -- KLI's B_u = Psi_u = phi_u/pi_u
#                                          (Theorem 5; `mgp_filter.jl`'s
#                                          `apply_move!` docstring: "boost
#                                          B = Psi_u = phi_u/pi_u (Theorem
#                                          5)"; `mers_filter_suite.tex` lines
#                                          427-428: "B_u(y-,y+) =
#                                          Phi_u(y-,y+)/pi_u(y-,y+)").
#
# NOT done here (explicitly out of scope, per this milestone's spec): deriving
# pi_u itself for soft/guided/hard. Those are legitimate, already-tested,
# hand-designed importance kernels (`seir_soft.jl`/`seir_guided.jl`/
# `seir_hard.jl` and their MERS mirrors) that the master project's decision
# (M00) says to KEEP, not replace with a generic derivation. This file gives
# them a named place in the compiler's IR (`ProposalStrategy` and its four
# subtypes) so a later milestone (M08+) can tag which concrete kernel
# implements which strategy, and gives the driver/boost composition that
# EVERY strategy -- however pi_u is built -- must satisfy generically.
#
# Primary sources: `mers_filter_suite.tex` lines 417-432 ("Coloring
# convention" -- the Phi_u/boost definition) and lines 846-865 ("Assembled
# filter equation and the driver beta"); `mgp_filter.jl`'s `kli_select`/
# `apply_move!` docstrings (the concrete landing spots this scaffolding feeds,
# still unfilled stubs, untouched by this milestone).
# =============================================================================

export ProposalStrategy, NaiveProposal, SoftProposal, GuidedProposal, HardProposal,
       boost, driver

"""
    ProposalStrategy

Abstract supertype tagging WHICH of the project's four hand-designed
importance-kernel families (`pi_u` construction) a proposal belongs to --
`NaiveProposal`, `SoftProposal`, `GuidedProposal`, `HardProposal`. Per M00's
decision (`handoffs/M00_reconnaissance.md`), this is the project's
proposal-design vocabulary going forward, not a placeholder for a to-be-
invented generic derivation: `pi_u` itself is NOT computed generically
anywhere in this compiler (that is the eight hand-coded kernels'
already-tested job, one per {SEIR,MERS} x {naive,soft,guided,hard}) --
these four singleton types exist only so a `ReducedTransition`/`FilterTerm`
processing pipeline (M08+) can carry, inspect, and dispatch on WHICH
strategy produced a given `pi_u`, alongside the `driver`/`boost`
composition below that is strategy-INDEPENDENT (every strategy's `pi_u`
composes with `alpha_u`/`Phi_u` the same way).
"""
abstract type ProposalStrategy end

"Tag for the `naive` proposal family (`seir_naive.jl`/`mers_naive.jl`):
lineage-count-proportional, non-anticipatory selection (see `NaiveSEIR`'s
module docstring, `seir_naive.jl:1-9`)."
struct NaiveProposal <: ProposalStrategy end

"Tag for the `soft` proposal family (`seir_soft.jl`/`mers_soft.jl`)."
struct SoftProposal <: ProposalStrategy end

"Tag for the `guided`/`relhaz` proposal family (`seir_guided.jl`/
`mers_guided.jl`) -- an anticipatory `pi_u` built from a reverse-time sweep
over the deme process (`mgp_filter.jl:153-155`'s `regular_step!` docstring:
\"builds an anticipatory pi from a reverse-time sweep ... to 'borrow
information from future events'\")."
struct GuidedProposal <: ProposalStrategy end

"Tag for the `hard` proposal family (`seir_hard.jl`/`mers_hard.jl`)."
struct HardProposal <: ProposalStrategy end

"""
    driver(α::Real, π::Real) -> Real

KLI's per-mark REGULAR driver rate `β_u^reg = α_u^reg · π_u` (Eq. 45;
`mers_filter_suite.tex` line 849's notebox), the rate at which the filter's
continuous-time SMC scheme proposes a jump of mark `u`: the population
hazard `α_u` (`event.hazard(x,θ)`, KLI §2.5) times the selection probability
`π_u` (`kli_select`'s target, KLI §4/Eqs. 44-46) that this particular
REGULAR reduced outcome is what gets proposed.

Trivial by design (`α*π`) -- the point of naming and testing it is making
the driver composition an explicit, inspectable part of the IR rather than
leaving it implicit inside eight different hand-coded `alpha.*pi` products
(`seir_naive.jl:142`, `mers_naive.jl:183`, and their `_soft`/`_guided`/
`_hard` siblings).
"""
driver(α::Real, π::Real) = α * π

"""
    boost(Φ::Rational{Int}, π::Real) -> Real

KLI's per-mark boost `B_u = Ψ_u = Φ_u/π_u` (Theorem 5; `mgp_filter.jl`'s
`apply_move!` docstring: "boost B = Ψᵤ = ϕᵤ/πᵤ (Theorem 5)";
`mers_filter_suite.tex` lines 427-428: `B_u(y-,y+) = Φ_u(y-,y+)/π_u(y-,y+)`,
with the resulting identity `β_uB_u = α_uπ_u·(Φ_u/π_u) = α_uΦ_u` given at
line 432).

`Φ` is a reduced transition's target compatibility mass (M04's
`ReducedTransition.Φ`, an exact `Rational{Int}`); `π` is that SAME reduced
transition's TOTAL proposal probability under whichever `ProposalStrategy`
built it -- i.e. the full probability of proposing this reduced outcome,
already marginalized over any within-outcome degrees of freedom a concrete
kernel introduces (e.g. `seir_naive.jl`'s per-lineage uniform choice `q`;
see this file's header and `test/kli_proposal_test.jl`'s numerical
cross-check for why `π` must be this FULLY-marginalized probability, not
just the top-level categorical `pi[k]` a hand-coded kernel happens to name
`pi`).

`pi_u` itself is NOT derived here for any of the four `ProposalStrategy`
families -- see this file's header. This function is the trivial `Φ/π`
composition, made explicit and tested rather than left implicit inside eight
hand-coded files.
"""
boost(Φ::Rational{Int}, π::Real) = Φ / π

"""
    boost(rt::ReducedTransition, π::Real) -> Real

Convenience overload: `boost(rt.Φ, π)`, so a `ReducedTransition` (M04) or,
transitively, a `FilterTerm` (M05/M06, via its `.reduced`/`.Φ` fields) can
be boosted directly without unpacking `.Φ` first.
"""
boost(rt::ReducedTransition, π::Real) = boost(rt.Φ, π)
