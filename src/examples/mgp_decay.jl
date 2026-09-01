# mgp_decay.jl
# =============================================================================
# M07 (Part A): the decay term lambda(t,x,y) -- KLI Eq. 47 / Appendix B Eq.
# B2 -- generalized from mers_filter_suite.tex's MERS-specific closed form to
# an arbitrary DEATH/SAMPLE `Event`.
#
#   M03 (`mgp_transitions.jl`) explicitly scoped DEATH/SAMPLE/NEUTRAL events
#   OUT of `full_transitions` because their KLI compatibility is
#   decay/rate-driven, not saturation-driven (production vector r=(0,...,0)
#   or, for SEIR's `sampling`, a production vector that still does not
#   participate in the same binomial-ratio compatibility story -- see below).
#   That decay math is implemented here for the first time, generically.
#
# Primary source: `mers_filter_suite.tex`
#   - "Decay lambda" (subsection, grep `\subsection{Decay \$\lambda\$}`):
#       lambda(t,x,y) = [sampling hazard: chi_C I_C + chi_H I_H]
#                     + [sub-threshold removal: gamma_C I_C*1{I_C<=ell_C}
#                                              + gamma_H I_H*1{I_H<=ell_H}]
#   - "RC and RH: removal (death)" (event-specific derivation): "the
#     *strict* I_d>ell_d is the driver-outflow condition ... When I_d=ell_d
#     the outflow has no compatible target and becomes the sub-threshold
#     decay gamma_d I_d 1{I_d<=ell_d} in lambda (Sec 8)."
#   - "SC and SH: destructive sampling": "the full camel sampling rate
#     chi_C I_C enters lambda because the exact likelihood conditions on no
#     sample at t not in evZ. SH is symmetric."
#   Also matches `mgp_filter.jl`'s (still-unfilled) `kli_decay` stub
#   docstring verbatim: "lambda = integral(alpha_Sample dx') +
#   integral(alpha_Recov*1_{I<=ellI} dx'); general definition in Appendix B,
#   Eq. B2."
#
# =============================================================================
# THE GENERALIZATION (from the two-deme MERS formula to an arbitrary Event)
# =============================================================================
#
# For a DEATH event `u` with hazard alpha_u(x,theta) acting on deme
# `d = event.from` (the only deme a DEATH event touches -- `_decode_move`'s
# `:chop` branch sets `into = Int[]`, no production, single source deme):
#
#     decay_u(x,theta,ell,n) = alpha_u(x,theta) * 1{ n[d] <= ell[d] }
#
# `n[d]` is the POST-EVENT occupancy of deme `d` and `ell[d]` its pruned
# (tracked-lineage) count, in the same convention M02-M06 use throughout
# (`n_d`/`ell_d`, post-event). Removing an individual needs the destination
# occupancy to satisfy `n[d]-1 >= ell[d]`, i.e. `n[d] >= ell[d]+1`; the
# complement `n[d] <= ell[d]` is exactly the sub-threshold condition where no
# compatible outflow exists at all, so the *entire* hazard decays. Given the
# model invariant `ell[d] <= n[d]` (never violated post-event, asserted
# throughout M02-M06), `n[d] <= ell[d]` collapses to `n[d] == ell[d]` exactly
# -- matching the tex's own parenthetical at line 841:
# `(ind_{I_C<=ell_C} === ind_{I_C=ell_C})`.
#
# For a SAMPLE event `u` (regardless of production vector -- see the
# "SEIR sampling vs MERS sample_remove" note below), its FULL hazard
# contributes to lambda UNCONDITIONALLY:
#
#     decay_u(x,theta,ell,n) = alpha_u(x,theta)
#
# Reasoning (from the tex, generalized): a SAMPLE event is singular in every
# shipped model (`event.regular == false`) -- by construction it only ever
# fires at the fixed time(s) the observed data dictates (KLI's `evZ`). Any
# occurrence of the underlying population jump on an open interval
# (`t not in evZ`) is therefore, by definition, NOT one of those fixed
# observed times, so it cannot correspond to a real sample in the data --
# the exact likelihood must condition on "no sample-type jump occurred here",
# and the full rate at which that (never-observed) jump COULD have occurred
# is exactly what decays the weight. No indicator/threshold gates this: it is
# not a question of whether enough tracked lineages remain (a DEATH-style
# compatibility question), it is that the event class itself is entirely
# absent from the regular-time driver (mgp_filter.jl:166-167 already zeroes
# `alpha`/`pi` for any `!ev.regular` event on the regular driver) and must
# reappear somewhere -- it reappears in full, here, in lambda.
#
# NEUTRAL events (`waning`, MERS's `birth_c/birth_h/death_c/death_h`) do NOT
# contribute to lambda at all -- confirmed against both naive filters'
# `event_rates!`/decay closures (mers_naive.jl:159-161, seir_naive.jl:120),
# neither of which references B_c/B_h/omega anywhere in their decay
# expression. This matches KLI's own definition: a NEUTRAL event never
# touches the coloring (`type` has "no coloring op" even conceptually), so
# there is no compatibility question to fail and hence nothing to decay.
# `total_decay` below simply never visits NEUTRAL/BIRTH/MIGRATION events.
#
# =============================================================================
# THE "SEIR sampling vs MERS sample_remove" STRUCTURAL DIFFERENCE (resolved)
# =============================================================================
#
# M03 found (handoffs/M03_full_kli_lowering.md Sec.3) that `r_u=(0,...,0)` is
# NOT universal for SAMPLE events: MERS's `sampling_c`/`sampling_h`
# (`move=sample_remove(...)`, `mgp_mers.jl:17-18`) have `r=(0,0)` (the host is
# consumed, not born into a slot -- tex line 527), but SEIR's `sampling`
# (`move=sample(I)`, `mgp.jl:105-106` / `mgp_macro.jl:97-100`) has `r=(0,1)`,
# a genuine production slot at the sampled deme (the inserted leaf), because
# SEIR's sample is NON-destructive: `Event(:sampling, ..., Pair{Symbol,Int}[],
# ...)` -- `Δ=[]`, the host is NOT removed from the population.
#
# RESOLUTION: this difference changes what happens AT the singular (observed)
# event time -- SEIR's sample inserts a leaf onto a lineage that continues to
# exist in the population afterward (r=(0,1), a real coloring slot for the
# new leaf's parent identity, handled by `singular_update!`/Eq.14/39-41 at
# M07's later end-to-end milestone), whereas MERS's sample_remove simultaneously
# inserts the leaf AND removes the host from circulation (r=(0,0), a compound
# sample/death mark per tex line 388's footnote) -- but it does NOT change the
# shape of the LAMBDA (decay) contribution derived above. Both are SAMPLE-type,
# `event.regular == false`, marks: their decay contribution is the full hazard,
# UNCONDITIONALLY, regardless of `event.r`. The production vector only matters
# for what the event does TO THE GENEALOGY when it fires at its dictated
# singular time (a saturation/compatibility question, M02/M03 territory,
# already explicitly out of `full_transitions`'s scope by `event.type`, not by
# `r`) -- it plays no role in the between-observations decay accumulation,
# which only ever asks "could this jump type have fired here, unobserved" and
# answers with the bare hazard. Confirmed directly against the actual code
# (not just argued from the tex): `seir_naive.jl`'s `event_rates!` decay
# expression (`seir_naive.jl:120`) is `ψ*I + χ*I + ...` -- the `ψ*I` term is
# EXACTLY `PhyloPOMP.SEIR`'s `sampling` event's hazard `θ.ψ*x.I`
# (`mgp.jl:197`/`mgp_macro.jl`'s `_rewrite`), evaluated and added
# UNCONDITIONALLY, the identical shape as MERS's `chi_c*Ic`/`chi_h*Ih` terms
# in `mers_naive.jl:159` despite the differing `r`.
#
# ONE GENUINE GAP found while confirming this (not fabricated, reported
# per M07's instructions): `seir_naive.jl`'s decay expression ALSO includes a
# `chi*I` term (`seir_naive.jl:120`, using a SECOND rate constant `chi`,
# distinct from `psi`) for a "destructive sample" jump-type
# (`seir_naive.jl:54-64`'s `k==2` branch inside `singular_part!`) that has NO
# corresponding `Event` anywhere in `PhyloPOMP.SEIR`'s `MGPModel` -- the `@mgp
# SEIR` block (`mgp.jl:188-198`) declares `chi` as a model PARAMETER but never
# uses it in any `@event`'s `rate=`. In every shipped test configuration
# `chi = 0.0` (`seir_funs.jl:182,190`), so this never numerically fires, but
# structurally it means `total_decay(SEIR, ...)` (built only from `SEIR`'s
# actual `Event`s) cannot currently reproduce the `chi*I` term of
# `seir_naive.jl`'s decay even in principle -- it is not a decay-formula
# disagreement (the `psi*I` term this file DOES compute matches exactly), it
# is a Population-IR representation gap (a jump type the naive filter
# simulates that the `@mgp` model never declared as its own `Event`). Fixing
# that gap would mean adding a second SAMPLE-type `Event` to `mgp.jl`'s
# `SEIR`/`SEIR_REFERENCE` tables -- out of scope here (read-only per this
# milestone's constraints on `mgp.jl`), flagged for a later milestone instead.
#
# =============================================================================
# A DISCREPANCY FOUND AGAINST THE NAIVE FILTERS' ACTUAL RUNNING CODE
# (reported, not silently resolved -- see handoffs/M07_driver_boost_decay.md)
# =============================================================================
#
# `mers_naive.jl:159-161` and `seir_naive.jl:120` do NOT literally compute
# `gamma_d*I_d*1{I_d<=ell_d}` for their DEATH-type decay contribution. Given
# the models' own invariant `I_d >= ell_d`, their actual expression
#     gamma_d*ell_d + indicator(I_d <= ell_d, gamma_d*(I_d-ell_d))
# algebraically reduces to `gamma_d*ell_d`, UNCONDITIONALLY (the indicator
# term is identically zero: `I_d<=ell_d` forces `I_d==ell_d`, at which point
# `I_d-ell_d==0` regardless). This EQUALS this file's
# `gamma_d*I_d*1{I_d<=ell_d}` only at the boundary `I_d==ell_d` (where both
# reduce to `gamma_d*ell_d`); it EXCEEDS this file's formula by exactly
# `gamma_d*ell_d` whenever `I_d>ell_d` (where this file's indicator makes the
# whole term 0, but the naive code still adds `gamma_d*ell_d`). This module
# deliberately implements the tex's literal "Decay lambda" formula (and
# `mgp_filter.jl`'s `kli_decay` docstring, which cites the identical
# `1_{I<=ellI}` form) rather than silently matching the naive filters'
# code, because (a) that IS the formula this milestone was asked to
# generalize, and (b) the naive filters are single-PARTICLE (one simulated
# population trajectory) importance samplers, not a direct implementation of
# the PDE-level "Assembled regular filter" (tex lines 791-830) that the boxed
# lambda formula is stated for -- the PDE's matched inflow/outflow pairs
# (e.g. tex line 803's `gamma_C(I_C+1)*ind_{I_C>=ell_C}*w(t,x+e_C,y)` inflow
# against the FULL `gamma_C*I_C*ind_{I_C>ell_C}` outflow, tex line 804) can
# cancel the "removing a tracked lineage" mass against a matching inflow term
# from a NEIGHBORING population state `x+e_C`; a single simulated trajectory
# has no such neighboring state to borrow probability mass from, so its
# per-particle importance weight must instead decay continuously at the
# EXTRA rate `gamma_d*ell_d` (whether or not `I_d>ell_d`) to correctly
# condition against "one of the currently-tracked lineages died of an
# unobserved background death" -- a real possibility at every instant,
# independent of whether an untracked, compatible removal is ALSO possible.
# This is a plausible, internally-consistent reading, not a proven one: it is
# flagged here, and in the M07 handoff, as requiring advisor/M08 attention
# (M08's end-to-end numerical comparison against the naive filters' total
# log-likelihood is the concrete gate that would settle it either way).
# =============================================================================

export DecayContribution, decay_contribution, total_decay

"""
    DecayContribution

Structured, provenance-carrying result of `decay_contribution` for a single
DEATH/SAMPLE `Event`.

Fields
- `event`  : the source `Event` (`event.type in (DEATH, SAMPLE)`).
- `alpha`  : the evaluated population hazard `alpha_u(x,theta) =
             event.hazard(x,theta)`.
- `gate`   : the `{0,1}`-valued (`Rational{Int}`, so it composes exactly with
             `alpha` when `alpha` is itself a `Rational`) compatibility
             indicator multiplying `alpha` -- `1//1` unconditionally for
             SAMPLE events; `n[d]<=ell[d]` (`d=event.from`) for DEATH events.
- `value`  : `alpha * gate`, this event's contribution to `lambda`.
- `reason` : `:sampling_hazard_full` or `:sub_threshold_removal`, documenting
             which branch of the generalization (see this file's header)
             produced `gate`.
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

The per-event lambda(t,x,y) contribution (KLI Eq. 47 / Appendix B Eq. B2),
generalized from `mers_filter_suite.tex`'s MERS-specific closed form (see
this file's header for the full derivation and citations) to an arbitrary
`event.type in (DEATH, SAMPLE)`.

- `event.type == SAMPLE`: `gate = 1//1` unconditionally -- the full hazard
  `event.hazard(x,θ)` enters lambda (tex: "the full camel sampling rate
  chi_C I_C ... enters lambda"). Independent of `event.r` -- see this file's
  header for why the SEIR (`r=(0,1)`) vs MERS (`r=(0,0)`) production-vector
  difference does not change this.
- `event.type == DEATH`: `gate = (n[d] <= ℓ[d]) ? 1//1 : 0//1`, where
  `d = event.from` is the (only) deme a DEATH event touches. Requires
  `ℓ[d] <= n[d]` (the model invariant); throws `ArgumentError` if violated,
  matching M02/M03's "fail loud on invalid input" convention.

`x`/`θ` are passed through unchanged to `event.hazard` (KLI's `alpha_u`,
`mgp.jl`'s `Event.hazard` field, `(x,θ) -> Float64` in general, but exact
`Rational` arithmetic is preserved end-to-end if `x`/`θ`'s relevant fields
are themselves `Rational` -- used throughout `test/kli_decay_test.jl` for
exact hand-verification, per this milestone's "exact arithmetic, not
floating point" convention).

Throws `ArgumentError` for any `event.type` other than `DEATH`/`SAMPLE` --
`BIRTH`/`MIGRATION` events belong to `mgp_filter_ir.jl`'s `RegularFlow`/
`SingularFlow`/`OutflowImbalanceTerm` classification, not lambda; `NEUTRAL`
events never contribute to lambda at all (see this file's header).
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

`lambda(t,x,y) = Σ_u decay_contribution(u, x, θ, ℓ, n).value`, summed over
every `event ∈ model.events` with `event.type in (DEATH, SAMPLE)`
(`BIRTH`/`MIGRATION`/`NEUTRAL` events are skipped -- see this file's header
for why NEUTRAL never contributes).

This is the model-level generalization of `mers_filter_suite.tex`'s
`lambda(t,x,y) = chi_C I_C + chi_H I_H + gamma_C I_C*1{I_C<=ell_C} +
gamma_H I_H*1{I_H<=ell_H}` and of `mgp_filter.jl`'s (still-unfilled)
`kli_decay` stub docstring, to an arbitrary model's event table.
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
