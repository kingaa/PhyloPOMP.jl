# mgp_filter_ir.jl
# =============================================================================
# M05: Filter IR — classifying M04's Vector{ReducedTransition} into the
# proposal-independent TARGET filter structure KLI actually integrates:
#
#   RegularFlow         -- contributes to the continuous-time regular-event
#                   filter update (KLI's beta^reg_u = alpha^reg_u * pi_u
#                   driver, Eq. 45; mers_filter_suite.tex Sec. "Assembled
#                   regular filter", lines 791-870).
#   SingularFlow  -- only fires at a fixed, OBSERVED data-event time
#                   (KLI's beta^ev_{e,u} = alpha_u*pi_u, supported on evZ;
#                   mers_filter_suite.tex lines 871-873).
#   OutflowImbalanceTerm -- a PLACEHOLDER for probability mass that a
#                   REGULAR event's outcome would be KLI-incompatible with
#                   the observed pruned tree if it were realized as regular
#                   flow (an unobserved coalescence/branch point). This file
#                   does NOT compute the actual decay/driver weight (Eq. 47 /
#                   Appendix B, Eq. B2) -- that is explicitly the LATER
#                   "driver/boost/decay" milestone's job (M07 per the master
#                   project plan). Each `OutflowImbalanceTerm` here carries
#                   only PROVENANCE (the ReducedTransition it came from and
#                   its raw Phi) so the bucket is inspectable, not
#                   numerically meaningful yet.
#
# =============================================================================
# CORRECTION (M06, Part 0): this bucket was originally named `DecayTerm`.
# That name is WRONG and has been corrected here. `mers_filter_suite.tex`
# lines 760-762 and 841-844 are explicit that the fork/branch-point mass
# this bucket carries is accounted for through a THIRD mechanism -- the
# structural imbalance between the full-hazard regular outflow (alpha_u,
# unreduced) and the non-fork-only regular inflow (Eqs. 7-12, i.e. what this
# file calls `RegularFlow`) -- and that this is explicitly NOT the same
# mechanism as KLI's actual `lambda(t,x,y)` (sampling hazard + sub-threshold
# removal, `mers_filter_suite.tex` lines 833-836, the "Decay lambda"
# section). Naming this bucket `DecayTerm` invited exactly the double-
# counting mistake the tex warns against (lines 841-844: "Adding either to
# lambda would double-count"): a future milestone (M07) implementing lambda
# might see a `DecayTerm` in the Filter IR and naively `sum` it into its
# `lambda` total. Renamed to `OutflowImbalanceTerm` (and `FilterSpec.decay`
# to `FilterSpec.outflow_imbalance`) so that mechanism is structurally
# unmistakable: there is no type or field anywhere in this file with
# "decay" in its name, so a future implementer of `lambda`/`kli_decay`
# (`mgp_filter.jl:87-88`, still an unverified stub, untouched by this
# milestone) cannot casually `sum(FilterSpec.decay)` into their `lambda`
# computation without first noticing that field does not exist and reading
# why. See `handoffs/M06_audit_provenance.md`'s "Correction to M05" section
# for the full writeup.
# =============================================================================
#
# WHY Fork is excluded from RegularFlow (the empirical M03/M04 finding,
# formalized here):
#   M03 (handoffs/M03_full_kli_lowering.md, "Known failures") found, and M04
#   (handoffs/M04_reduce_m.md) corroborated, that `seir_naive.jl`'s REGULAR
#   birth proposal (`regular_part!`, e.g. k==2 for infection,
#   seir_naive.jl:150-156) unconditionally performs the CrossDeme swap and
#   NEVER realizes InlineSameDeme or Fork; `mers_naive.jl`'s `regular_part!`
#   (k==1 for transmission_cc, mers_naive.jl:187-190; k==3/4 for
#   transmission_hc, mers_naive.jl:195-205) shows the identical pattern.
#   Fork (`fork!`) is constructed ONLY inside `singular_part!`
#   (seir_naive.jl:88,90; mers_naive.jl's Node branch, lines 60-125), i.e.
#   at the FIXED times the input tree dictates a branch point.
#
#   `mers_filter_suite.tex` explains WHY this is principled, not just an
#   implementation quirk of the naive proposal:
#     - lines 760-762 (TCC/THH): "The missing mass corresponds to unobserved
#       within-camel branch points and is accounted for through the regular
#       birth inflow/outflow balance (Sec 8), not through lambda."
#     - lines 833-836 (the explicit "Decay lambda" formula): lambda(t,x,y) =
#       [sampling hazard: chi_C*I_C + chi_H*I_H] + [sub-threshold removal:
#       gamma_C*I_C*1{I_C<=ell_C} + gamma_H*I_H*1{I_H<=ell_H}] -- note this
#       formula has NO fork/branch-point term in it at all.
#     - lines 841-844 ("No separate branch-point hazard"): "The within-deme
#       unseen branch-point loss arises automatically from Eqs 7-8 inflow
#       against the TCC/THH share of the (full-birth) outflow ... Adding
#       either to lambda would double-count."
#     - lines 846-865 (the driver box): births are kept FULL in the reduced
#       regular exit rate (alpha_TCC+alpha_THH+alpha_THC+alpha_TCH, no
#       ell-dependent reduction), while the INFLOW terms I_7..I_12 (Eqs
#       7-12, lines 735-789) are built ONLY from the non-fork
#       (noop/cross) reduced cases. The imbalance between the full outflow
#       and the non-fork-only inflow IS where the branch-point/fork mass
#       goes -- it is real probability mass that must be accounted for as an
#       outflow/inflow imbalance the driver-boost-decay algebra resolves,
#       but it is NOT realized as a regular-time genealogy-changing move
#       (doing so would require an unobserved coalescence at a time the
#       fixed, observed tree says nothing branches), and it is NOT part of
#       lambda (lines 833-836's formula has no branch-point term, and lines
#       841-844 says adding one would double-count).
#
#   In this milestone's Filter IR vocabulary we call that excluded mass an
#   `OutflowImbalanceTerm` -- a bucket name that names the ACTUAL mechanism
#   the tex identifies (the inflow/outflow imbalance of Secs 7-8, NOT
#   lambda) -- WITHOUT asserting it is literally identical to the tex's
#   explicit lambda(t,x,y) formula (Sec 8.1, lines 833-844). Working out the
#   exact algebraic form of how this imbalance resolves (e.g. whether it
#   nets out purely via the driver's outflow term, or needs an explicit
#   correction elsewhere in Appendix B's general construction) is exactly
#   what the later driver/boost/decay milestone (M07) must derive; this file
#   only records that the mass exists, where it came from, and that it must
#   NOT be double-counted as RegularFlow, and must NOT be silently summed
#   into a future `lambda`/`kli_decay` computation either.
#
# Explicitly OUT of scope for this milestone (M05) and left out of scope by
# M06's Part 0 correction (naming only, no new math):
#   - The actual decay WEIGHT formula (Eq. 47 / B2, KLI's lambda) --
#     `OutflowImbalanceTerm`'s `Φ` field is only the raw
#     `ReducedTransition.Φ` it was built from, clearly documented as a
#     placeholder, not a computed decay/imbalance contribution.
#   - RC/RH/SC/SH (removal/sampling, r=(0,0)) -- M03 already scoped these
#     out (`full_transitions` throws on DEATH/SAMPLE/NEUTRAL events), so
#     `reduced_transitions` never produces input for them; this file's
#     classifier only ever receives BIRTH/MIGRATION-event ReducedTransitions.
#   - Proposal logic (pi_u itself, `apply_move!`'s `q`) -- FilterSpec
#     represents the proposal-INDEPENDENT target filter structure only.
#   - mgp_filter.jl's stubs and the 8 hand-coded filter modules (read-only
#     oracles, not touched).
#
# Primary sources: mers_filter_suite.tex (cited by line range above);
#   src/examples/seir_naive.jl, src/examples/mers_naive.jl (cited by
#   file:line above); src/examples/mgp_reduce.jl (M04, consumed unmodified);
#   src/examples/mgp_filter.jl's docstrings (kli_select/kli_decay/
#   apply_move!/singular_update!, UNVERIFIED STUBS, cited only for the KLI
#   equation numbers they attribute to these concepts: Eqs 44-46 (select),
#   47/B2 (decay), 9/37/38/46 (move/boost), 14/39-41/B3 (singular)).
# =============================================================================

export FilterTerm, RegularFlow, SingularFlow, OutflowImbalanceTerm, FilterSpec,
       classify_filter_terms, filter_spec

"""
    FilterTerm

Abstract supertype for the three kinds of Filter IR term a `ReducedTransition`
(M04) can be classified into: `RegularFlow`, `SingularFlow`,
`OutflowImbalanceTerm`. Every concrete subtype carries, at minimum, the
source `Event`, the `ReducedTransition` it was built from (provenance back to
M04's `Φ` and, transitively, M03's `Vector{KLITransition}`), and the raw
`Φ::Rational{Int}` value copied from that `ReducedTransition` for convenient
access without unpacking `reduced`.
"""
abstract type FilterTerm end

"""
    RegularFlow(event, reduced, Φ)

A `ReducedTransition` that contributes to the continuous-time REGULAR filter
update (KLI's `beta^reg_u = alpha^reg_u * pi_u` driver, Eq. 45) -- i.e. it is a
genealogy outcome the filter's between-events proposal may actually realize.
Per this file's header and the M03/M04 findings it formalizes, only `:noop`
and `:cross` `ReducedTransition`s of a REGULAR (`event.regular == true`)
BIRTH/MIGRATION event are classified here; `:fork` never is.
"""
struct RegularFlow <: FilterTerm
    event   :: Event
    reduced :: ReducedTransition
    Φ       :: Rational{Int}
end

"""
    SingularFlow(event, reduced, Φ)

A `ReducedTransition` that only fires at a fixed, OBSERVED data-event time
(KLI's `beta^ev_{e,u} = alpha_u * pi_u`, supported on `evZ` only,
`mers_filter_suite.tex` lines 871-873). Two ways a `ReducedTransition` lands
here:

  1. It came from a SINGULAR event (`event.regular == false`): every one of
     its `ReducedTransition`s (regardless of `:noop`/`:cross`/`:fork` kind)
     can only ever be realized at the fixed time(s) the observed tree
     dictates for that event, since `event.regular == false` already means
     `mgp_filter.jl`'s `regular_step!` zeroes its contribution to the
     continuous driver at every OTHER time (`alpha = [ev.regular ?
     kli_hazard(...) : 0.0 ...]`, `mgp_filter.jl:166-167`). Since this is
     singular, whichever outcome the observed tree structure actually
     dictates at that instant is exactly consistent by construction -- no
     imbalance bucket is needed for a singular event's own reduced
     transitions.

  NOTE: no currently-defined SEIR/MERS event exercises this case through
  this milestone's pipeline -- see this file's "Gap" note below and
  `handoffs/M05_filter_ir.md`.
"""
struct SingularFlow <: FilterTerm
    event   :: Event
    reduced :: ReducedTransition
    Φ       :: Rational{Int}
end

"""
    OutflowImbalanceTerm(event, reduced, Φ, mechanism, reason)

PLACEHOLDER for the probability mass carried by a `ReducedTransition` that
CANNOT be realized as `RegularFlow` because doing so would require an
unobserved genealogy event (a coalescence/branch point) inconsistent with the
fixed, observed pruned tree between data-event times -- see this file's
header comment for the `mers_filter_suite.tex` line citations (760-762,
833-836, 841-844, 846-865) motivating why `:fork` reduced transitions of a
REGULAR event land here.

**Renamed from `DecayTerm` in M06 (Part 0 correction).** The old name
implied this bucket IS (or feeds) KLI's `lambda(t,x,y)` decay term -- it does
NOT. Per `mers_filter_suite.tex` lines 760-762 and 841-844, this fork/
branch-point mass is resolved through the structural imbalance between the
full-hazard regular outflow (`alpha_u`) and the non-fork-only regular inflow
(Eqs. 7-12, `RegularFlow` in this file's vocabulary) -- a THIRD mechanism,
distinct from both `RegularFlow` and `lambda`. The tex's own `lambda` formula
(lines 833-836: sampling hazard + sub-threshold removal) has no branch-point
term, and lines 841-844 explicitly warn: "Adding either to lambda would
double-count."

**This struct carries ONLY provenance, not a computed weight of any kind.**
`Φ` is copied VERBATIM from the source `ReducedTransition.Φ` (the raw sum of
`phi_u` over the collapsed full transitions, M04) -- it is NOT KLI's
`lambda(t,x,y)` (Eq. 47/B2, a rate-integrated quantity involving the
population hazard `alpha_u`, not just `Φ_u`), and it is NOT yet the resolved
inflow/outflow-imbalance correction either (that numeric derivation is
EXPLICITLY OUT OF SCOPE for M05/M06 -- see the "Next milestone" section of
`handoffs/M06_audit_provenance.md`, M07's job).

**`mechanism` (required field, structural guard):** a `Symbol` tagging WHICH
resolution mechanism this term's mass is understood to belong to.
  - `:inflow_outflow_imbalance` -- the ONLY value this milestone (M05/M06)
    ever constructs: the fork/branch-point mass resolved via the Eqs. 7-12
    inflow vs. full-hazard outflow imbalance (tex lines 760-762, 841-844).
  - A hypothetical `:lambda` value is DELIBERATELY NEVER produced here --
    it would mean this term's mass belongs to KLI's actual decay term
    instead, which per the tex citations above is a DIFFERENT, disjoint
    pool of mass (sampling hazard + sub-threshold removal only). If a
    future milestone ever needs to represent that pool, it should NOT reuse
    `OutflowImbalanceTerm` with `mechanism = :lambda` -- it should be a
    genuinely different type/bucket (e.g. a future `DecayTerm` computing
    KLI's `lambda`, distinct from this type), so that `mechanism` staying
    fixed at `:inflow_outflow_imbalance` here remains a true, checkable
    invariant of this type, not a discriminated union pretending to be one.

`reason` is a `Symbol` documenting WHY this term was excluded from
`RegularFlow`, for audit/inspection purposes (unchanged from M05):
  - `:fork_unobserved_at_regular_time` -- a `:fork` `ReducedTransition` of a
    REGULAR event (the only case this milestone constructs).
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

Container for one event's classified Filter IR: the full set of `RegularFlow`,
`SingularFlow`, and `OutflowImbalanceTerm` terms built from its
`ReducedTransition`s (M04). This is the proposal-INDEPENDENT target filter
structure -- it says WHICH bucket each reduced outcome belongs to, not HOW a
proposal kernel `pi_u` would sample among them (that is `apply_move!`'s job in
`mgp_filter.jl`, explicitly out of scope here).

**Renamed from `FilterSpec.decay` in M06 (Part 0 correction)** -- see
`OutflowImbalanceTerm`'s docstring and `handoffs/M06_audit_provenance.md`'s
"Correction to M05" section: the field is named `outflow_imbalance`, not
`decay`, precisely so that a future milestone implementing KLI's `lambda`
cannot find a field literally called `decay` on `FilterSpec` and assume it is
(or belongs in) `lambda` without reading why it is named otherwise.
"""
struct FilterSpec
    event             :: Event
    regular           :: Vector{RegularFlow}
    singular          :: Vector{SingularFlow}
    outflow_imbalance :: Vector{OutflowImbalanceTerm}
end

"""
    classify_filter_terms(event::Event, reduced::AbstractVector{ReducedTransition}) -> FilterSpec

Classify each `ReducedTransition` in `reduced` (M04's `reduced_transitions`
output for `event` at some `(ℓ, n)`) into `RegularFlow`, `SingularFlow`, or
`OutflowImbalanceTerm`, per this file's header rationale:

- If `event.regular` (a REGULAR BIRTH/MIGRATION event, between fixed data
  times): `:noop` and `:cross` reduced transitions -> `RegularFlow`; `:fork`
  -> `OutflowImbalanceTerm` (`mechanism = :inflow_outflow_imbalance`,
  `reason = :fork_unobserved_at_regular_time`). Matches
  `seir_naive.jl:146-167`'s `regular_part!` (k==1/2 for `infection`, k==3/4
  for `progression` -- none of the four branches ever calls `fork!`) and
  `mers_naive.jl`'s `regular_part!` (k==1 for TCC, mers_naive.jl:187-190;
  k==3/4 for THC, mers_naive.jl:195-205 -- never calls `fork!` either).
- If `!event.regular` (a SINGULAR BIRTH/MIGRATION event, fixed observed time
  only): every reduced transition -> `SingularFlow`, regardless of kind --
  see `SingularFlow`'s docstring for why no imbalance bucket applies to a
  singular event's own outcomes. NOTE: no current SEIR/MERS BIRTH/MIGRATION
  event has `regular == false` (SEIR's `sampling` and MERS's
  `sampling_c`/`sampling_h` are the only `regular == false` events in either
  model, and all three are `SAMPLE`-type, out of `full_transitions`'s scope
  per M03 -- so this branch is validated only against a synthetic
  hand-built `Event`, see `test/kli_filter_ir_test.jl` and
  `handoffs/M05_filter_ir.md`'s "Gap" note).

A lightweight provenance-integrity check is performed: every `KLITransition`
inside every `reduced[i].transitions` must have been produced FOR `event`
(`t.event === event`), so a caller cannot silently classify one event's
`ReducedTransition`s under a different `event` argument by mistake.
"""
function classify_filter_terms(event::Event, reduced::AbstractVector{ReducedTransition})
    for rt in reduced, t in rt.transitions
        t.event === event ||
            throw(ArgumentError("classify_filter_terms: a ReducedTransition's " *
                                 "provenance references event `$(t.event.name)`, " *
                                 "not the `event` argument `$(event.name)` -- " *
                                 "refusing to classify mismatched provenance."))
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

Convenience composition from an event/state pair straight to its classified
Filter IR: `reduced_transitions(event, ℓ, n; Q=Q) |> rts ->
classify_filter_terms(event, rts)`.
"""
filter_spec(event::Event, ℓ::AbstractVector{<:Integer}, n::AbstractVector{<:Integer};
            Q::Real = 1) =
    classify_filter_terms(event, reduced_transitions(event, ℓ, n; Q = Q))
