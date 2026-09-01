# Handoff — M07: Driver/Boost/Decay

## Status
COMPLETE

## Objective
Two parts, Part A first (more novel/riskier). **Part A**: implement KLI's
decay term `λ(t,x,y)` (Eq. 47 / Appendix B Eq. B2) generically, for the
first time — `mers_filter_suite.tex`'s MERS-specific closed form generalized
to an arbitrary `DEATH`/`SAMPLE` `Event`, the class M03 explicitly scoped
out of `full_transitions` because their compatibility is decay/rate-driven,
not saturation-driven. **Part B**: formalize the existing
naive/soft/guided/hard proposal taxonomy (M00's decision) as an explicit
`ProposalStrategy` IR vocabulary, plus the generic `driver`
(`β_u=α_u·π_u`) / `boost` (`B_u=Φ_u/π_u`) composition every proposal
kernel — however `π_u` is built — must satisfy. `π_u` itself is NOT derived
generically for soft/guided/hard (explicitly out of scope, per-model
hand-designed kernels the project keeps). `OutflowImbalanceTerm`'s weight was
NOT implemented (confirmed structural, not a formula — see below).
`mgp_filter.jl`'s stubs and all 8 hand-coded filter modules remain untouched,
read-only oracles.

## What changed
- Added `src/examples/mgp_decay.jl`: `DecayContribution` (struct: `event`,
  `alpha`, `gate::Rational{Int}`, `value`, `reason::Symbol`),
  `decay_contribution(event, x, θ, ℓ, n) -> DecayContribution` (DEATH/SAMPLE
  only, throws `ArgumentError` otherwise), `total_decay(model, x, θ, ℓ, n) ->
  Real` (sums `.value` over every `DEATH`/`SAMPLE` event; BIRTH/MIGRATION/
  NEUTRAL are skipped, never visited).
- Added `src/examples/mgp_proposal.jl`: `abstract type ProposalStrategy end`
  and its four singleton subtypes `NaiveProposal`/`SoftProposal`/
  `GuidedProposal`/`HardProposal`; `driver(α::Real, π::Real) = α*π`;
  `boost(Φ::Rational{Int}, π::Real) = Φ/π`, plus a
  `boost(rt::ReducedTransition, π::Real)` convenience overload.
- Wired into `src/examples/Examples.jl`: `mgp_decay.jl` included immediately
  after `mgp_mgpaudit.jl`, `mgp_proposal.jl` immediately after that (needs
  `ReducedTransition` from `mgp_reduce.jl`, already loaded earlier), both
  before `mgp_filter.jl` (untouched).
- Added `test/kli_decay_test.jl` (43 new `@test`s) and
  `test/kli_proposal_test.jl` (18 new `@test`s), both registered in
  `test/runtests.jl` right after `mgpaudit_test.jl`.

## Files added
- `src/examples/mgp_decay.jl`
- `src/examples/mgp_proposal.jl`
- `test/kli_decay_test.jl`
- `test/kli_proposal_test.jl`
- `handoffs/M07_driver_boost_decay.md` (this file)

## Files modified
- `src/examples/Examples.jl` — two lines added (the two new includes).
- `test/runtests.jl` — two lines added (the two new test includes).
- No existing struct, macro, filter, or likelihood-bearing file touched.
  `mgp_filter.jl`, `mgp_filter_ir.jl`, `mgp_reduce.jl`, `mgp_transitions.jl`,
  `mgp_phi.jl`, and all 8 hand-coded filter modules are byte-for-byte
  unchanged (read-only oracles, cited by file:line, never edited).

## Mathematical decisions

### 0. `OutflowImbalanceTerm` — confirmed structural, not a formula
Re-read `mers_filter_suite.tex`'s "Assembled regular filter" (lines
791-830) and "Decay λ" (lines 833-844) sections directly (not just trusted
from M06's handoff). Confirmed the task's framing is correct:
- The boxed λ formula (lines 833-836) has **no branch-point/fork term at
  all**: `λ = χ_C I_C + χ_H I_H + γ_C I_C·1{I_C≤ℓ_C} + γ_H I_H·1{I_H≤ℓ_H}`.
- Line 841-844 ("No separate branch-point hazard"): the fork/branch-point
  loss "arises automatically from Eqs. 7-8 inflow against the TCC/THH share
  of the (full-birth) outflow ... Adding either to λ would double-count."
- The "Assembled regular filter" equation (lines 800-814) and the
  "driver box" (846-865) confirm the mechanism is purely structural: births
  keep the FULL, unreduced hazard `α_u` as the outflow term
  (`α_TCC+α_THH+α_THC+α_TCH`, no ℓ-dependent reduction), while the inflow
  side is built ONLY from the non-fork (`RegularFlow`) reduced cases
  (`I_7..I_12`). There is no separate closed-form "imbalance weight" —
  the imbalance falls out automatically once the filter equation is
  assembled with full-hazard outflow vs. non-fork-only inflow, which is
  explicitly M08's end-to-end wiring job, not M07's.
- **No formula for `OutflowImbalanceTerm`'s weight was written anywhere in
  this milestone.** `mgp_filter_ir.jl` was not touched.

### 1. Decay λ generalization (DEATH events)
For a DEATH event `u` with hazard `α_u(x,θ)` acting on deme `d=event.from`
(the only deme a DEATH event touches — `_decode_move`'s `:chop` branch sets
`into=Int[]`):
```
decay_u(x,θ,ℓ,n) = α_u(x,θ) · 1{n[d] ≤ ℓ[d]}
```
`n[d]` is post-event occupancy, `ℓ[d]` the post-event pruned count (same
M02-M06 convention throughout). Given the model invariant `ℓ[d]≤n[d]`, this
collapses to `n[d]==ℓ[d]` exactly — matching the tex's own parenthetical
(line 841: `ind_{I_C≤ℓ_C} ≡ ind_{I_C=ℓ_C}`).

**Reduces EXACTLY to the tex's MERS formula.** Applying the generic formula
to `removal_c`/`removal_h` (α=γ_C I_C / γ_H I_H, d=1/2) reproduces
`γ_C I_C·1{I_C≤ℓ_C} + γ_H I_H·1{I_H≤ℓ_H}` term-for-term.

**Hand-verified MERS instance** (γ_c=1/2, γ_h=2/5, χ_c=1/10, χ_h=3/10, all
exact `Rational{Int}`; I_C=5,ℓ_C=2 above threshold, I_H=3,ℓ_H=3 exactly at
threshold):
```
removal_c:   α=γ_c·I_C=5/2, gate=0 (5>2)     -> value=0
removal_h:   α=γ_h·I_H=6/5, gate=1 (3=3)     -> value=6/5
sampling_c:  α=χ_c·I_C=1/2, gate=1 (SAMPLE)  -> value=1/2
sampling_h:  α=χ_h·I_H=9/10, gate=1 (SAMPLE) -> value=9/10
total_decay(MERS,...) = 13/5   (matches the tex formula exactly, since the
                                 tex formula IS what was implemented)
```
Boundary tests: `I_d==ℓ_d` (gate=1) and `I_d==ℓ_d+1` (gate=0, just above
threshold) both verified in `test/kli_decay_test.jl`.

### 2. Sampling hazard (SAMPLE events): unconditional, `r`-independent
For a SAMPLE event, its FULL hazard contributes to λ unconditionally
(`gate=1//1` always), regardless of `event.r`. Confirmed against both naive
filters' actual code, not just argued from the tex:
- MERS `sampling_c`/`sampling_h` (`r=(0,0)`, `mers_naive.jl:159`):
  `chi_c*Ic + chi_h*Ih`, unconditional.
- SEIR `sampling` (`r=(0,1)` — M03's finding, `mgp.jl:197`):
  `seir_naive.jl:120`'s decay expression includes `ψ*I`, unconditional —
  EXACTLY `PhyloPOMP.SEIR`'s `sampling` event's hazard `θ.ψ*x.I`, evaluated
  and added unconditionally, the identical shape as MERS's `chi_d*I_d`
  terms despite `r=(0,1)` vs `r=(0,0)`.

### 3. SEIR sampling `r=(0,1)` vs MERS `sample_remove` `r=(0,0)` — resolved
The `r`-vector difference (SEIR: non-destructive sample, host stays,
genuine production slot for the inserted leaf; MERS: destructive
`sample_remove`, host consumed, no production slot) affects only what
happens **at the singular (observed) event time** — a saturation/coloring
question, already out of `full_transitions`'s scope by `event.type`, not by
`r` (M03's finding). **It does not change the shape of the λ contribution**:
both are SAMPLE-type, `event.regular==false` marks whose decay contribution
is the full hazard, unconditionally, independent of `r`. Confirmed
numerically (§2 above), not just argued.

**One genuine gap found (reported, not silently absorbed)**:
`seir_naive.jl:120`'s decay expression is `ψ*I + χ*I + ...` — a SECOND rate
constant `χ` (destructive sample, `seir_naive.jl:54-64`'s `k==2` branch)
with **no corresponding `Event`** in `PhyloPOMP.SEIR`'s `MGPModel`. `@mgp
SEIR` (`mgp.jl:188-198`) declares `χ` as a model parameter but no `@event`
ever uses it in a `rate=`. `χ=0.0` in every shipped test config
(`seir_funs.jl:182,190`), so this never numerically fires, but
`total_decay(SEIR,...)` cannot reproduce the `χ*I` term even in principle —
a Population-IR representation gap (a jump type the naive filter simulates
but `@mgp` never declared), not a decay-formula disagreement (the `ψ*I`
term this milestone's code DOES compute matches `seir_naive.jl` exactly).
Fixing it would mean adding a second SAMPLE-type `Event` to `mgp.jl` —
out of scope (read-only per this milestone's constraints), flagged for later.

### 4. A discrepancy found against the naive filters' actual code (flagged, not silently resolved)
`mers_naive.jl:159-161` and `seir_naive.jl:120` do **not** literally compute
`γ_d·I_d·1{I_d≤ℓ_d}` for DEATH marks. Given the models' own invariant
`I_d≥ℓ_d`, their actual expression
`γ_d·ℓ_d + indicator(I_d≤ℓ_d, γ_d·(I_d-ℓ_d))` algebraically reduces to
`γ_d·ℓ_d` **unconditionally** (the indicator term is identically zero: given
the invariant, `I_d≤ℓ_d` forces `I_d==ℓ_d`, at which point `I_d-ℓ_d==0`
regardless). This equals this milestone's formula only at the boundary
`I_d==ℓ_d`; it **exceeds** it by exactly `γ_d·ℓ_d` whenever `I_d>ℓ_d` (where
this milestone's indicator gives 0 but the naive code still adds `γ_d·ℓ_d`).

Concrete numbers (MERS instance above): tex-formula total λ=13/5; the naive
filter's actual `decay` variable at the same state computes 18/5 (difference
=1=exactly `γ_c·ℓ_C`, the removal_c excess). Verified by hand and confirmed
in the SEIR instance too (same pattern, `seir_naive.jl:120`'s
`γ*ellI` unconditional term).

**This milestone deliberately implements the tex's literal boxed formula**
(matching `mgp_filter.jl`'s own `kli_decay` stub docstring, which cites the
identical `1_{I≤ellI}` form), not the naive filters' code, because that IS
the formula the milestone was asked to generalize, and because — a
documented, reasoned hypothesis, not a proven resolution — the naive
filters are single-PARTICLE (one simulated trajectory) importance samplers,
not a direct implementation of the PDE-level "Assembled regular filter";
the PDE's matched inflow/outflow cancellation (a real neighboring population
state `x+e_C` to borrow probability mass from) has no single-trajectory
analog, so a particle's importance weight may need to decay continuously at
the extra rate `γ_d·ℓ_d` regardless of `I_d` vs `ℓ_d`, to correctly
condition against "a currently-tracked lineage died of an unobserved
background death" — a real possibility whether or not an untracked,
compatible removal is ALSO possible. **Flagged explicitly for advisor/M08
attention** — M08's end-to-end numerical comparison against the naive
filters' total log-likelihood is the concrete gate that would settle this
either way. Full derivation and citations in `mgp_decay.jl`'s header
comment.

### 5. Boost/driver composition (Part B)
Trivial by design: `driver(α,π)=α·π`, `boost(Φ,π)=Φ/π`, matching
`mers_filter_suite.tex` lines 427-432 and `mgp_filter.jl`'s `apply_move!`
docstring ("boost B = Ψᵤ = ϕᵤ/πᵤ"). The substantive finding is in the
**numerical cross-check**, not the (intentionally trivial) functions
themselves — see below.

## Architecture decisions
- `DecayContribution` carries `event`/`alpha`/`gate`/`value`/`reason` —
  provenance-preserving, matching M04-M06's structured-value convention
  (`ReducedTransition`, `FilterTerm`), so a future `mgpaudit` extension can
  report `λ` per-event without redoing the derivation.
- `decay_contribution` takes `x, θ` explicitly (not just `ℓ, n` as the
  milestone's own pseudocode signature suggested) because `event.hazard` is
  `(x,θ)->Float64` — evaluating `α_u` genuinely requires the population
  state and parameters, not just lineage counts. Documented deviation from
  the literal suggested signature, per the milestone's own "name flexible"
  allowance.
- `total_decay` skips BIRTH/MIGRATION/NEUTRAL events entirely (never even
  calls `decay_contribution` on them) rather than having `decay_contribution`
  return zero for them — matches the established "fail loud on
  structurally-wrong calls" convention (M02/M03): a caller that accidentally
  asks for a BIRTH event's decay contribution gets an `ArgumentError`, not a
  silent 0.
- `boost(rt::ReducedTransition, π)` convenience overload added (not in the
  literal spec) so a future M08 caller can boost a `ReducedTransition`
  directly, mirroring the `event`/`kli_binomial_ratio` two-arity pattern M02
  established.
- `mgp_proposal.jl` is placed after `mgp_decay.jl`, both before
  `mgp_filter.jl`, and after `mgp_reduce.jl` (needed for the
  `ReducedTransition` overload's type).

## Numerical cross-check (Part B, required)
Target: `seir_naive.jl`'s `infection` event, `k==2` branch (tracked-parent,
CrossDeme I→E), same concrete instance M04/M05 verified
(`reduced_transitions(infection,[2,2],[6,5])` gives `Φ(cross)=1/10`).

Traced `seir_naive.jl:111,115-116` (`pi[2]=ellI/I`, pre-event `ellI`),
`:145` (shared `ll -= decay*step+log(pi[k])`), and `:150-156` (the k==2
branch: `ll += log(ellI)`; `swap!` reassigns `ellI` to its post-event value;
`ll += log(1-ellI/I)-log(E)`). With `a=3` (pre-event `ellI`, so post-event
`ellI=2`, matching `ℓ_I=2` in the M04 instance), `I=5`, `E_post=6`:
```
naive_factor = (1/pi[2]) · a · (1-ellI_post/I) · (1/E_post)
             = (5/3)·3·(3/5)·(1/6) = 1/2                      (exact Rational)
```
Separately: the FULLY-marginalized proposal probability of the reduced
("cross") outcome is `π_u(cross) = pi[2]·q = (a/I)·(1/a) = 1/I = 1/5`
(`q=1/a` is `seir_naive.jl:152`'s `rand(cols[Infec])` per-lineage choice,
which algebraically cancels `a` — a general, parameter-independent fact:
`π_u(cross)=1/n_I` always for this branch). Then:
```
boost(1//10, 1//5) == 1//2  ==  naive_factor      (exact match, confirmed in code)
```
This confirms `π_u` in the `boost(Φ,π)` formula must be the FULLY
proposal-marginalized probability of the reduced outcome (`pi[k]·q`
combined), not just the top-level categorical `pi[k]` a hand-coded kernel
happens to name `pi` — documented explicitly in `boost`'s docstring.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before (M06 baseline): **4153/4153 passed**.
- After this milestone: **4214/4214 passed** (`4153 + 43 + 18`), ~44s wall
  time, 0 failures, 0 errors. The delta is exactly the two new testsets:
  `"Decay lambda (M07 Part A)"` (43 tests, `test/kli_decay_test.jl`) and
  `"Proposal backend (M07 Part B)"` (18 tests, `test/kli_proposal_test.jl`).
  No pre-existing testset's pass count changed.

## Verification status
| Layer | Status | Oracle |
|---|---|---|
| `OutflowImbalanceTerm` weight — confirmed structural, no formula written | Confirmed (re-read tex directly) | `mers_filter_suite.tex` lines 791-844 |
| Decay λ, DEATH generalization | Verified (new), exact `Rational{Int}` | `test/kli_decay_test.jl`, hand-computed against tex's boxed formula |
| Decay λ, SAMPLE generalization | Verified (new) | same, MERS + SEIR |
| MERS cross-check vs `mers_filter_suite.tex`'s formula | Exact match (by construction) | `test/kli_decay_test.jl` "MERS concrete instance" |
| MERS cross-check vs `mers_naive.jl`'s actual `decay` code | Discrepancy found and documented (not silently resolved) | `mgp_decay.jl` header §"A DISCREPANCY FOUND"; this handoff §4 |
| SEIR cross-check vs `seir_naive.jl`'s actual `decay` code | Same discrepancy pattern found (DEATH); `ψ*I` term matches exactly (SAMPLE); `χ*I` representation gap found | this handoff §3-4 |
| Boundary `I_d==ℓ_d` / `I_d==ℓ_d+1` | Verified (new) | `test/kli_decay_test.jl` "MERS boundary" |
| `ProposalStrategy` taxonomy | Implemented (vocabulary only, no `π_u` derivation) | `test/kli_proposal_test.jl` |
| `driver`/`boost` generic composition | Verified (new), trivial by design | `test/kli_proposal_test.jl` |
| Boost numerical cross-check vs `seir_naive.jl` k==2 | Verified (new), exact match | `test/kli_proposal_test.jl`, this handoff |
| `π_u` derivation for soft/guided/hard | NOT implemented — explicitly out of scope | N/A |
| End-to-end filter wiring | NOT implemented — explicitly out of scope, M08 | N/A |

## Known failures / unresolved issues
- **Flagged for advisor**: the naive-filter-vs-tex-formula discrepancy for
  DEATH-mark decay (§4 above) is unresolved — a documented hypothesis
  (particle-level SMC weight decay ⊇ PDE-level λ), not a proven fact. This
  milestone implemented the tex's literal formula (matching
  `mgp_filter.jl`'s own `kli_decay` docstring), not the naive filters' code.
  If the advisor's reading differs, `mgp_decay.jl`'s DEATH branch is the
  single place to revise — it does not touch SAMPLE (which matches exactly)
  or Part B.
- **Flagged for advisor**: a genuine SEIR Population-IR representation gap
  found (§3) — `seir_naive.jl` simulates a second, `χ`-rated destructive
  sample jump type with no corresponding `Event` in `mgp.jl`'s `SEIR`
  model. Currently harmless (`χ=0` in all shipped configs) but structurally
  real. Out of scope to fix here (would require editing `mgp.jl`, a
  read-only oracle per this milestone's constraints).
- No test-authoring bugs found during this milestone (unlike M02/M03/M05,
  which each found and fixed one) — all hand-derivations matched on first
  computation, cross-checked against Julia itself before being written into
  the permanent test file.

## Git state
- branch: `atpabuser-devel`
- commit if created: none
- uncommitted files (cumulative, M00-M07; nothing committed or pushed by any
  milestone so far, per instructions):
  - `docs/compiler/architecture_before.md`, `docs/compiler/compiler_roadmap.md` (M00)
  - `handoffs/M00_reconnaissance.md` .. `handoffs/M06_audit_provenance.md` (M00-M06)
  - `src/examples/mgp_audit.jl`, `test/population_ir_test.jl` (M01)
  - `src/examples/mgp_phi.jl`, `test/kli_phi_test.jl` (M02)
  - `src/examples/mgp_transitions.jl`, `test/kli_full_transitions_test.jl` (M03)
  - `src/examples/mgp_reduce.jl`, `test/kli_reduce_test.jl` (M04)
  - `src/examples/mgp_filter_ir.jl`, `test/kli_filter_ir_test.jl` (M05, M06 rename)
  - `src/examples/mgp_explain.jl`, `src/examples/mgp_mgpaudit.jl`, `test/mgpaudit_test.jl` (M06)
  - `handoffs/M07_driver_boost_decay.md` (new, this file)
  - `src/examples/mgp_decay.jl` (new)
  - `src/examples/mgp_proposal.jl` (new)
  - `test/kli_decay_test.jl` (new)
  - `test/kli_proposal_test.jl` (new)
  - `src/examples/Examples.jl` (modified — two more `include` lines added)
  - `test/runtests.jl` (modified — two more `include` lines added)

## Resume instructions
1. Advisor reviews `src/examples/mgp_decay.jl` and `src/examples/mgp_proposal.jl`,
   in particular: (a) the DEATH decay discrepancy against the naive filters'
   actual code (§4) — confirm or correct the "particle vs. PDE" hypothesis;
   (b) the SEIR `χ`-rate representation gap (§3); (c) whether `π_u` in
   `boost(Φ,π)` should be documented even more explicitly as "the fully
   proposal-marginalized probability", given how easy it would be for a
   future caller to pass just a top-level `pi[k]` instead.
2. M08 ("End-to-End SEIRS/SEIR Compiler") should treat everything M01-M07
   built as ready to wire together into `mgp_filter.jl`'s actual
   `kli_select`/`kli_decay`/`apply_move!`/`singular_update!` stubs: `kli_decay`
   is now `total_decay` (this milestone) plus the still-unresolved
   inflow/outflow-imbalance algebra for the fork mass (`OutflowImbalanceTerm`,
   confirmed structural in this milestone, not yet computed); `apply_move!`
   is `boost`/`driver` (this milestone) composed with the still-undeprived
   `π_u` from the naive/soft/guided/hard proposal families (M07's
   `ProposalStrategy` tags them, does not derive them). M08 is the first
   major project gate: assemble the full regular/singular filter equation and
   compare its log-likelihood against the trusted hand-coded `seir_naive.jl`/
   `mers_naive.jl`, over many states/parameters/seeds, per the master plan's
   5 verification gates (structural, full KLI, reduced KLI, filter structure,
   numerical equivalence) — the first 4 are substantially built by M01-M07;
   #5 (numerical equivalence) is M08's actual, still-unattempted payload, and
   is exactly the gate that would resolve this milestone's flagged DEATH-decay
   discrepancy empirically.
3. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any M08
   change and confirm the count does not regress below 4214/4214 (this
   milestone's new baseline).

## Next milestone
M08 — End-to-End SEIRS/SEIR Compiler: the first major project gate. Wire
`mgp_filter.jl`'s `kli_select`/`kli_decay`/`apply_move!`/`singular_update!`
stubs using M01 (Population IR), M02 (φ_u), M03 (full transitions), M04
(Φ_u reduction), M05/M06 (Filter IR + provenance/audit), and M07 (decay λ +
driver/boost composition + proposal-strategy tags) as the building blocks,
producing an actual executable filter. Then run the master plan's 5th
verification gate — numerical equivalence — comparing the assembled filter's
log-likelihood against `seir_naive.jl`/`mers_naive.jl` over many
states/parameters/seeds. This is also the concrete test that would settle
M07's flagged, unresolved DEATH-decay discrepancy (§4 above) empirically:
if the assembled filter using this milestone's tex-literal λ formula matches
the naive filters' total log-likelihood, the "particle vs. PDE" hypothesis
is confirmed (the discrepancy nets out via some other term M08 assembles);
if it does not match, `mgp_decay.jl`'s DEATH branch needs revision toward
the naive filters' actual per-particle formula instead.

## Context note
The most important thing to preserve if this conversation were compacted:
`src/examples/mgp_decay.jl`'s `decay_contribution`/`total_decay` implement
KLI's λ(t,x,y) generically for the first time, generalizing
`mers_filter_suite.tex`'s MERS-specific boxed formula
(`λ = χ_C I_C+χ_H I_H + γ_C I_C·1{I_C≤ℓ_C}+γ_H I_H·1{I_H≤ℓ_H}`) to any
DEATH event (`α_u·1{n[d]≤ℓ[d]}`, `d=event.from`) and any SAMPLE event
(`α_u`, unconditional, independent of `event.r` — confirmed this does NOT
differ between MERS's `r=(0,0)` `sample_remove` and SEIR's `r=(0,1)`
`sample`). **A genuine, confirmed discrepancy was found** against the naive
filters' own running `decay` code (`mers_naive.jl:159-161`,
`seir_naive.jl:120`): their actual per-step DEATH contribution is
`γ_d·ℓ_d` UNCONDITIONALLY, not gated by `I_d≤ℓ_d` — matching this
milestone's tex-literal formula only at the `I_d==ℓ_d` boundary, exceeding
it by `γ_d·ℓ_d` whenever `I_d>ℓ_d`. This milestone implemented the tex's
literal formula anyway (matching `mgp_filter.jl`'s own `kli_decay`
docstring) and documented the discrepancy as an open, flagged question for
M08's numerical-equivalence gate to settle, rather than silently picking one
formula and hiding the disagreement. `OutflowImbalanceTerm` was
**confirmed** (by re-reading the tex directly) to be structural — resolved
by the assembled filter equation's full-hazard-outflow-vs-non-fork-inflow
imbalance, with no separate closed-form weight — and was correctly left
unimplemented, per scope. `src/examples/mgp_proposal.jl`'s
`ProposalStrategy`/`NaiveProposal`/`SoftProposal`/`GuidedProposal`/
`HardProposal` tag the existing four-kernel taxonomy (M00's decision, kept
not replaced); `driver(α,π)=α·π` and `boost(Φ,π)=Φ/π` are intentionally
trivial, but were cross-checked numerically against `seir_naive.jl`'s actual
`infection` `k==2` (tracked-parent, CrossDeme) branch and matched exactly
(`boost(1//10,1//5)==1//2`, confirmed to equal the naive filter's actual
combined log-likelihood correction there) — the key finding being that
`π_u` in this formula must be the FULLY proposal-marginalized probability
(`pi[k]·q` combined), not a hand-coded kernel's top-level categorical `pi[k]`
alone. Baseline is now 4214/4214 (`4153` M06 baseline `+ 43` Part A tests
`+ 18` Part B tests).
