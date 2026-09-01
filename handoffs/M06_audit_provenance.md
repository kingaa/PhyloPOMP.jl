# Handoff — M06: Audit Provenance (Part 0 Correction, `explain`, `@mgpaudit`)

## Status
COMPLETE

## Objective
Two parts. **Part 0 (small, first):** fix a mislabeling the human supervisor
found in M05's `DecayTerm` — the `:fork` bucket M05 built is NOT KLI's
`λ(t,x,y)` decay term (a different mechanism per `mers_filter_suite.tex`),
and the name invited exactly the double-counting mistake the tex warns
against. **Part 1:** implement `explain(term)`, an ordinary function walking
the full provenance chain `FilterTerm -> ReducedTransition -> KLITransition
-> Event` M01-M05 built, following M01's `audit_model`/`ModelAudit`
structured-value-plus-`Base.show` convention. **Part 2:** implement
`mgpaudit`/`@mgpaudit`, an end-to-end per-mark derivation report walking
`audit_model` (M01) through `classify_filter_terms` (M05/M06) for every
event of a model, with DEATH/SAMPLE/NEUTRAL events explicitly (never
silently) marked out of scope.

## What changed

### Part 0 — Correction to M05
Renamed `DecayTerm` (`src/examples/mgp_filter_ir.jl`) to
`OutflowImbalanceTerm`, and `FilterSpec.decay` to
`FilterSpec.outflow_imbalance`. Added a required `mechanism::Symbol` field
(always `:inflow_outflow_imbalance` in this and M05's milestone; the
docstring is explicit that a hypothetical `:lambda` value must never be
produced by this type — see "Mathematical decisions" below). Updated
`test/kli_filter_ir_test.jl`'s every reference (imports, field accesses,
comments) and added one new testset asserting the mechanism-vs-λ distinction
structurally (`"OutflowImbalanceTerm mechanism is not lambda (M06
correction)"`, 11 new `@test`s: mechanism is always
`:inflow_outflow_imbalance`, never `:lambda`; `FilterSpec` has no field
literally named `decay`; `mgp_filter.jl`'s real λ stub, `kli_decay`, is a
structurally separate, still-unimplemented code path).

### Part 1 — `explain`
Added `src/examples/mgp_explain.jl`: `explain(t::KLITransition) ->
TransitionExplanation`, `explain(rt::ReducedTransition) ->
ReducedExplanation`, `explain(term::RegularFlow|SingularFlow|
OutflowImbalanceTerm) -> TermExplanation`. Each returned struct has a
`Base.show(io, MIME"text/plain", ...)` method (M01's convention). Walks the
full chain down to the source `Event` (mark name, semantic `EventType`,
production vector `r`), the saturation `s` and exact `φ_u` of every
contributing full transition, the full KLI transition kind
(`Identity`/`InlineSameDeme`/`CrossDeme`/`Fork`), which full transitions
collapsed into a `ReducedTransition` and its `Φ_u`, and (for a `FilterTerm`)
the `regular_flow`/`singular_flow`/`outflow_imbalance` classification bucket
plus `mechanism`/`reason` where applicable.

### Part 2 — `mgpaudit` / `@mgpaudit`
Added `src/examples/mgp_mgpaudit.jl` (a SEPARATE file from `mgp_audit.jl` —
see "Architecture decisions" for why): `default_audit_state(model) ->
(ℓ, n)` (one illustrative default state per model, documented); `mgpaudit(model;
ℓ=nothing, n=nothing) -> MgpAuditReport`, an ordinary function walking
`audit_model` (M01) → `full_transitions` (M03, which internally calls M02's
`enumerate_saturations`/`kli_binomial_ratio`) → `reduce_event_indicator`
(M04) → `classify_filter_terms` (M05/M06) → `explain` (Part 1) for every
BIRTH/MIGRATION event, and an explicit `MarkAuditOutOfScope` note (reusing
`full_transitions`'s/`ChopTransition`'s established out-of-scope reasoning)
for every DEATH/SAMPLE/NEUTRAL event — never a silent omission. `@mgpaudit
ModelName [key=value ...]` is a THIN macro forwarding its arguments to
`mgpaudit`, containing no derivation logic of its own (per the project's
macro-hygiene rule).

## Files added
- `src/examples/mgp_explain.jl` (Part 1)
- `src/examples/mgp_mgpaudit.jl` (Part 2)
- `test/mgpaudit_test.jl`
- `handoffs/M06_audit_provenance.md` (this file)

## Files modified
- `src/examples/mgp_filter_ir.jl` — Part 0 rename (`DecayTerm` ->
  `OutflowImbalanceTerm`, `FilterSpec.decay` -> `.outflow_imbalance`, new
  `mechanism` field), header/docstrings rewritten to cite the exact tex
  lines the human supervisor pointed at.
- `test/kli_filter_ir_test.jl` — every `DecayTerm`/`spec.decay` reference
  updated to `OutflowImbalanceTerm`/`spec.outflow_imbalance`; one new
  testset added (see Part 0 above).
- `src/examples/mgp_audit.jl` — header comment updated to point to
  `mgp_mgpaudit.jl` for `@mgpaudit`/`mgpaudit` (the M01 `audit_model`/
  `validate_model` code itself is BYTE-FOR-BYTE unchanged; only the file's
  own header prose was edited).
- `src/examples/Examples.jl` — two lines added: `include("mgp_explain.jl")`
  (after `mgp_filter_ir.jl`) and `include("mgp_mgpaudit.jl")` (after
  `mgp_explain.jl`), both before `mgp_filter.jl` (untouched).
- `test/runtests.jl` — one line added: `include("mgpaudit_test.jl")`.
- `mgp_filter.jl` and all 8 hand-coded filter modules untouched (read-only
  oracles, never edited, per every milestone's standing constraint).

## Mathematical decisions

### Correction to M05 (Part 0) — full writeup
The human supervisor independently re-read `mers_filter_suite.tex` lines
750-865 and found M05's `DecayTerm` name misleading. The exact tex lines
(re-verified against the actual file, not just the milestone-spec's
paraphrase — line numbers below are `grep -n` confirmed against
`src/examples/mers_filter_suite.tex` as checked out for this milestone):

- **Line 760-762** (end of the TCC derivation, Eq. 7): *"The bracketed
  factor is the reduced non-fork contribution for TCC. The missing mass
  corresponds to unobserved within-camel branch points and is accounted for
  through the regular birth inflow/outflow balance (§8), not through λ; the
  same holds for THH (Eq. 8)."*
- **Lines 833-836** (the "Decay λ" section, Sec. "Filter Components"):
  ```
  λ(t,x,y) = [sampling hazard: χ_C I_C + χ_H I_H]
           + [sub-threshold removal: γ_C I_C·1{I_C≤ℓ_C} + γ_H I_H·1{I_H≤ℓ_H}]
  ```
  — note this formula has NO branch-point/fork term in it at all.
- **Lines 841-844** ("No separate branch-point hazard"): *"The within-deme
  unseen branch-point loss arises automatically from Eqs. 7-8 inflow against
  the TCC/THH share of the (full-birth) outflow; the cross-deme loss
  likewise from Eqs. 9-12 against the THC/TCH share. Adding either to λ
  would double-count."*

So the fork/branch-point mass is resolved by a THIRD mechanism — the
structural imbalance between the full-hazard regular outflow (`α_u`,
unreduced — lines 846-865's "driver box" keeps births FULL in the reduced
regular exit rate) and the non-fork-only regular inflow (Eqs. 7-12, i.e.
M05's/M06's `RegularFlow`) — which is textually and mathematically DISTINCT
from KLI's `λ` (lines 833-836's formula). M05's `DecayTerm` name implied
this bucket IS (or directly feeds) `λ`; it does not, and the tex explicitly
warns that treating it as such double-counts.

**The fix**: renamed `DecayTerm` → `OutflowImbalanceTerm`
(`src/examples/mgp_filter_ir.jl`), and `FilterSpec.decay` →
`FilterSpec.outflow_imbalance`. Chose the RENAME (over keeping `DecayTerm`
with just a `mechanism` tag) because the milestone spec's own framing —
*"the end state MUST make it structurally awkward... for a future milestone
to accidentally sum this bucket into a λ total without noticing"* — is best
served by removing the string "decay" from both the type name and the field
name entirely: a future M07 implementer grepping for "decay" or literally
writing `sum(spec.decay)` will find NOTHING to accidentally sum, and must
instead discover `outflow_imbalance` and read why it is named that way. As
an additional, independent guard (not a replacement for the rename), also
added a REQUIRED `mechanism::Symbol` field to `OutflowImbalanceTerm`, fixed
at `:inflow_outflow_imbalance` for every instance this milestone
constructs, with a docstring stating a hypothetical `:lambda` value must
NEVER be produced by this type — if a future milestone needs to represent
KLI's actual `λ` pool, it should be a genuinely different type (so
`mechanism == :inflow_outflow_imbalance` remains a true, checkable
invariant of `OutflowImbalanceTerm` specifically, not a discriminated union
pretending to unify two different quantities). `reason` (`:fork_
unobserved_at_regular_time`) is unchanged from M05 — it answers "why was
this excluded from RegularFlow", which the rename doesn't affect.

This is a NAMING/TYPING fix only — no new numeric derivation. `Φ` on
`OutflowImbalanceTerm` is still the raw, uncomputed `ReducedTransition.Φ`
copied verbatim (M05's original design decision, unchanged); working out the
actual algebraic resolution of the inflow/outflow imbalance remains
EXPLICITLY OUT OF SCOPE, deferred to M07 (see "Next milestone" below).

### `explain`/`mgpaudit` (Parts 1-2): no new KLI math
Both parts are pure presentation/composition layers over the already-tested
M01-M05(+M06 Part 0) pipeline — no new saturation, φ_u, Φ_u, classification,
or decay/λ logic was written. `mgpaudit`'s only genuine judgment call is
`default_audit_state` (see "Architecture decisions").

## Architecture decisions
- **`mgp_mgpaudit.jl` is a SEPARATE file from `mgp_audit.jl`**, even though
  the milestone spec's natural-home suggestion was to extend `mgp_audit.jl`
  directly. Reason: `MarkAuditInScope` (Part 2) has fields typed
  `Vector{ReducedExplanation}`/`FilterSpec`/`Vector{TermExplanation}`, all
  defined by files (`mgp_reduce.jl`, `mgp_filter_ir.jl`, `mgp_explain.jl`)
  included in `Examples.jl` AFTER `mgp_audit.jl` (which M01 placed
  immediately after `mgp.jl`, before any IR/math file, since `audit_model`/
  `validate_model` genuinely have no such dependency). Julia requires a
  struct's field TYPES to already exist at struct-DEFINITION (parse/eval)
  time — unlike an ordinary function body, which only needs its referenced
  globals to exist by CALL time — so a struct with those field types cannot
  live in a file included before `mgp_explain.jl`. Rather than reordering
  `mgp_audit.jl` itself (which would be an unforced, invasive change to a
  file M01 placed deliberately early, and would have no benefit since
  `mgpaudit` is genuinely a LATER-layer concept depending on M02-M05 IR),
  `mgp_mgpaudit.jl` is a new file included near the end of the chain (after
  `mgp_explain.jl`, before `mgp_filter.jl`). `mgp_audit.jl` itself is
  otherwise BYTE-FOR-BYTE unchanged (only its header comment now points
  readers to `mgp_mgpaudit.jl`).
- **`explain` return types are separate structs per chain layer**
  (`TransitionExplanation`/`ReducedExplanation`/`TermExplanation`), not one
  polymorphic "explanation" type — mirrors the existing
  `KLITransition`/`ReducedTransition`/`FilterTerm` layering (M03/M04/M05)
  rather than collapsing it. Each has its own `Base.show`, so
  `explain`-ing any layer in isolation (e.g. a bare `KLITransition` a future
  caller has on hand, without a wrapping `ReducedTransition`) still prints
  something useful, per the milestone's "ideally also for a bare
  ReducedTransition or KLITransition if that's easy" ask.
- **`mgpaudit`'s default `(ℓ, n)` state is picked per model, not per
  event**, per the milestone's own instruction ("pick a reasonable default
  state per model"). Chose SEIR's default (`ℓ=[2,2], n=[6,5]`) to be the
  EXACT instance M04/M05 already hand-verified for `infection` (reused, not
  invented); chose MERS's default (`ℓ=[2,2], n=[5,5]`) to be large enough in
  BOTH demes that all four BIRTH marks' `Fork` saturations (same-deme for
  TCC/THH, needing `ℓ_d≥2`; cross-deme for THC/TCH, needing `ℓ_C,ℓ_H≥1`) are
  simultaneously reachable, so the default report illustrates all three
  Filter IR buckets for every in-scope MERS mark rather than showing an
  empty `OutflowImbalanceTerm` list purely as an artifact of an
  under-powered default state. Verified this reasoning empirically (see
  "Tests run" — the default report's printed output was inspected and
  confirmed to show non-empty `RegularFlow`+`OutflowImbalanceTerm` for
  `infection` and for all four MERS transmission marks).
- **`mgpaudit`'s out-of-scope note is a single shared closure
  (`OUT_OF_SCOPE_NOTE`)**, not a per-event-type hand-written string, since
  M03's own scope-boundary reasoning (`full_transitions`'s `ArgumentError`,
  `ChopTransition`'s docstring) is already generic over DEATH/SAMPLE/NEUTRAL
  — reusing that reasoning verbatim (per the milestone's explicit
  instruction to reuse the established framing) rather than re-deriving it
  per event type.
- **`@mgpaudit`'s keyword-forwarding is validated at macro-expansion time**
  (each trailing argument must be an `Expr(:(=), ...)`) so a malformed call
  like `@mgpaudit SEIR ℓ` (missing `=value`) fails fast with a clear
  `ArgumentError` at expansion time rather than producing a confusing
  downstream `MethodError` from `mgpaudit` itself — this is argument-shape
  validation, not derivation logic, so it does not violate the
  "macro contains no derivation logic" rule.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before Part 0 (M05 baseline): **3983/3983 passed**.
- After Part 0 (rename + new mechanism-guard testset): **3994/3994 passed**
  (`3983 + 11`), 0 failures, 0 errors. The delta is exactly the new
  `"OutflowImbalanceTerm mechanism is not lambda (M06 correction)"` testset
  in `test/kli_filter_ir_test.jl`. No pre-existing testset's pass count
  changed (confirming the rename touched only names, not behavior).
- After Parts 1-2 (this milestone's full delivery): **4153/4153 passed**
  (`3994 + 159`), ~44s wall time, 0 failures, 0 errors. The delta is exactly
  the new `"Provenance / explain and mgpaudit (M06)"` testset in
  `test/mgpaudit_test.jl` (159 `@test`s across 8 sub-testsets). No
  pre-existing testset's pass count changed.
- Manually verified (not just via `@test`) at the REPL: `mgpaudit(SEIR)` and
  `@mgpaudit MERS`'s printed (`Base.show`) output — eyeballed to cover every
  event by name, with BIRTH/MIGRATION marks showing a non-trivial reduced-
  transition table and RegularFlow/OutflowImbalanceTerm split, and
  DEATH/SAMPLE/NEUTRAL marks showing the explicit out-of-scope note text.

## Verification status
| Layer | Status | Oracle |
|---|---|---|
| Part 0: `DecayTerm` -> `OutflowImbalanceTerm` rename, no behavior change | Verified (new) | `test/kli_filter_ir_test.jl`, all pre-existing testsets' pass counts unchanged |
| Part 0: `mechanism` field always `:inflow_outflow_imbalance`, never `:lambda` | Verified (new) | `test/kli_filter_ir_test.jl` `"OutflowImbalanceTerm mechanism is not lambda"`; `test/mgpaudit_test.jl` `"explain: OutflowImbalanceTerm"` |
| Part 1: `explain(::KLITransition)` | Verified (new) | `test/mgpaudit_test.jl` `"explain: KLITransition / ReducedTransition"` |
| Part 1: `explain(::ReducedTransition)` | Verified (new) | same testset |
| Part 1: `explain(::RegularFlow)` correctly IDs event/kind | Verified (new) | `"explain: RegularFlow (SEIR infection)"` |
| Part 1: `explain(::OutflowImbalanceTerm)` correctly IDs event/kind/mechanism | Verified (new) | `"explain: OutflowImbalanceTerm (MERS transmission_cc)"` |
| Part 1: `explain(::SingularFlow)` | Implemented, exercised only structurally (no real singular BIRTH/MIGRATION event exists in SEIR/MERS — same M05-documented gap) | not separately `@test`ed against a real event; code path shares `explain`'s generic dispatch, would be exercised by M05's synthetic-`Event` test if extended |
| Part 2: `mgpaudit(SEIR)`/`@mgpaudit SEIR` run without error, cover all 5 events | Verified (new) | `"mgpaudit(SEIR) completeness"` |
| Part 2: `mgpaudit(MERS)`/`@mgpaudit MERS` run without error, cover all 12 events | Verified (new) | `"mgpaudit(MERS) completeness"` |
| Part 2: DEATH/SAMPLE/NEUTRAL -> explicit out-of-scope note, never omitted | Verified (new) | both completeness testsets, per-event-type assertions |
| Part 2: `@mgpaudit` keyword forwarding (`ℓ=...`, `n=...`) | Verified (new) | `"@mgpaudit with explicit state"` |
| Part 2: `mgpaudit` argument validation (shape, `ℓ<=n`) | Verified (new) | `"mgpaudit argument validation"` |
| `default_audit_state` matches documented per-model choices | Verified (new) | `"default_audit_state"` |
| KLI's actual `λ`/decay weight math | NOT implemented — explicitly out of scope, deferred further (M07) | N/A |
| Inflow/outflow-imbalance resolution (the actual algebra) | NOT implemented — this milestone only renamed the placeholder so it can't be mistaken for `λ`; the algebra itself is still M07's job | N/A |
| Proposal logic (`π_u`) | NOT implemented — explicitly out of scope | N/A |

## Known failures / unresolved issues
- None outstanding in the delivered code.
- Carried forward from M05, unchanged by this milestone (not re-litigated,
  since M06's Part 0 was a naming/typing fix only, not a re-derivation): no
  real singular (`event.regular == false`) BIRTH/MIGRATION event exists in
  either shipped model, so `SingularFlow`'s classification rule (and, in
  this milestone, `explain(::SingularFlow)`) remains validated only against
  M05's synthetic hand-built `Event`, not a real one. Flagging again for the
  advisor per M05's own carried-forward note.
- New, minor, flagged for the advisor: `mgpaudit`'s default `(ℓ, n)` choice
  is a judgment call (see "Architecture decisions") — reasonable and
  documented, but arbitrary in the sense that any other `ℓ≥2` state in every
  fork-capable deme would have served the "illustrate all three buckets"
  goal equally well. Not a correctness question (the derivation itself is
  exactly as tested throughout M02-M05 at whatever `(ℓ,n)` is supplied), but
  worth a second look if a future milestone wants `mgpaudit`'s default
  output to be canonical/citable rather than illustrative.

## Git state
- branch: `atpabuser-devel`
- commit if created: none
- uncommitted files (cumulative, M00-M06; nothing committed or pushed by any
  milestone so far, per instructions):
  - `docs/compiler/architecture_before.md`, `docs/compiler/compiler_roadmap.md` (M00)
  - `handoffs/M00_reconnaissance.md` .. `handoffs/M05_filter_ir.md` (M00-M05)
  - `src/examples/mgp_audit.jl` (M01; header comment updated by M06, logic
    byte-for-byte unchanged), `test/population_ir_test.jl` (M01)
  - `src/examples/mgp_phi.jl` (M02), `test/kli_phi_test.jl` (M02)
  - `src/examples/mgp_transitions.jl` (M03), `test/kli_full_transitions_test.jl` (M03)
  - `src/examples/mgp_reduce.jl` (M04), `test/kli_reduce_test.jl` (M04)
  - `src/examples/mgp_filter_ir.jl` (M05; Part 0 rename applied by M06),
    `test/kli_filter_ir_test.jl` (M05; updated by M06's Part 0)
  - `handoffs/M06_audit_provenance.md` (new, this file)
  - `src/examples/mgp_explain.jl` (new, M06 Part 1)
  - `src/examples/mgp_mgpaudit.jl` (new, M06 Part 2)
  - `test/mgpaudit_test.jl` (new)
  - `src/examples/Examples.jl` (modified — two more `include` lines added)
  - `test/runtests.jl` (modified — one more `include` line added)

## Resume instructions
1. Advisor reviews: (a) the Part 0 rename itself
   (`src/examples/mgp_filter_ir.jl`) against the exact tex line citations
   above — confirm the rename resolves the mislabeling without overclaiming
   what the new bucket IS (it still carries no computed weight, per M05's
   original, unchanged design); (b) `explain`'s output shape
   (`src/examples/mgp_explain.jl`) — confirm it is the right shape for a
   future milestone (M07+) to build richer diagnostics on top of, rather
   than needing to be redesigned; (c) `mgpaudit`'s `default_audit_state`
   choice and the `mgp_audit.jl`/`mgp_mgpaudit.jl` file split — confirm the
   split's justification (Julia's struct-field-type-at-definition-time
   requirement) is accepted as a good reason, not a workaround that should
   instead have reordered `Examples.jl`'s earlier includes.
2. M07 ("Proposal Backend and Generic Driver/Boost/Decay") should treat
   `OutflowImbalanceTerm` (this milestone's corrected name) as the thing to
   FINALLY resolve numerically: derive the actual inflow/outflow-imbalance
   algebra (`mers_filter_suite.tex`'s `I_7`-`I_12` vs. full-hazard outflow,
   lines 735-865) that resolves this bucket's mass, and separately, derive
   KLI's actual `λ(t,x,y)` (Eq. 47/B2; the sampling hazard + sub-threshold
   removal decay, lines 833-836) as a GENUINELY DIFFERENT computation/type —
   the whole point of this milestone's Part 0 rename is to make it
   structurally awkward to conflate the two, so M07 should introduce
   whatever new type/field it needs for `λ` itself rather than reusing
   `OutflowImbalanceTerm`'s `Φ` or `mechanism` field for it.
3. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any M07
   change and confirm the count does not regress below 4153/4153 (this
   milestone's new baseline).

## Next milestone
M07 — Proposal Backend and Generic Driver/Boost/Decay: implement `λ` itself
(the actual Eq. 47/B2 decay weight — sampling hazard + sub-threshold
removal, `mers_filter_suite.tex` lines 833-836), the actual inflow/outflow-
imbalance resolution for the fork mass this milestone's Part 0 carefully
kept DISTINCT from `λ` (`OutflowImbalanceTerm`, `mechanism =
:inflow_outflow_imbalance`), and formalize the existing naive/soft/guided/
hard proposal taxonomy (`seir_naive.jl`/`seir_soft.jl`/`seir_guided.jl`/
`seir_hard.jl` and their MERS mirrors) as KLI's `π_u`. `mgp_filter.jl`'s
`kli_decay`/`kli_select`/`apply_move!`/`singular_update!` stubs (still
untouched, still `error(...)`-only) are the concrete slots this work fills
in; `explain`/`mgpaudit` (this milestone's Parts 1-2) should be extended,
not redesigned, to also report `λ`/`π_u` once they exist, following the same
structured-value-plus-`Base.show` convention.

## Context note
The most important thing to preserve if this conversation were compacted:
M05's `DecayTerm` was RENAMED to `OutflowImbalanceTerm` in this milestone
(Part 0), with `FilterSpec.decay` becoming `FilterSpec.outflow_imbalance`
and a new required `mechanism::Symbol` field (`:inflow_outflow_imbalance`,
the only value ever constructed) — this bucket is NOT KLI's `λ` and never
was; the RENAME (not just a comment) is the load-bearing fix, per
`mers_filter_suite.tex` lines 760-762 and 841-844 (`"Adding either to λ
would double-count"`). `src/examples/mgp_explain.jl`'s `explain(x)` now
walks the full provenance chain `FilterTerm -> ReducedTransition ->
KLITransition -> Event` for any node in that chain, returning a structured,
`Base.show`-able value (`TermExplanation`/`ReducedExplanation`/
`TransitionExplanation`), following M01's `audit_model` convention.
`src/examples/mgp_mgpaudit.jl`'s `mgpaudit(model; ℓ=..., n=...) ->
MgpAuditReport` and the thin `@mgpaudit ModelName [key=value...]` macro
wrapper walk the FULL M01-M06 pipeline for every event of a model at one
concrete `(ℓ, n)` state, with BIRTH/MIGRATION events getting a full
derivation and DEATH/SAMPLE/NEUTRAL events getting an explicit (never
silent) out-of-scope note. `mgp_mgpaudit.jl` is a SEPARATE file from
`mgp_audit.jl` (M01) purely because of Julia's struct-field-type-at-
definition-time ordering requirement, included near the end of
`Examples.jl`'s chain (after `mgp_explain.jl`). Baseline is now 4153/4153
(`3983` M05 baseline `+ 11` Part 0 tests `+ 159` Part 1/2 tests). No KLI
math (λ, the imbalance algebra, or `π_u`) was implemented anywhere in this
milestone — that is still M07's job.
