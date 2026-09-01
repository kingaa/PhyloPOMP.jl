# Handoff — M11: Model Cross-Validation / Additional Known Models (BDSS)

## Status
COMPLETE

## Objective
Per the master project plan's own instruction for this milestone: "Do not
force tests for models that do not exist in the repository or for which no
reliable oracle is available. Clearly label: exact oracle / handwritten
implementation comparison / property test only / not yet independently
verified." M09's "Next milestone" section and M10's own re-check both
already confirmed BDSS does not exist in this repository. This milestone
re-confirms that (see "Verification" below) and produces the required
labeling table honestly, from what M00-M10 actually established — not
copied blindly from the task prompt's illustrative example.

## Re-confirmation: BDSS does not exist in this repository
`grep -rli "bdss\|birth.death.sampling.sampling\|bd-ss"` across every
`.jl`/`.tex`/`.md` file returns zero matches. `grep -rn "@mgp "` across
`src/` returns exactly two model definitions: `@mgp SEIR` (`src/examples/mgp_macro.jl:188`)
and `@mgp MERS` (`src/examples/mgp_mers.jl:3`). No third model exists at any
milestone. Per the master plan's own explicit instruction quoted above, no
test is forced for BDSS.

## The labeling table

```
exact oracle:                            none
handwritten implementation comparison:   SEIR, MERS
property test only:                      none (see note below)
not yet independently verified:          BDEI, BDSS  (neither model exists in this repository)
```

### Why "exact oracle" is empty, not SEIR/MERS

The task prompt's own illustrative example asks specifically whether
anything rises above "handwritten implementation comparison" to "exact
oracle" status — a genuinely independent, non-project-authored source of
truth (e.g. a closed-form analytical result). Re-checking M00's findings
directly (`handoffs/M00_reconnaissance.md`'s "Known failures" section,
"Not investigated in this milestone"):

- `mers_filter_suite.tex`, the hand-derivation document that M02-M04's
  generic `enumerate_saturations`/`kli_binomial_ratio`/`full_transitions`/
  `reduce_event_indicator` were checked against term-by-term, is itself
  **project-authored** — a hand-worked derivation written for this project,
  not an independent external source. M00 explicitly flagged that its
  formulas were checked for internal self-consistency with the codebase,
  not independently re-derived from the external KLI paper
  (`StructuredMGPs.pdf`, which lives outside the repo and was never read in
  depth by any milestone).
- `seir_naive.jl`/`mers_naive.jl` (the Gate-5 comparison targets for M08/M09)
  are hand-coded, single-particle importance-sampling filters written for
  this same project, independently of the M02-M09 compiler pipeline in the
  sense that neither consumes the other's code — but both are still
  project-authored implementations of the same KLI theory, not an external
  ground truth. This is precisely why the task labels this tier
  "handwritten implementation comparison", distinct from "exact oracle".
- The one genuinely external, independent check in this project is the
  forward simulator's cross-validation against R's `phylopomps::runSEIR`
  (`scripts/seir_crossvalidate.{jl,R}`, `results_7.23.md`, cited in M00's
  verification-status table). This validates the forward Gillespie
  simulator's demographic/event-timing correctness — it does NOT touch the
  KLI genealogy filter, φ_u/Φ_u math, or log-likelihood computation, which
  is this project's actual object of study. It cannot be cited as an "exact
  oracle" for the filter math specifically.
- No closed-form analytical likelihood result (e.g. an exactly-solvable
  toy case with a known analytic log-likelihood) exists anywhere in this
  repository for either SEIR or MERS.

Conclusion, stated plainly rather than inflated: **no model in this
repository currently has "exact oracle" status.** The strongest evidence
tier actually achieved is "handwritten implementation comparison"
(Gate 5 numerical equivalence, exact bit-for-bit `Rational{Int}`-derived
matches), for SEIR and MERS only.

### SEIR, MERS — handwritten implementation comparison

Both models pass every one of the master plan's 5 verification gates
(structural, full-KLI, reduced-KLI, filter structure, numerical
equivalence), most concretely:
- SEIR: `test/kli_seir_compiled_test.jl`, Gate 5 — 289/289 finite
  log-likelihood comparisons against `seir_naive.jl` matched exactly
  (permanent test), 1913/1913 in a broader development sweep (M08).
- MERS: `test/kli_mers_compiled_test.jl`, Gate 5 — 692/692 finite
  log-likelihood comparisons against `mers_naive.jl` matched exactly
  (permanent test), 2903/2903 in a broader development sweep (M09).

Both comparisons are against a hand-coded, independently-written (not
compiler-generated) reference implementation of the same underlying KLI
theory — the "handwritten implementation comparison" tier, exactly as
defined by the master plan, not "exact oracle" (see above) and not merely
"property test only" (the comparisons are exact numerical matches against a
second, independent implementation, not self-consistency checks on the
compiler's own output).

### "property test only" — none, and why that's not a gap

M12 (this project's next milestone, executed immediately after this one)
adds broad randomized property tests over SEIR and MERS
(`test/kli_properties_test.jl`) — but since both models already carry the
strictly stronger "handwritten implementation comparison" evidence, the
property tests are additive verification on top of an already-verified
model, not a model's ONLY evidence. No model in this repository sits in the
"property test only" tier (that tier is reserved for a model that has
randomized structural checks but no independent handwritten comparison to
cross-check numerical correctness against — that does not describe SEIR or
MERS here).

### BDEI, BDSS — not yet independently verified

Neither model exists anywhere in this repository:
- BDEI: confirmed absent by M09's own grep and re-confirmed structurally by
  M10 (`handoffs/M10_bdei.md`) — no `@mgp BDEI` block, no `bdei_naive.jl`,
  no test file, at any milestone. M10 substituted a synthesis exercise over
  SEIR's existing `progression`/birth events as an ACCEPTABLE partial
  substitute for the master plan's specific "migration recoloring not
  collapsible with identity" stress test (an advisor-approved scope
  decision, not a claim that a BDEI model itself was built or verified).
- BDSS: confirmed absent by this milestone's own re-check above (zero grep
  matches, only two `@mgp` blocks in the entire codebase).

Both are labeled "not yet independently verified" rather than any stronger
tier — there is no BDEI or BDSS model, hand-coded reference, or test to cite
evidence from. This is a factual statement about what does not exist, not a
statement that the compiler pipeline itself would fail on such a model; M10
argued (without building a new model) that the compiler's generic
machinery — `full_transitions`/`reduce_event_indicator`, which are provably
event-agnostic and demography-agnostic per M04's "Architecture decisions" —
does not special-case SEIR/MERS in a way that would obviously break on a
differently-shaped model, but this is a design argument, not verification,
and is not represented in the table above.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before this milestone: **4270/4270 passed** (M10's unchanged baseline,
  confirmed live).
- After this milestone: **4270/4270 passed**, unchanged — no code was
  added or modified (this milestone is pure labeling/bookkeeping, per its
  own scope: "Keep this short — it's a labeling/bookkeeping milestone, not
  new derivation work").

## Known failures / unresolved issues
- None newly introduced. Carried over from M09/M10, unresolved by this
  milestone (out of scope): whether to fix `mgp_mers.jl`'s `death_c`/
  `death_h` hazard-formula mismatch further (M09 Finding 3 already
  resolved the specific mismatch flagged at that time — see the pre-M00
  git history's `931af3b`/subsequent fix noted in this project's task
  preamble); whether `mgp_filter.jl`'s generic stubs should ever be filled
  in generically.

## Git state
- branch: `atpabuser-devel`
- commit: none (nothing committed or pushed by any milestone so far, per
  instructions)
- Files changed by this milestone: `handoffs/M11_model_crossvalidation.md`
  (new, this file) only. `git status` is otherwise identical to M10's end
  state.

## Resume instructions
1. Advisor reviews the labeling table above, in particular the "exact
   oracle: none" conclusion — confirm or push back on whether the R
   `phylopomps::runSEIR` cross-validation (forward-simulator-only) or
   `mers_filter_suite.tex` (project-authored) should be treated
   differently than argued here.
2. Proceed to M12 (already executed alongside this milestone in the same
   dispatch — see `handoffs/M12_property_testing.md`).

## Next milestone
M12 — Property-Based Verification: broad randomized invariant tests over
SEIR/MERS's KLI compiler machinery, generalizing the targeted-instance
pattern M02-M09 already used into explicit, named, hundreds-of-trials
PROPERTY tests. See `handoffs/M12_property_testing.md` (executed in the
same dispatch as this milestone).

## Context note
The most important thing to preserve if this conversation were compacted:
this milestone made **no code changes** — it is pure bookkeeping, populating
the master plan's own required oracle-labeling scheme honestly from what
M00-M10 actually established. The key, potentially surprising conclusion:
**no model in this repository has "exact oracle" status** — even SEIR/MERS's
strong Gate-5 numerical-equivalence results are "handwritten implementation
comparison" (two independent project-authored implementations agreeing),
not verification against a genuinely external, non-project-authored source
of truth (the external KLI paper itself was never read in depth by any
milestone; the one true external check, R's `phylopomps::runSEIR`, only
validates the forward simulator, not the filter/likelihood math). BDEI and
BDSS remain "not yet independently verified" because neither model exists
anywhere in this codebase — re-confirmed by a fresh grep in this milestone,
not merely copied from M09/M10's prior findings. Baseline unchanged at
4270/4270.
