# Handoff — M00: Reconnaissance

## Status
COMPLETE

## Objective
Map what actually exists in PhyloPOMP.jl today — the model DSL, the
`MGPModel`/`Event` scaffold, the hand-coded SEIR/MERS filters, tests, and
KLI-adjacent source material — so that later milestones can plan a
compiler (Model DSL -> Population IR -> Full Genealogical/KLI IR ->
Reduction over m -> Reduced KLI IR -> Filter IR -> {code, audit, LaTeX,
verification}) against ground truth rather than assumption. Read-only
investigation; no likelihood math or existing source touched.

## What changed
- Explored `src/`, `test/`, `src/examples/`, root-level `.md` docs,
  `mers_filter_suite.tex`, and located the (out-of-repo) KLI source PDF.
- Ran the full test suite (`RUN_HEAVY_TESTS=no julia --project=.
  test/runtests.jl`) as a read-only sanity check — no files modified by
  this run.
- Wrote `docs/compiler/architecture_before.md`: file:line-referenced map
  of the DSL -> model representation -> filter execution path as it
  actually exists, including explicit note of which pieces are stubs.
- Wrote `docs/compiler/compiler_roadmap.md`: maps each desired compiler
  pass (parse_model, validate_model, lower_population,
  enumerate_saturations, derive_full_compatibility, compute_phi,
  reduce_event_indicator, build_filter_ir, build_proposal, compile_filter,
  audit, verify) onto exists/partial/not-present with file references.
- Wrote this handoff.

## Files added
- `docs/compiler/architecture_before.md`
- `docs/compiler/compiler_roadmap.md`
- `handoffs/M00_reconnaissance.md`

## Files modified
- None. This milestone is read-only recon; `git status --short` before
  writing the three new files showed a clean tree (see Git state below).

## Mathematical decisions
- None yet (M00 is recon only). One pre-existing convention worth
  recording: the repository already distinguishes *static* per-event data
  (`Event.Δ`, `.hazard`, `.r`, `.from`/`.into`, `.regular`, `.observed` —
  `src/examples/mgp.jl:59-69`) from *dynamic* genealogy-dependent
  quantities (saturation, ℓ, the coloring operator) explicitly in the
  `Event` docstring, matching the static/dynamic split described in this
  project's background. This convention should be preserved, not
  reinvented, in M01+.

## Architecture decisions
- Confirmed the MGP compiler scaffold (`src/examples/mgp.jl`,
  `mgp_macro.jl`, `mgp_filter.jl`, `mgp_mers.jl`) is a separate,
  currently-non-functional-past-step-1 code path, loaded *after* all eight
  hand-coded SEIR/MERS filter modules in `src/examples/Examples.jl` — it
  does not replace or feed into the working filters yet.
- Confirmed the reduced coloring `Y=(d,m)` -> `d`-only conflation flagged
  as a risk in the task background is real: `src/coloring.jl`'s
  `Coloring{D,N}` tracks only deme membership; no code anywhere in the
  repository represents the full coloring including the event-indicator
  mark `m`. All eight hand-coded filters and the MGP scaffold alike work
  directly in the reduced space, with the marginalization over `m`
  performed on paper (in `mers_filter_suite.tex`) rather than in code.
- Confirmed `compute_phi`/`Φ_u` is hand-derived and hand-transcribed per
  model per event in all eight filter modules (`seir_naive.jl`,
  `mers_naive.jl`, and their `_funs.jl`/`_soft.jl`/`_guided.jl`/
  `_hard.jl` siblings) — no generic φ_u/Φ_u computation exists. This is
  the single largest gap relative to the target compiler architecture; see
  `compiler_roadmap.md` items 4-7.
- Confirmed the epsilon-floor proposal hack flagged as a risk in the task
  background has already been fully removed (commit `931af3b`, current
  HEAD) — `grep -rln "floored\|proposal_floor\|epsilon" src/ test/`
  returns only LaTeX build artifacts and one stale descriptive comment.
  No target-math floor hack exists at HEAD, and none ever touched Φ_u
  specifically (it was always a π-only, importance-kernel-side mechanism).
- Noted a milestone-numbering collision: `check_milestone{1,2,3}.md` at
  the repo root document a *prior, unrelated* milestone sequence (the
  forward simulator work). This project's M00/M01/... is a fresh sequence
  living under `handoffs/`; future milestones should keep using
  `handoffs/M0N_*.md` to avoid confusion with the root-level files.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Result: **3248/3248 passed**, ~39s wall time, no failures, no errors.
  Read-only run (working tree unmodified by test execution). Note: no test
  in this suite exercises `mgp_filter.jl`'s `regular_step!`/`kli_select`/
  `kli_decay`/`apply_move!`/`singular_update!` — confirmed via grep across
  `test/*.jl` (no reference to any of those five names) — so this pass
  count validates the DSL/lowering/simulator/hand-coded-filter layers, not
  the MGP filter scaffold, which cannot currently run to completion.

## Verification status
| Layer | Status | Oracle |
|---|---|---|
| Structure (this milestone's own output) | N/A | N/A — docs only |
| DSL -> Event/MGPModel lowering | Verified (pre-existing) | `test/seir_macro_equivalence.jl` vs. hand-written `SEIR_REFERENCE` table |
| Forward simulator vs. hand-coded filter | Verified (pre-existing) | `test/seir_simulate.jl`, `test/mers_simulate.jl` (finiteness only) |
| Forward simulator vs. R `phylopomp::runSEIR` | Verified (pre-existing, external, not in CI) | `scripts/seir_crossvalidate.{jl,R}`, `results_7.23.md` |
| Hand-coded filter proposal-kernel agreement (naive/soft/guided/hard) | Informally verified (pre-existing) | `mers_filter_suite.tex` benchmark table; `test/mers_soft.jl` Soft-vs-Guided divergence check |
| MGP filter scaffold (`mgp_filter.jl`) | Not testable — stubs `error()` | None exists |
| Full/reduced KLI-coloring distinction (Y=(d,m) vs d) | Not represented in code at all | None exists |
| Generic φ_u/Φ_u computation | Not implemented anywhere | None exists |

## Known failures / unresolved issues
- **Design fork: where does the split between `parse_model` and
  `lower_population` go?** `@mgp`/`@event` currently fuse parsing and
  lowering into one macroexpansion pass (`mgp_macro.jl`). M01+ needs to
  decide whether to keep this fusion (simpler, but no intermediate form to
  validate/audit before lowering) or split it into a real two-stage
  pipeline. Not resolvable from the repo alone — depends on how much
  intermediate inspectability the compiler project wants.
- **Design fork: is the existing four-kernel taxonomy
  (naive/soft/guided/hard) the right `build_proposal` design, or should
  M01+ design π_u construction from scratch?** The existing taxonomy is
  real, tested, and works for two models, but it was hand-designed
  per-model rather than derived from a general procedure. A human/advisor
  call is needed on whether to formalize this taxonomy as the compiler's
  proposal-design vocabulary or treat it as one example among several
  possible proposal strategies.
- **Genuine mathematical gap, not just an engineering gap**:
  `enumerate_saturations`, `derive_full_compatibility` (Q_u),
  `compute_phi` (φ_u/Φ_u), and `reduce_event_indicator` (marginalize m)
  have *no* generic implementation anywhere in the codebase and, as far as
  this recon could determine, no in-repo document derives the *general*
  procedure — only `mers_filter_suite.tex`'s per-event, per-model
  hand-derivations exist. The external KLI paper
  (`StructuredMGPs.pdf`, located outside the repo at
  `/home/dislam/University of Michigan Dropbox/Deepan Islam/Summer-2026/
  phylodynamics_companion/Literature/StructuredMGPs.pdf`) presumably has
  the general formulas (φ_u, Φ_u, Q_u as given in the task background),
  but it was not read in depth during this recon (out of scope: recon was
  codebase-focused, per the task's Skill/tool set, and the PDF is not
  git-tracked). M01 should plan to either (a) bring the KLI paper (or the
  relevant excerpted equations) into the repo as tracked reference
  material, or (b) treat `mers_filter_suite.tex`'s hand-derivations plus
  the task background's own restated formulas as the working source of
  truth. This choice affects how "verified" any M01+ derivation can claim
  to be, and is worth an explicit advisor decision.
- **Ambiguity**: `Event.r` (production vector) is populated by the DSL but
  never read by any downstream code (`mgp_filter.jl`, `simulate.jl`) —
  confirmed by grep. It's unclear whether this is simply "not wired up
  yet" (trivial, M01+ work) or whether the current field's shape/semantics
  need to change once `compute_phi` is actually implemented against it.
  Flagging so M01 doesn't assume the field is already correct-and-just-
  unused.
- **Not investigated in this milestone** (explicitly out of scope for M00
  per the task's focus on codebase structure): the actual mathematical
  content of the KLI paper itself, and whether `mers_filter_suite.tex`'s
  per-event formulas are *correct* derivations from it (as opposed to
  merely present and self-consistent, which is what was checked). M01
  should not assume the tex document's formulas have been independently
  re-derived/verified by this recon — only that they exist, are cited,
  and pass the codebase's own internal self-consistency tests.

## Git state
- branch: `atpabuser-devel`
- HEAD commit: `931af3b` ("Remove epsilon-floor proposal machinery from
  MERS kernels")
- commit if created: none — docs not committed (per task instructions,
  M00 deliverables are left uncommitted for advisor review)
- uncommitted files:
  - `docs/compiler/architecture_before.md` (new)
  - `docs/compiler/compiler_roadmap.md` (new)
  - `handoffs/M00_reconnaissance.md` (new)
  - No existing tracked file was modified.

## Resume instructions
1. Advisor reviews the three deliverables (`docs/compiler/
   architecture_before.md`, `docs/compiler/compiler_roadmap.md`, this
   handoff) and resolves the design forks listed under "Known failures /
   unresolved issues" above, especially the `parse_model`/
   `lower_population` split and the proposal-taxonomy question.
2. Decide the KLI-paper-as-source-of-truth question (bring `StructuredMGPs
   .pdf` into the repo as tracked reference material, or excerpt the
   needed equations into a new tracked doc) before M01 starts deriving
   `enumerate_saturations`/`derive_full_compatibility`/`compute_phi`
   against it.
3. M01 ("Population IR / structural audit") should start from
   `src/examples/mgp.jl`'s existing `Event`/`MGPModel` structs (already a
   working Population IR for two models) rather than redesigning from
   scratch — the roadmap doc's item 1-3 assessment is that parsing/
   lowering/validation are close to done and mostly need restructuring
   (making validation a callable/inspectable pass; deciding on the
   parse/lower split), not new mathematics.
4. Commit the three M00 files (or let the advisor commit them) once
   reviewed — they were deliberately left uncommitted per the milestone
   instructions.

## Next milestone
M01 — Population IR / structural audit

## Context note
The most important thing to preserve if this conversation were compacted:
the codebase already has a working, tested "Population IR" in
`src/examples/mgp.jl`'s `Event`/`MGPModel` (verified against `SEIR` and
`MERS`), and a working, tested forward simulator/hand-coded-filter layer —
so M01+ is not starting from zero. What's genuinely missing is the KLI
mathematical core (saturation enumeration, Q_u, φ_u/Φ_u, marginalization
over m) as a *generic* computation; today it exists only as eight
independent, hand-derived, hand-transcribed closed-form implementations
(SEIR × {naive,soft,guided,hard}, MERS × {naive,soft,guided,hard}) plus one
hand-written LaTeX derivation document (`mers_filter_suite.tex`) that could
serve as a check target. `mgp_filter.jl`'s four stub functions
(`kli_select`, `kli_decay`, `apply_move!`, `singular_update!`) are exactly
the intended landing spots for that generic math, already docstring-
annotated with the KLI equation numbers they must satisfy — filling them in
correctly (and proving it, per the "Validation gates" checklist at
`mgp_filter.jl:190-201`) is effectively the whole compiler project's
mathematical payload. The KLI source PDF itself lives outside the git repo
on this machine's Dropbox path (see "Known failures" above) and was not
deeply read during this recon.
