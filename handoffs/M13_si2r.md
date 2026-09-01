# Handoff — M13: SI2R (Superspreading) Compiler Generalization Test

## Status
COMPLETE

## Objective
Answer the project's recurring circularity concern (raised explicitly in
M11, addressed partially in M12b for the Fork/Moran slice) with a second,
independent kind of evidence: run the M02-M07 compiler pipeline
(`enumerate_saturations`, `kli_binomial_ratio`, `full_transitions`,
`reduce_event_indicator`, `total_decay`) on a model that has **never
existed in this repository or been used during the compiler's
development** — SI2R, a two-deme (low-rate/high-rate spreader)
superspreading model from `/home/dislam/Desktop/Projects/efficiency/si2r/
si2r_model.qmd` — and check every result against that document's
independently hand-derived filter table (18 lines) and equations, not
against anything written by this project.

## What changed
- Added `src/examples/mgp_si2r.jl`: declares `SI2R` via `@mgp`, 9 events
  (`TL`, `TH`, `L`, `H`, `RL`, `RH`, `W`, `SL`, `SH`), matching the QMD's
  compartments (`S`,`I_L`,`I_H`,`R`), demes (`I_L`,`I_H`), and `alpha_u`
  column verbatim. `mgp.jl` already `include`s it (line 121); `SI2R` is
  exported alongside `SEIR`/`MERS` (`mgp.jl:21`).
- Added `test/kli_si2r_test.jl` (9276 `@test`s), registered in
  `test/runtests.jl`. Structure:
  - Part 0/0b: structural fields (`type`,`from`,`r`,`Δ`,`regular`) for all
    9 marks against the QMD's production column, plus all 9 marks' actual
    `hazard(x,θ)` closures against the QMD's `alpha_u` column at a
    concrete state (hazards are otherwise never exercised by Parts 1-7,
    which are all algebraic/combinatorial and hazard-free).
  - Part 1: `kli_binomial_ratio` vs. a from-scratch `qmd_binomial_ratio`
    reimplementation, swept for all 9 marks (not just the 4 BIRTH/
    MIGRATION marks), with a non-vacuity guard on
    `enumerate_saturations`'s length.
  - Part 2/2b: `full_transitions` classification (Identity/
    InlineSameDeme/CrossDeme/Fork) at a concrete state for TL/TH/L/H,
    PLUS the specific "QMD lines 1&2 collapse, lines 4&5 collapse" claim
    checked two ways: grouping structure (`reduce_event_indicator`
    produces exactly 2/3/2/2 groups, by type-set) and the QMD's weighted
    closed forms (Eq. 132/135), reproducing M09 Finding 2's
    `phi_id + C(ell,1)*phi_inline` weighting (NOT `reduce_event_indicator`'s
    raw, unweighted `Phi_noop`, which is a different, smaller quantity).
  - Part 3: Chu-Vandermonde identity, `sum phi_u(s)*C(ell,s) == 1`.
  - Part 4: Phi_u invariant + per-member phi correctness + no-loss
    conservation across `reduce_event_indicator`'s groups.
  - Part 5: `total_decay` vs. the QMD's Eq. (decay), with a non-vacuity
    guard forcing the `I<=ell` leftover branch to fire on some trials
    (it's a ~6%/trial event under pure random sampling, easy to leave
    accidentally untested).
  - Part 6/6b: TL's same-deme Fork = Kingman/Moran `1/C(I_L,2)`
    (independent of ell, mirrors M12b's MERS TCC/THH result on a THIRD
    model) vs. TH's CROSS-deme Fork = ordered-pair `1/(I_L*I_H)` — the
    two production slots sit in different demes, so there's no "which
    slot" symmetry to collapse, a genuinely different shape from TL's
    that a compiler merely pattern-matching MERS's TCC would get wrong.
  - Part 7: `kli_binomial_ratio` against the QMD's algebraically
    SIMPLIFIED boost-equation closed forms (lines 158-162: TH cross,
    L cross/noop, H cross/noop) — independent verification content, since
    Part 1's `qmd_binomial_ratio` is otherwise just a second
    implementation of the same product formula `kli_binomial_ratio`
    itself computes.

## Key structural findings (all confirmed exactly as the task predicted)
- **TL** (`r=(2,0)`, from I_L): same shape as MERS's TCC. Lines 1-3 ->
  Identity+InlineSameDeme collapse into one `:noop` group, Fork alone.
  Fork phi = `1/C(I_L,2)`, independent of `ell_L`/`I_H` -- Kingman/Moran,
  now confirmed on a THIRD model (SEIR has no same-deme-fork event; MERS
  TCC/THH were the first two).
- **TH** (`r=(1,1)`, from I_H): same shape as MERS's THC/TCH, EXCEPT its
  Fork is genuinely new: TL's Fork is unordered (`1/C(I_L,2)`, one deme),
  TH's Fork is ordered across two different demes (`1/(I_L*I_H)`), since
  MERS's THC/TCH never have `s` summing to 2 at all (each is `r=(1,1)`
  with only ONE slot ever fillable per saturation-limited-by-`r`, never
  a 2-slot Fork). SI2R's TH is the first model in this project with a
  cross-deme `r=(1,1)` BIRTH whose Fork case (`s=(1,1)`, both slots
  filled) is actually reachable and combinatorially distinct from the
  same-deme case.
- **L, H** (migration, `r=(0,1)`/`(1,0)`): same shape as SEIR's
  `progression` -- Identity + CrossDeme only, no InlineSameDeme possible
  (the "from" deme has `r_from=0`, so no slot ever exists in that deme to
  collapse with Identity).
- **RL, RH**: DEATH, no `phi_u` -- decay-only, confirmed by `total_decay`
  Part 5.
- **W**: NEUTRAL, no genealogical effect, `r=(0,0)`, confirmed.
- **SL, SH**: SAMPLE, singular only (`regular=false`), `r=(1,0)`/`(0,1)`
  derived from `move=sample(I_L)`/`sample(I_H)` (not written explicitly
  in the DSL -- Part 0's `sl.r == [1,0]` / `sh.r == [0,1]` checks pin this
  down, since nothing else in the DSL declaration states it directly).

## Verification against the QMD's actual equations, not just its table
The task flagged a real risk during review (caught by the advisor before
the first version of this test was accepted as complete): a naive
`qmd_binomial_ratio(n,ell,r,s) = prod C(n-ell,r-s)/C(n,r)` helper is
mathematically the SAME formula `kli_binomial_ratio` implements, so a
sweep comparing the two is two implementations of one product formula
agreeing -- real, but weaker evidence than it first appears. Two things
were added specifically to get independent content:
1. Part 7 checks against the QMD's algebraically SIMPLIFIED closed forms
   (its boost-equation lines 158-162), not the product form.
2. Part 2b's weighted-closed-form check reproduces the QMD's OWN stated
   simplification path (Eq. 132/135: `phi_id + ell*phi_inline` collapsing
   the other deme's factor to exactly 1 via the binomial-sum identity)
   and confirms the *algebra*, not just that both sides evaluate to the
   same product.

## Deviation, documented (not a bug)
`N` (total population) is declared as a free parameter in `SI2R`'s
`params=` tuple, matching `mgp.jl`'s `SEIR` precedent (`θ.N`), even though
the QMD's own parameter table (si2r_model.qmd, lines 61-69) omits `N`
entirely (implicitly `S+I_L+I_H+R`). This follows the existing SEIR/MERS
convention in this codebase rather than the QMD's textual parameter list;
flagged here for visibility, not treated as a discrepancy to fix.

## Tests run
`RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`:
- Before this milestone (M12b baseline, confirmed by re-running): **21266/21266**.
- After: **30542/30542** (`21266 + 9276`), ~50s wall time, 0 failures, 0 errors.

## Scope: what this validates, and what it does not
This validates the M02-M07 pipeline's **generality** across a genuinely
new model shape (cross-deme `r=(1,1)` fork, a case no prior model in this
repo reaches) against an externally-authored reference document. It does
**not** validate:
- Any filter/proposal/driver machinery for SI2R -- no `mgp_si2r_filter.jl`
  was written (not requested; the task's scope was the M02-M07 IR pipeline
  and the QMD's filter-DESIGN table, not an executable compiled filter).
  SI2R therefore has no naive/soft/guided/hard filter, no simulator
  cross-check, and no Gate-5-style numerical-equivalence test, unlike
  SEIR/MERS.
- The QMD's driver assembly (Eq. driver) or boost (Eq. boost) equations in
  full -- Part 7 checks selected `phi_u` VALUES the boost equation states
  as intermediate quantities, not the boost/driver formulas themselves as
  assembled expressions.
- Anything about parameter regimes or dynamics beyond the states directly
  swept (`I_L`,`I_H` up to ~200, `ell` up to `I`).

## Files added
- `src/examples/mgp_si2r.jl`
- `test/kli_si2r_test.jl`
- `handoffs/M13_si2r.md` (this file)

## Files modified
- `src/examples/mgp.jl` -- one line (`include("mgp_si2r.jl")`), `SI2R`
  added to the `export` list.
- `test/runtests.jl` -- one line (`include("kli_si2r_test.jl")`).

No `seir_naive.jl`/`mers_naive.jl`/hand-coded filter files touched.

## Next milestone
M10's flagged question (does BDEI need a genuine new `@mgp` model, or does
SEIR's `progression` already cover the "recoloring not collapsible with
identity" stress test?) remains open and unrelated to this milestone.
Separately, a natural follow-up here would be a `mgp_si2r_filter.jl` +
`si2r_pomp.c`-style naive filter (the QMD's own Task 3) to get a Gate-5
numerical-equivalence check for SI2R the way M08/M09 did for SEIR/MERS --
not attempted here since the task's explicit scope was the IR pipeline
against the QMD's filter-design table, not an executable filter.
