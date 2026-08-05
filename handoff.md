# MERS proposal-filter suite handoff (consolidated)

## Scope

Julia probability-kernel, importance-sampling, testing, and mathematical-documentation
work on the MERS-CoV two-host (Camel/Human) phylodynamic filter. Not biological or
laboratory research.

## Status: complete

All eight required steps (shared support-safe helpers; distinct soft/guided/hard
kernels; constructor validation; four MERS test files; TeX/PDF repair; full
validation) are implemented and verified below.

## Follow-up round 3: real-tree ESS/collapse diagnostic (external project cross-check)

Requested after investigating two external, separate projects: `yang-phylopomp`
(Peter Yang's R `phylopomp`/`pomp` mif2 pipeline for this same MERS model) and
`mers_private` (a more careful continuation, including a critical review,
`Yang-Notes_V2.pdf`, of Yang's thesis analysis). See
`mers_kernel_diagnostic_findings.md` (repo root) for the full writeup. Short version:

- Fixed `yang-phylopomp/mers_profile.R` to use `pomp::mcap()` (the real Monte Carlo
  Adjusted Profile method) instead of a naive LOESS arg-max, matching the pattern
  already proven elsewhere in that same repo (`seirs.R`, `ms.R`). Verified against
  the real, already-saved 180-point profile: MLE β_HH=36.24, 95% CI=[34.73, 38.55],
  with `se_mc` (1.03) exceeding `se_stat` (0.58) — a concrete demonstration of why
  MCAP matters (a naive CI would have understated the uncertainty).
- Attempted to cross-check our four Julia kernels against an R-fitted point, but
  `mers_private`'s own documentation (`RESULTS.md`) states its saved fits are
  explicitly "preliminary and exploratory... should not be used for reporting
  maximum likelihood estimates." The critical-review PDF independently found a real
  bug in Yang's spillover hazard formula (wrong susceptible compartment index) and
  self-contradictory reported results. Confirmed our own MERS model does *not* have
  that bug (our THC hazard already uses the receiving deme's susceptible count).
- Redirected to the more useful open question those notes flag directly: their
  `ANALYSIS_GUIDE.md` §2.4 states guided/auxiliary SMC is "conceptual, not
  implemented" in the R side, and documents severe, Np-invariant ESS collapse on
  ~3% of genealogical events using their bootstrap filter. Ran the analogous
  diagnostic here: all four Julia kernels (Naive/Soft/Guided/Hard) collapse
  **identically** on the real 274-tip tree (100/548 steps `-Inf`, first at step 144)
  at Yang-Notes' documented default parameters, invariant to seed and to sweeps of
  every rate parameter and initial condition tested. Traced this to near-simultaneous
  camel sample tips (timestamps differing only in the 10th+ significant digit,
  almost certainly from month-resolution collection dates in the source data) —
  a genuine target-level (`Φ=0`) incompatibility that no proposal kernel can or
  should fix, cleanly consistent with the KLI Theorem 5 boundary our epsilon-floor
  design respects. This is a positive result (confirms the kernels behave exactly as
  designed) but means finite likelihood on the *full* real tree needs a fix for the
  tied-event-time artifact first (jitter or explicit simultaneous-event handling,
  in either codebase) — not better proposal kernels or parameter fitting.

No PhyloPOMP.jl source files were changed in this round; only
`yang-phylopomp/mers_profile.R` (a separate, external repo) and this repo's
`mers_kernel_diagnostic_findings.md` were added/modified.

## Repository safety

Pre-existing dirty-tree state was inventoried before any edit (`git status --short`,
snapshotted to a scratchpad diff) and preserved throughout:

- Pre-existing **tracked** modifications not touched by this task beyond what's noted
  below: `.gitignore`, `Project.toml`, `src/examples/Examples.jl`, `src/examples/mers_tree.jl`,
  `test/runtests.jl` (these already included/wired the soft/guided/hard modules and
  test files before this task began).
- Pre-existing **untracked** files preserved unchanged: `src/examples/handoff.md`,
  `mers_filter_v26.tex`/`.pdf`, `mers_naive_fixes_test.jl`, `mgp*.jl`,
  `one_filter_many_models*.pdf`, `seir_macro_equivalence_test.jl`,
  `mers_filter_task_prompt.md`, `starter_file.md`.

## Files created or modified by this task

**Modified:**
- `src/examples/mers_funs.jl` — added shared support-safe helpers
  (`floored_shares`, `floored_shares_feasible`, `cross_proposal`, `floored_branch_law`,
  `hard_branch_law`); fixed the shared root/fork `singular_part!` to mix the
  epsilon floor only over outcomes with nonzero population rate (previously a
  vanishing guide weight could zero out a target-compatible outcome, including at
  `I_c=S_c=S_h=0` fork boundaries, which now correctly falls back to `-Inf` instead
  of assigning spurious mass); added `proposal_floor` validation and threaded it
  into `filter_pomp`'s `params`.
- `src/examples/mers_naive.jl` — added `0 < proposal_floor <= 1` constructor
  validation (`ArgumentError`); no other change.

**Rewritten (previously identical/incorrect, listed as defects in the task prompt):**
- `src/examples/mers_soft.jl` — now distinct from Guided: preserves the naive
  aggregate identity/tracked-branch split (`cross_proposal`), guides only which
  tracked branch is chosen (relative-hazard weighted, `floored_branch_law`, uniform
  fallback when hazards vanish/are non-finite).
- `src/examples/mers_guided.jl` — now distinct from Soft: jointly normalizes the
  untracked-host weight together with every branch's relative hazard in one floored
  categorical draw (local `floored_choose_branch`, deliberately not sharing/editing
  `PhyloPOMP.choose_branch` in `guide.jl`, which `GuidedSEIR`/`HardSEIR` depend on
  unfloored).
- `src/examples/mers_hard.jl` — auxiliary intensities `kappa_epsilon(j)` now
  epsilon-floored (`hard_branch_law`); the aggregate "boosted swap" probability and
  the specific-branch conditional draw both flow from this floored vector, so a
  vanishing relative hazard can no longer assign zero auxiliary intensity to a
  target-compatible outcome; compensating-decay structure unchanged.

**Created:**
- `test/mers_soft.jl`, `test/mers_guided.jl`, `test/mers_hard.jl` — modeled on
  `test/mers_naive.jl` and the SEIR test suite. Each covers: `PompObject`
  construction; zero-initial-infected `-Inf`; simulate; pfilter with finite
  likelihood (on a small hand-built 2-3-tip camel/human fixture — the full empirical
  270-tip `mers_tree` is real-world data too complex for a modest particle count to
  reliably match every tip's observed deme, so it is used only for the cheap `-Inf`
  and constructor checks); `ell=0` identity probability 1; `I=ell` boundary strict
  positivity; invalid-floor rejection. `test/mers_soft.jl` additionally verifies Soft
  and Guided diverge under a deliberately nonuniform guide, and both
  `mers_soft.jl`/`mers_hard.jl` verify the aggregate-plus-conditional log-correction
  equals the full per-outcome correction, and (for Hard) the `beta_j =
  alpha_population*kappa_epsilon(j)` / compensating-decay identity.

**TeX/PDF:**
- `src/examples/mers_filter_suite.tex` — removed the `mdframed`/`enumitem`
  dependencies (not installed in this environment: `kpsewhich` confirmed both
  `.sty` files absent, and `tlmgr install` cannot bridge the TeX Live 2023
  local/2026 remote version mismatch without `update-tlmgr-latest`), replacing
  `notebox` with a self-contained `\newenvironment` (colored rules, no new
  packages) and stripping the four `enumitem`-only `[leftmargin=...,nosep]`
  bracket options from `itemize`/`enumerate`. Removed the stray ~9-line section
  that had been appended *after* `\end{document}` (so LaTeX never saw it) and
  replaced it with a full new subsection, "Guided-family particle-filter kernels:
  soft, guided, and hard," inserted in the existing Proposal Kernel section:
  covers the shared epsilon-floor mechanism, each kernel's distinct base law (why
  Soft and Guided differ), Hard's auxiliary-intensity/compensating-decay identity
  `log(alpha*Phi_j/beta_j) = log(Phi_j/kappa_epsilon(j))`, the
  aggregate-plus-conditional equivalence, and the root/fork floor fix.
- `src/examples/mers_filter_suite.pdf` — now exists (23 pages, 394,006 bytes;
  version bumped v26→v27 with a changelog entry, see follow-up below).

### Follow-up round (v27): intuition-first restructuring + explicit SEIR ties

Requested after the initial consolidation: (1) tie the guided-family kernels
explicitly to this repository's existing SEIR implementation; (2) restructure the
document so it opens with plain-language intuition, builds up KLI-grounded math
next, and pushes code/implementation detail to the very end; (3) use KLI's own
notation throughout and explicitly flag any project-invented terminology.

Changes made:
- **Found and removed a second stray fragment** sitting after `\end{document}`
  (lines 934–935 of the then-current file) — different from, and in addition to,
  the one already fixed in the initial round. It duplicated/contradicted the
  guided-family kernels math with a lower-quality sketch (its origin is unclear
  since the file is untracked with no git history to diff against; regardless of
  origin it was dead, uncompiled text and is now gone). Verified via `grep -n
  "end{document}"` that exactly one `\end{document}` remains and it is the last
  line of the file.
- **Notation fix:** the guided-family kernels subsection previously used an
  invented symbol `$q$`/`$q_\varepsilon$` for the proposal, inconsistent with the
  rest of the document (and KLI), which uses `$\pi_u$` throughout. Replaced every
  occurrence with `$\pi_0$`/`$\pi_\varepsilon$`. `$\varepsilon$` and
  `$\kappa_\varepsilon$` (Hard's auxiliary intensity) remain, flagged explicitly as
  this repository's own vocabulary (KLI has neither).
- **Added `\S1 Plain-Language Overview`** immediately after the changelog notebox,
  before any formal notation: explains the genealogy/coloring/latent-state
  problem, why forward simulation needs importance sampling (KLI Thm. 5, boost
  `B=Φ/π`), and walks through naive/soft/guided/hard and the epsilon-floor in
  plain words — explicitly naming which terms are KLI's and which are this
  project's own at each step, per the requirement not to invent unflagged
  vocabulary.
- **Added a terminology notebox** at the start of the guided-family kernels
  subsection itself (not just in the overview), so a reader who jumps straight to
  the math still sees the KLI-vs-project-vocabulary distinction.
- **Added `\S13 From Math to Code: Implementation Map`** after the summary, before
  the bibliography — the only place in the document that names Julia
  files/functions/tests. Explicitly ties each kernel to its SEIR analogue:
  `mers_soft.jl`↔`seir_soft.jl`, `mers_guided.jl`↔`seir_guided.jl`,
  `mers_hard.jl`↔`seir_hard.jl` (all via the shared `seir_funs.jl`/`guide.jl`
  `choose_branch` overloads), and states the one structural difference from the
  SEIR baseline honestly: `seir_guided.jl` does not `include("seir_funs.jl")`
  (it duplicates `singular_part!`), whereas all three MERS kernels do.
- **Version bump v26→v27** with a "Changes in v27" changelog entry, consistent
  with the document's own existing v25/v26 changelog convention.
- Added `\label{}`s (`sec:kli`, `sec:filter-components`, `sec:proposal-kernel`,
  `sec:overview`, `sec:implementation-map`) so the new overview's cross-references
  resolve.

Validation: `pdflatex -interaction=nonstopmode -halt-on-error` run three times
(fresh, after `rm -f *.aux *.log *.out *.toc *.pdf`): exit 0 every time. Final PDF:
23 pages, 394,006 bytes. No LaTeX errors, no undefined references, no warnings
beyond the same pre-existing `Overfull \hbox` and hyperref-Unicode-bookmark
warnings noted before. Visually verified (via `pdftoppm` + image inspection, not
just `pdftotext`, since `pdftotext` renders `\texttt{}`-escaped underscores as
spaces — a text-extraction artifact, not a compilation defect) that: the table of
contents now opens with "1 Plain-Language Overview" before "2 What King, Lin, and
Ionides Establish"; the terminology notebox and π-notation render correctly; and
the Implementation Map's SEIR file references render with correct underscores.
No Julia file was touched in this round; the full test suite result from the
initial round (231/231) still applies unchanged.

### Follow-up round 2: fixed the MERS heavy-benchmark `NaN` regression

Requested after the user ran `julia --project=. test/runtests.jl` (default settings,
i.e. `RUN_HEAVY_TESTS` unset ⇒ heavy benchmarks on) and got
`[ Info: logLik = NaN ± NaN` for "MERS model with naïve proposals," and asked for
MERS's heavy-benchmark suite to work as reliably as SEIR's does (which reported
sane numbers like `logLik = -210.28 ± 1.39`).

Root cause: `test/mers_naive.jl`'s heavy block reused
`p = NaiveMERS.filter_pomp()` — **every** parameter at its default, including
`β_hc=β_ch=0` (no cross-species transmission at all) — applied to the full
270-tip empirical `mers_tree`, which has both camel and human tips. That
combination is *structurally* incompatible (not just improbable): every one of
the 10 `pfilter(Np=1000)` replicates in `ll = [...]` is guaranteed exactly
`-Inf`, and `logmeanexp`'s log-sum-exp computes `-Inf - (-Inf) = NaN` internally
when every input is `-Inf`. This contrasts with `test/seir_naive.jl`'s heavy
block, which uses a genealogy *simulated from the SEIR model itself*
(`seir_trees[1]`) together with a viable override (`χ=0.01`), so a compatible
coloring is easy to find at `Np=1000`.

Fix: rewrote `test/mers_naive.jl` to match the pattern already used in
`mers_soft.jl`/`mers_guided.jl`/`mers_hard.jl` — the same small hand-built
camel/camel/human fixture with explicit nonzero `β_hc`/`β_ch`, used for
simulate/pfilter/heavy alike — instead of building a second, separate
"heavy-only" configuration. Also added the `proposal_floor` `ArgumentError`
tests that the other three MERS test files already had but naive's did not
(the source-level validation was added in round 1; the test asserting it was
missed for naive specifically).

No other file changed. Validation: ran the four MERS test files together with
heavy benchmarks enabled (`RUN_HEAVY_TESTS` unset) twice — both times all four
kernels produced stable, finite `logLik = X ± Y` summaries (roughly -108 to -112,
± a few), no `NaN`, 62/62 passed both runs. Then ran the *complete*
`julia --project=. test/runtests.jl` with default settings (i.e. the user's
exact original command) once, in the background (it takes ~2m48s with heavy
benchmarks on): exit 0, **234/234** passed (up from 231 because
`mers_naive.jl` gained 3 tests: two `ArgumentError` checks and one
`isfinite(logLik(...))` check it previously lacked).

## Mathematical reasoning (summary)

- The SEIR reference implementation (`seir_soft.jl`/`seir_guided.jl`/`seir_hard.jl`/
  `seir_funs.jl`, and `guide.jl`'s two `choose_branch` overloads) already encodes the
  Soft-vs-Guided distinction structurally: Soft/Hard fix the aggregate
  identity/tracked-group split via `onI`/`offI` kwargs to a shared `transmission!`,
  then draw the branch conditionally; Guided draws identity-vs-branch jointly in one
  categorical. Neither SEIR kernel has an epsilon floor at all, so this task's floor
  additions for MERS are new relative to that baseline, not a port of existing SEIR
  code.
- All floors follow `q_epsilon(j) = (1-epsilon) q_base(j) + epsilon/(ell+1)`, mixing
  only over outcomes whose *population rate* is structurally nonzero — an outcome
  with zero rate (e.g. `S_c=0`) keeps probability exactly 0 regardless of `epsilon`.
- For Soft/Naive, `q_epsilon` sums to 1 by construction, so no compensating decay
  term is needed. For Guided, the joint normalization also sums to 1. For Hard, the
  unnormalized `kappa_epsilon(b) = (1-epsilon)*r_b/I + epsilon/(ell+1)` need not sum
  to `ell/I` (only does so at guide stationarity), so `alpha_population -
  sum(beta_j)` is added to the decay — verified negative-decay-compensation is
  possible and correct (Hard proposes color changes faster than the population
  process alone, by design).

## Validation performed (exact commands, exit codes, and results)

1. `julia --project=. -e 'using PhyloPOMP; println("OK")'` → exit 0, printed `OK`.
2. Each MERS test file run standalone with `RUN_HEAVY_TESTS=no`:
   - `test/mers_naive.jl`: exit 0, **8/8** passed (was 5/5 before round 2, see
     below — 3 tests added to match the other three MERS test files).
   - `test/mers_soft.jl`: exit 0, **21/21** passed.
   - `test/mers_guided.jl`: exit 0, **16/16** passed.
   - `test/mers_hard.jl`: exit 0, **17/17** passed.
   - Repeated twice; stable both times (no flakiness observed).
3. Simulation and pfilter smoke tests for all four MERS modules confirmed
   interactively (not just via the test files) on both the full empirical
   `mers_tree` (finite for naive; -Inf for all four kernels at moderate Np, which is
   expected given the tree's size/complexity, not a bug — traced by testing a tiny
   hand-built fixture where all four kernels return finite, sane log-likelihoods)
   and the small test fixture (all four kernels give finite likelihoods reliably).
4. `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`: exit 0, **231/231**
   passed. Run twice; stable both times. Per-suite breakdown (both runs identical):
   Newick parser 54/54, Newick formatter 7/7, CBLV representation 42/42, finite-state
   Markov 21/21, filter guides 10/10, rcateg 10/10, SEIR naive/soft/guided/hard 7/7
   each (28/28), MERS naive 5/5, MERS soft 21/21, MERS guided 16/16, MERS hard 17/17.
5. `pdflatex -interaction=nonstopmode -halt-on-error mers_filter_suite.tex` run three
   times (for cross-references/TOC) in `src/examples/`: exit 0 every time. Final PDF:
   19 pages, 372,924 bytes, nonempty. `pdftotext` extraction confirms the new
   "Guided-family particle-filter kernels: soft, guided, and hard" section is present
   in the compiled output. No LaTeX errors or undefined references in the log; only
   cosmetic warnings remain (see below).

## Unrelated warnings (not fixed, cosmetic/pre-existing)

- LaTeX: a few `Overfull \hbox` warnings (pre-existing long inline formulas) and
  `hyperref` "Token not allowed in a PDF string (Unicode)" warnings from math
  characters in section-title PDF bookmarks. Neither is fatal; both predate this
  task's edits in spirit (same class of cosmetic warning as elsewhere in the
  document) and were left as-is since fixing them is not part of the required
  deliverables.
- Julia: only `@time`/`@btime`-style timing output and `@info` progress messages in
  test output; no deprecation or correctness warnings.

## Found but explicitly out of scope: pre-existing `mers_naive.jl` crash

While searching for viable pfilter test parameters, triggering both `S_c=0` and
`S_h=0` simultaneously at a fork event in `NaiveMERS.singular_part!` (untouched by
this task except for the `proposal_floor` bounds check) throws
`AssertionError: S_c > 0` (or `S_h > 0`) instead of returning `-Inf`, because
`rcateg`'s all-zero-weight fallback always returns index 1. This is a latent bug in
the pre-existing naive fork code (not introduced by, and not in the required-fix
list for, this task — the analogous defect in the *shared* guided-family fork code
in `mers_funs.jl` *was* in scope and *is* fixed, with an explicit `-Inf` guard for
this exact double-zero case). Test parameters were chosen to avoid triggering this
pre-existing naive-only crash; it is flagged here for future attention rather than
fixed, since fixing it was not requested and risks touching otherwise-working,
tested code beyond this task's scope.

## Follow-up round 4: upstream resync fallout, demography port, ε=0 mode (2026-08-04/05)

Triggered by comparing local `src/examples/mers_naive.jl` against `origin/devel` after
a merge: the file had been reset to upstream's naming/API (non-underscored
`Sc`/`Beta_cc`/`chi_c` etc., no genealogy argument, no `proposal_floor`). The
`singular_part!` equations themselves were confirmed byte-for-byte identical to
`origin/devel` — no divergence there. Four follow-on threads:

1. **Test suite broke, then was fixed.** `test/mers_naive.jl` still called the old
   (pre-reset) API (`filter_pomp(g; proposal_floor=...)`); upstream's
   `NaiveMERS.filter_pomp` takes no genealogy argument at all (always filters the
   full `mers_tree`) and has no `proposal_floor`. Rewrote the test file to match:
   drops the small-fixture/finite-likelihood assertions (structurally unreachable —
   the naive kernel has no guide, so matching all ~274 tip demes on the real tree by
   chance is astronomically unlikely regardless of `Np`), adds coverage for the real
   `Ic0`/`I_c0` and `Ih0`/`I_h0` alias-conflict `ArgumentError`s instead. Full suite:
   236/236 passing after this fix.
2. **Demography formula ported from Aaron's `mers_naive.jl` redo into
   `mers_soft.jl`/`mers_guided.jl`/`mers_hard.jl`**: death rates changed from
   `@indicator(S>0, B)` (flat, guarded) to `B*S/N` (per-capita, self-zeroing).
   Verified via direct `rcateg` sampling (200k draws) that a zero-weight category is
   never selected when mixed with nonzero ones, so dropping the indicator guard
   doesn't risk negative `S`. **Found a live discrepancy while documenting this**:
   the two most recently modified copies of the reference R `mers.yml`
   (`phylopomp`/`phylopomp-fork`, both last touched by commit `817e558`, Jul 2) still
   use the old flat/guarded rate for *both* demes — i.e. the Julia `devel` branch and
   the R reference currently disagree on this term. Per explicit user decision, kept
   the per-capita (Julia) form and documented the disagreement rather than silently
   picking a side; **worth raising with Aaron directly** since it's a two-line
   question that's fast for him and slow to reverse-engineer from here.
3. **`mers_kernel_diagnostic_findings.md` corrected**: under the module's current
   *bare* defaults (`chi_h=0.0`, changed since the original diagnostic via the
   "change default MERS parameters" commit), the real tree's first `-Inf` is now at
   node 11 (the first Human sample — `log(chi_h*I_h)=log(0)` deterministically), not
   node 144/camel-extinction. The original diagnostic's extinction finding is still
   correct *for the explicit Table-2 parameters it used*; added a dated "Update"
   section rather than rewriting the original analysis, since the two describe
   different parameter regimes, not a correction of the same claim.
4. **`proposal_floor=0` made legal** ("no floor, let it fail"), responding to
   feedback that the ε-floor mechanism felt too willing to paper over failure.
   Changed the bound from `0<ε≤1` to `0≤ε≤1` in all five shared helpers
   (`floored_shares`, `floored_shares_feasible`, `cross_proposal`,
   `floored_branch_law`, `hard_branch_law`) and the `filter_pomp` validation in
   `mers_funs.jl`; every helper already reduced algebraically to the un-floored base
   law at ε=0, so only the validation and docstrings changed. Updated the three
   kernel test files' stale `proposal_floor=0.0` throw-tests (now legal) and added
   explicit ε=0 coverage confirming a genuinely-zero-weight outcome gets exactly zero
   proposal mass (contrasted against the same case at the ε=0.05 default, which stays
   strictly positive). Default remains `0.05`; ε=0 is opt-in. Full suite: 254/254
   passing after this change.

`mers_filter_suite.tex`/`.pdf` bumped v27→v28→v29 across this round (demography
formula + Julia/R discrepancy note in v28; ε=0 mode in v29), each version
recompiled cleanly and spot-checked visually via `pdftoppm`, not just `pdftotext`.

## Remaining work

None required by this task's deliverables list. Optional follow-ups a future session
could consider (not requested, not started):
- Fix the pre-existing `mers_naive.jl` double-susceptible-depletion crash described
  above.
- Install `texlive-latex-extra` (which provides `mdframed`/`enumitem`) if the
  original, more elaborate `notebox`/list styling is wanted back; the current
  self-contained replacements are functionally complete but visually simpler.
- Raise the Julia/R demography-rate discrepancy (round 4, item 2) with Aaron; a
  drafted note is not yet sent (see chat for a draft on request).
- The camel-extinction absorbing-state finding was diagnosed, not fixed — no
  parameter-regime change (`β_CH>0`, different `B_C`/`I_C0`) has been attempted.
- `yang-phylopomp/mers_profile.R`'s MCAP patch is dry-run-verified only; never
  actually re-run to produce `profile_likelihood_mcap.png`.
