# Handoff — M08: End-to-End SEIR Compiler (first major project gate)

## Status
COMPLETE

## Decay discrepancy resolution (read this first)

M07 flagged, but did not resolve, a real numerical discrepancy: `mgp_decay.jl`'s
tex-literal DEATH formula (`λ_DEATH = γ_d·I_d·1{I_d≤ℓ_d}`) does not match
`seir_naive.jl`/`mers_naive.jl`'s actual running `decay` accumulator
(`γ_d·ℓ_d`, unconditional). This milestone resolved it with a targeted,
isolated numerical experiment, then a full algebraic reconciliation,
independently re-confirmed by a broad programmatic sweep.

**Isolated experiment** (`NaiveSEIR.event_rates!` called directly, not
through the stochastic simulator, β=σ=ω=ψ=χ=0, γ=1 only):
```
I=5, ellI=2 (I>ellI):  alpha[5] (naive's regular proposal rate) = 3.0
                        naive's actual decay return           = 2.0
                        mgp_decay.jl's tex-literal λ            = 0.0   <- disagrees
I=2, ellI=2 (I==ellI): alpha[5] = 0.0
                        naive's actual decay return           = 2.0
                        mgp_decay.jl's tex-literal λ            = 2.0   <- agrees (boundary)
```
This reproduces M07's finding exactly, from the real code (script preserved
conceptually in this milestone's derivation; see `mgp_seir_filter.jl`'s
header for the full write-up).

**Resolution.** `mgp_decay.jl`'s tex-literal λ formula is **correct and
unchanged** — it is an exact transcription of `mers_filter_suite.tex`'s
boxed "Decay λ" equation. `seir_naive.jl`'s raw `decay` variable is **not**
literally λ(t,x,y): it is λ *plus* an extra term that exists only because
`seir_naive.jl`'s single-particle proposal draws regular DEATH jumps at a
**reduced** rate (`α_reduced = γ_d(n_d−ℓ_d)`, i.e. `NaiveProposal`'s "an
untracked individual recovers"), whereas the tex's own "Assembled regular
filter" / "reduced regular exit rate" box (`mers_filter_suite.tex` lines
~846-858) prescribes the regular DEATH *exit rate* as the **full**,
threshold-gated hazard `γ_d·n_d·1{n_d>ℓ_d}`. The gap between what the exact
model's regular exit rate should charge and what the reduced proposal
actually charges is un-priced probability mass that must be charged to
decay:
```
leftover_u = γ_d·n_d·1{n_d>ℓ_d} − driver(α_u, π_u^naive),   π_u^naive=(n_d-ℓ_d)/n_d
```
Algebraic identity (writing `a = n_d−ℓ_d ≥ 0`):
```
λ_tex + leftover = γ_d·n_d·(1{a=0}+1{a>0}) − γ_d·a = γ_d·n_d − γ_d·a = γ_d·ℓ_d
```
— exactly `seir_naive.jl`'s unconditional `γ_d·ℓ_d`, for **every** `n_d≥ℓ_d`,
no case-split needed once expanded. Confirmed programmatically (not just by
hand) by sweeping `I` from `ellI` to `ellI+5` at `ellI=3`, and `ellI` from 0
to 4 with `I` up to `ellI+3`, calling `NaiveSEIR.event_rates!` directly at
each point: **all matched to floating-point precision**. Also re-verified
against M07's own MERS instance: `13/5 (λ_tex) + 1 (leftover, removal_c,
I_C=5>ℓ_C=2) + 0 (leftover, removal_h, I_H=3=ℓ_H=3, boundary) = 18/5`,
M07's cited naive total.

`mgp_decay.jl` is **not modified** by this milestone (still correct as pure
λ). The new `leftover` term (and its citation of the tex's own "reduced
regular exit rate" box) lives in `src/examples/mgp_seir_filter.jl`'s
`compiled_decay`, documented in full in that file's header comment.

## Objective
Assemble M01-M07's pieces into an executable, SEIR-specific compiled filter
and run the master project's 5 verification gates, in particular Gate 5
(numerical equivalence against `seir_naive.jl`) — the first major project
gate.

## What changed
- Added `src/examples/mgp_seir_filter.jl`: `compiled_decay` (λ_tex +
  DEATH-leftover, resolving the discrepancy above), `compiled_event_rates!`
  and `compiled_regular_part!` (structural mirrors of `NaiveSEIR.event_rates!`
  / `regular_part!` — identical control flow and RNG-call order, but every
  likelihood term sourced from M01-M07's IR: `event.hazard`,
  `full_transitions`/`IdentityTransition`/`CrossDemeTransition` (M03),
  `compiled_decay`), and `compiled_filter_pomp` (mirrors
  `NaiveSEIR.filter_pomp`, reusing `NaiveSEIR.singular_part!` verbatim for
  the singular/branch-point/root/chop update — a documented scoping decision,
  see below — with `compiled_regular_part!` standing in for the regular
  update).
- Added `test/kli_seir_compiled_test.jl` (22 new `@test`s): an isolated
  `compiled_regular_part!` vs. `NaiveSEIR.regular_part!` bit-exact check
  (400 randomized (state, ℓ, params, seed) trials, mixed event types), and
  a Gate-5 end-to-end sweep (20 randomized parameter draws × genealogies ×
  150 seeds each, seed-matched `Np=1` `pfilter` runs compared).
- `src/examples/Examples.jl`: one line added (`include("mgp_seir_filter.jl")`,
  after `mgp_filter.jl`).
- `test/runtests.jl`: one line added (`include("kli_seir_compiled_test.jl")`,
  after `seir_simulate.jl`).
- **No existing file modified**: `mgp_decay.jl`, `mgp_filter.jl` (stubs
  still `error()`, untouched — see scoping decision below), and all 8
  hand-coded filter modules are byte-for-byte unchanged (`git diff --stat`
  confirms zero diff on every one of them).

## Files added
- `src/examples/mgp_seir_filter.jl`
- `test/kli_seir_compiled_test.jl`
- `handoffs/M08_seirs_end_to_end.md` (this file)

## Files modified
- `src/examples/Examples.jl` — one line (the new include).
- `test/runtests.jl` — one line (the new include).

## Scoping decision: singular part reused, not re-derived
`compiled_filter_pomp` reuses `NaiveSEIR.singular_part!` **verbatim** for
root-planting / sample-chop / branch-point-fork, rather than re-deriving it
from `full_transitions`/`reduce_event_indicator`. Justification: M02-M06
already established (Gates 2/3 below) that the IR's Φ_u numbers match
`singular_part!`'s implicit derivations term-for-term; re-deriving the
singular/fork logic here would duplicate already-verified work without new
coverage, and risks introducing bugs in code this milestone was not tasked
to re-derive. `mgp_filter.jl`'s generic stubs (`kli_select`, `kli_decay`,
`apply_move!`, `singular_update!`) are **left untouched** for the same
reason stated in the task: filling them generically for every `EventType`
requires committing to and testing a canonical dispatch signature well
beyond SEIR — explicitly out of this milestone's scope ("you do NOT need to
build a fully general `compile_filter` pass").

## A second bug found and fixed during Gate 5 development
Isolated unit testing (`compiled_regular_part!` vs. `NaiveSEIR.regular_part!`
on fixed states, no genealogy) found a **second**, independent discrepancy:
initial code computed the "identity" (k=1) contribution for infection's
regular step as `log(Φ_noop)` using `reduce_event_indicator`'s collapsed
`:noop` `ReducedTransition` — but that group sums `IdentityTransition`
**and** `InlineSameDemeTransition` together, and M03 already documented
(`test/kli_full_transitions_test.jl`'s "SEIR infection" comment) that
`InlineSameDemeTransition`/`ForkTransition` are **never reachable during a
regular step** (`seir_naive.jl`'s naive proposal only ever realizes Identity
or CrossDeme regularly — a Q_u regular/singular-gating fact M03 explicitly
scoped out). Using the collapsed group silently double-counted mass the
proposal structurally never proposes. **Fix**: use the unreduced
`IdentityTransition` directly (`only(filter(t -> t isa IdentityTransition,
ts))`), and divide out `pi_[1]` (infection's `Φ_identity` bundles both the
I-side factor, which duplicates `pi_[1]` exactly, and the E-side factor —
progression's `Φ_identity` has no such duplication since `r_E=0` makes its
E-side factor trivially 1, so no division is needed there; both cases are
documented inline in `mgp_seir_filter.jl`). Verified: 3000/3000 isolated
trials with mixed event types matched after the fix (0/3000 before).

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before (M07 baseline): **4214/4214 passed**.
- After this milestone: **4236/4236 passed** (`4214 + 22`), ~47s wall time,
  0 failures, 0 errors. Delta is exactly the new `"Compiled SEIR filter
  (M08)"` testset. No pre-existing testset's pass count changed.

## The 5 verification gates

| Gate | Result | Evidence |
|---|---|---|
| 1. Structural (Δ, α, r, event semantics, wiring) | PASS (pre-existing, confirmed still passing) | `test/seir_macro_equivalence.jl` (3006 tests), `test/population_ir_test.jl` (348 tests) — both in this run's 4236 |
| 2. Full KLI per-event transitions match trusted derivation | PASS (M03, confirmed) | `test/kli_full_transitions_test.jl` "SEIR infection"/"SEIR progression" (14 tests); re-exercised live by `compiled_regular_part!`'s `full_transitions` calls, cross-checked bit-exact against `seir_naive.jl` in `test/kli_seir_compiled_test.jl`'s isolation testset |
| 3. Reduced KLI (Φ_u) matches | PASS (M04, confirmed) | `test/kli_reduce_test.jl` "SEIR infection"/"SEIR progression" (22 tests); same live re-exercise as Gate 2 |
| 4. Filter structure (regular/singular/decay classification) + this milestone's decay resolution | PASS | M05/M06's `test/kli_filter_ir_test.jl` (unchanged, cited); this milestone's decay reconciliation (`mgp_seir_filter.jl`'s `compiled_decay`), verified via the isolated numerical experiment above and the programmatic (I,ellI) sweep |
| 5. Numerical equivalence, `log L_compiled ≈ log L_reference` | PASS | See below |

### Gate 5 detail
Approach: both filters are single-particle (`Np=1`) importance samplers
whose `ll` is itself a random variable (depends on which silent/untracked
regular events the proposal draws). `compiled_regular_part!` was built to
issue the **identical sequence** of `rcateg`/`rand()` calls, in the same
order, as `NaiveSEIR.regular_part!` — so resetting the global RNG
(`Random.seed!(seed)`) to the same value immediately before each `pfilter`
call makes the two draws bit-identical trajectories (task option (b), not a
statistical comparison).

- **Isolated unit check** (`compiled_regular_part!` vs. `NaiveSEIR.regular_part!`,
  no genealogy, fixed/randomized states): 400 trials in the permanent test
  (mixed event types, random `(S,E,I,R,ℓ,β,σ,γ,ω)` and seeds), **0
  mismatches**. An ad hoc 3000-trial sweep during development (not
  preserved as a file, numbers quoted from the development session) also
  gave 0/3000 mismatches after the two bugs above were fixed.
- **End-to-end sweep** (permanent test, `test/kli_seir_compiled_test.jl`):
  20 randomized parameter draws (β∈[1,7], σ∈[0.3,3.3], γ∈[0.3,3.3],
  ω∈[0.1,2.1], ψ∈[0.02,0.32], pop∈{60,100,150}) × a genealogy simulated per
  combo (`simulate(PhyloPOMP.SEIR,...)`, target sample count 1-4,
  tmax∈[2,6]) × 150 seeds each (3000 seed-trials attempted). Actual run:
  19/20 combos successfully built their target genealogy (1 timed out
  within the retry budget, not a filter failure), yielding **289 finite
  log-likelihood comparisons, 289 exact matches, worst `|Δll|`≈8.9e-15**
  (floating-point-noise level, not a real discrepancy). Tolerance:
  `atol=1e-6, rtol=1e-8` (chosen
  because both filters, when correct, do the same `Rational{Int}`-derived
  arithmetic converted to `Float64` at slightly different points in the
  expression — real formula bugs, e.g. both bugs found during development,
  produced differences of order 0.1-10, far outside this tolerance, so it
  is not loose enough to mask a genuine error).
- **Development-time broader sweep** (not committed as a test, to keep CI
  runtime bounded, but run and recorded here for the advisor): 25 parameter
  combos × up to 400 seeds each, **1913/1913 finite-log-likelihood pairs
  matched exactly**, spanning target sample counts 1-4 and pop∈{60,100,150}.
- One case (target sample count = 4, deeper trees) had very few finite
  (non `-Inf`) single-particle realizations (as expected: matching a longer
  observed tree with `Np=1` naive proposals is rare), but every finite pair
  found still matched exactly — no silent gap at higher sample counts, just
  fewer data points there.

## Known failures / unresolved issues
- SEIR's `χ` (destructive-sample rate) representation gap, flagged in M07,
  remains unfixed (`χ=0` in every shipped/tested config, so harmless
  numerically; fixing it requires adding a second SAMPLE-type `Event` to
  `mgp.jl`, a read-only oracle here too).
- MERS end-to-end compilation is explicitly out of scope for this milestone
  (per the task) — M09's job.
- `mgp_filter.jl`'s generic stubs remain unfilled (`error()`); this
  milestone's compiled filter lives in a separate, SEIR-specific file
  instead (see "Scoping decision" above).
- Gate 5's sweep, while broad (20 combos × 150 seeds in the permanent test,
  25 × up to 400 during development), is not exhaustive over parameter
  space; no failures were found in >2000 total finite-likelihood
  comparisons across both runs.

## Git state
- branch: `atpabuser-devel`
- commit: none (nothing committed by any milestone so far, per instructions)
- new/modified files beyond M00-M07's (see M07 handoff for the cumulative
  list up to that point):
  - `src/examples/mgp_seir_filter.jl` (new)
  - `test/kli_seir_compiled_test.jl` (new)
  - `handoffs/M08_seirs_end_to_end.md` (new, this file)
  - `src/examples/Examples.jl` (modified — one more `include` line)
  - `test/runtests.jl` (modified — one more `include` line)

## Resume instructions
1. Advisor reviews `src/examples/mgp_seir_filter.jl`'s header (decay
   resolution derivation) and the two documented bug fixes (decay leftover;
   Q_u regular/singular gating for `IdentityTransition` vs.
   `InlineSameDemeTransition`).
2. M09 (MERS end-to-end) can reuse the SAME decay-leftover pattern
   (`compiled_decay`'s DEATH branch) directly — MERS's `removal_c`/
   `removal_h` are structurally identical DEATH events to SEIR's
   `recovery`. The Q_u regular/singular-gating fix (use `IdentityTransition`/
   `CrossDemeTransition` directly, not the collapsed `:noop`/`:cross`
   `ReducedTransition` groups, for REGULAR steps) is a general finding that
   will matter for MERS's TCC/THH/THC/TCH regular dynamics too — worth
   checking explicitly in M09, not assumed.
3. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any M09
   change; must not regress below 4236/4236.

## Next milestone
M09 — End-to-end MERS compiler, applying this milestone's decay-leftover and
Q_u-gating findings to MERS's two-deme event table.

## Context note
The most important thing to preserve if this conversation were compacted:
M07's flagged DEATH-decay discrepancy is **resolved**: `mgp_decay.jl`'s
tex-literal λ is correct and unchanged; `seir_naive.jl`'s raw `decay`
variable equals λ_tex plus a `leftover = γ_d·n_d·1{n_d>ℓ_d} −
driver(α_u,(n_d-ℓ_d)/n_d)` term that exists because `NaiveProposal` uses a
reduced (not full, threshold-gated) rate for regular DEATH jumps — the
algebraic identity `λ_tex+leftover=γ_d·ℓ_d` holds for every `n_d≥ℓ_d`,
confirmed both symbolically and by a direct call into `NaiveSEIR.event_rates!`.
`src/examples/mgp_seir_filter.jl` implements this (in `compiled_decay`) plus
a full structural mirror of `seir_naive.jl`'s regular-part Gillespie loop,
sourcing each term from M01-M07's IR (`full_transitions`,
`IdentityTransition`/`CrossDemeTransition` used **directly**, not via
`reduce_event_indicator`'s collapsed groups — a second bug found and fixed
here, since the collapsed `:noop` group incorrectly includes
`InlineSameDemeTransition` mass that's structurally unreachable during a
regular/unobserved step). `NaiveSEIR.singular_part!` is reused verbatim
(documented scoping choice) rather than re-derived. Gate 5 (numerical
equivalence) passes: seed-matched `Np=1` particle-filter runs of
`compiled_filter_pomp` vs. `NaiveSEIR.filter_pomp` agree bit-for-bit
(`atol=1e-6,rtol=1e-8`) across every finite comparison in a 20-combo×150-seed
permanent test (289/289 finite pairs matched, worst |Δll|≈8.9e-15) and a
broader 25-combo×400-seed development sweep (1913/1913 matched). Baseline is
now 4236/4236 (`4214` M07 baseline `+ 22` new tests).
