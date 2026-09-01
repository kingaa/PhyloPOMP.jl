# Handoff — M12: Property-Based Verification

## Status
COMPLETE

## Objective
Per the master project plan, implement generic randomized invariant
("property") tests over the models that actually exist (SEIR, MERS),
generalizing the targeted-instance-comparison pattern M02-M09 already
started (e.g. M08/M09's 400-trial isolated `compiled_regular_part!`
checks, 289-2903-trial Gate-5 finite-comparison sweeps) into explicit,
named PROPERTY tests: broad randomized sweeps of `(ℓ,n,s)` state
combinations, not just numerical-equivalence re-runs, checking structural
invariants of M02-M04's compiler core (`mgp_phi.jl`, `mgp_transitions.jl`,
`mgp_reduce.jl`) that must hold for EVERY reachable state, not only the
handful of hand-picked instances M02-M04's own acceptance tests already
cover.

## Dependency check (done first, per the task's explicit instruction)
`test/Project.toml` and the root `Project.toml` were both read directly.
Neither lists any QuickCheck-style property-testing package (no
`PropCheck`, `Supposition`, `QuickCheck`, or similar under `[deps]` in
either file). Per the task's instruction, no new dependency was added;
randomized-trial loops were hand-written using Julia's standard `Test` and
`Random` (`MersenneTwister`) machinery, matching the exact pattern already
established in `test/kli_seir_compiled_test.jl`/`test/kli_mers_compiled_test.jl`
(a seeded `MersenneTwister` per testset, `rand(rng, ...)` draws inside a
`for _ in 1:ntrials` loop, `@test` per draw).

## What changed
- Added `test/kli_properties_test.jl`: five named `@testset`s, each a
  broad randomized sweep (400-600 trials per property, ~14489 individual
  `@test`s in total — trial counts vary because some properties assert
  multiple invariants per trial) over every BIRTH/MIGRATION event that
  actually exists in `PhyloPOMP.SEIR`/`PhyloPOMP.MERS` (`infection`,
  `progression`, `transmission_cc`, `transmission_hh`, `transmission_hc`,
  `transmission_ch` — collected generically by filtering `event.type`, not
  hand-listed by name, so a future BIRTH/MIGRATION event added to either
  model is automatically covered without editing this file).
- `test/runtests.jl`: one line added (`include("kli_properties_test.jl")`,
  at the end of the file, after the MERS proposal-kernel test includes).
- **No existing file modified** beyond that one `include` line. `mgp_phi.jl`,
  `mgp_transitions.jl`, `mgp_reduce.jl`, `mgp_seir_filter.jl`,
  `mgp_mers_filter.jl`, `mers_naive.jl`, `seir_naive.jl`, and all 8
  hand-coded filter modules are byte-for-byte unchanged — no bug was found
  that required a fix (see "Findings" below).

## Files added
- `test/kli_properties_test.jl`
- `handoffs/M12_property_testing.md` (this file)

## Files modified
- `test/runtests.jl` — one line (the new `include`).

## The 5 properties

| # | Property | Trials | Result |
|---|---|---|---|
| 1 | Population/genealogy support: `ℓ_d <= n_d`, and no saturation ever implies `s_d > ℓ_d` (checked on inputs AND on `enumerate_saturations`/`full_transitions`/`reduce_event_indicator` output) | 400 trials (1972 `@test`s) | PASS, 0 findings |
| 2 | Full-to-reduced conservation: `Σφ_u == ΣΦ_u`, exact `Rational{Int}` equality, plus a provenance partition check | 500 trials (1500 `@test`s) | PASS, 0 findings |
| 3 | Proposal support: `Φ_u(z) >= 0` (and `φ_u >= 0`) for every BIRTH/MIGRATION event, broad `(ℓ,n)` sweep | 500 trials (2457 `@test`s) | PASS, 0 findings |
| 4 | Impossible transitions have zero compatibility: `n_d < r_d`, `ℓ_d > n_d`, `s_d < 0`, `s_d > r_d` cases all degrade to `φ_u=0`/`Φ_u=0`, never an error or a nonsensical value | 400 trials (1327 `@test`s) | PASS, 0 findings |
| 5 | Genealogy-operation structural invariants: `CrossDemeTransition.ancestral_deme != .target_deme`; `ForkTransition.slot_demes` length `== sum(s)`; `IdentityTransition`/`InlineSameDemeTransition`/`CrossDemeTransition`/`ForkTransition` correspond to `sum(s)` of 0/1/1/`>=2` respectively; reduced-key/kind consistency; no duplicate reduced groups | 600 trials (7233 `@test`s) | PASS, 0 findings |

Total new tests: **14489** (`18759 - 4270`, exactly matching the
before/after delta below).

### Design notes per property
- **P1** checks the invariant on two levels, as the task specified: (a) a
  sanity check that the randomized generator itself only ever produces
  `ℓ_d <= n_d` states (`random_valid_state`'s construction is `n` first,
  then `ℓ ~ Uniform(0,n)` per deme, so this is true by construction, but
  asserted explicitly rather than merely assumed); (b) that
  `enumerate_saturations`/`full_transitions`/`reduce_event_indicator` never
  emit an `s`/reduced-Fork slot-deme count exceeding `ℓ` in the
  corresponding deme, checked directly on their actual output.
- **P2** is the randomized generalization of the exact invariant already
  spot-checked at M04's six hand-picked instances (TCC/THH/THC/TCH/SEIR
  infection/SEIR progression, `test/kli_reduce_test.jl`) and explicitly
  documented in `mgp_reduce.jl`'s `reduce_event_indicator` docstring — this
  property test confirms it holds broadly (500 randomized `(event,ℓ,n)`
  combinations across both models), not just at those six points.
- **P3** checks non-negativity of `Φ_u`/`φ_u` directly (the quantity the
  compiled SEIR/MERS proposals assign probability mass proportional to,
  either directly for THC/TCH-style events or via the `C(ℓ,s)`-weighted
  aggregate M09 Finding 2 documents for TCC/THH-style events) — a basic
  sanity property that should hold trivially given `φ_u` is a product of
  ratios of non-negative binomial coefficients, but confirmed broadly
  rather than assumed, per the task's explicit instruction to verify
  carefully rather than reason from the formula alone.
- **P4** exercises three distinct "structurally infeasible" scenarios per
  the task's examples: `n_d < r_d` (denominator-zero case), `ℓ_d > n_d`
  (the "should never be a valid model state" case, confirmed to degrade to
  `φ_u=0` for every saturation rather than crash — traced to
  `safe_binomial`'s `a<0` branch zeroing `n_d-ℓ_d < 0`), and out-of-range
  `s` (`s_d<0`, `s_d>r_d`) fed directly to `kli_binomial_ratio`, bypassing
  `enumerate_saturations` entirely (the defensive convention M02 already
  found and fixed one real bug in — this property test confirms that fix
  holds broadly, not just at the one instance M02's own test caught it at).
- **P5** is a structural (not numerical) fuzz test: for every
  `KLITransition` produced across 600 randomized `(event,ℓ,n)` draws,
  confirm the type-level invariants (`sum(s)` matches the classification,
  `ancestral_deme == event.from` uniformly, `CrossDemeTransition`'s two
  deme fields always differ, `ForkTransition.slot_demes` length matches
  `sum(s)` exactly) plus reduced-level invariants (`kind == key[1]`,
  `cross`/`fork` key components well-formed, no duplicate reduced groups
  within one event/state's output) — generalizing M03/M04's per-event
  named-instance structural assertions (which checked these facts at
  TCC/THH/THC/TCH/infection/progression specifically) to arbitrary
  randomized states.

## Findings
**None.** All five properties passed on every trial across all 14489
individual `@test`s, with no test-authoring bug found or fixed during this
milestone (unlike M02, which found and fixed a real `kli_binomial_ratio`
bug during its own test-writing). This is a genuine, not merely convenient,
result: the properties were written from the task's own specification and
from M02-M04's already-documented invariants (the docstrings in
`mgp_phi.jl`/`mgp_reduce.jl` already state several of these as claimed
invariants — e.g. `reduce_event_indicator`'s docstring literally states
`sum(rt.Φ for rt in result) == sum(t.phi for t in transitions)` always
holds), so a clean pass here is corroborating evidence that M02-M04's
existing targeted tests were not missing a broader-scope violation, rather
than evidence the properties were too weak to find one. No property's
random generator was narrowed after the fact to avoid a failing case — the
generators used are the ones described above, run once, and all passed.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before this milestone: **4270/4270 passed** (M11's unchanged baseline).
- After this milestone: **18759/18759 passed** (`4270 + 14489`), ~51s wall
  time, 0 failures, 0 errors. The delta is exactly the new
  `"Property-based verification (M12)"` testset in
  `test/kli_properties_test.jl`. No pre-existing testset's pass count
  changed.

## Verification status
| Property | Status | Trials |
|---|---|---|
| P1: population support `ℓ<=n` | PASS | 400 (1972 assertions) |
| P2: full-to-reduced conservation | PASS | 500 (1500 assertions) |
| P3: proposal support `Φ_u>=0` | PASS | 500 (2457 assertions) |
| P4: impossible transitions -> zero | PASS | 400 (1327 assertions) |
| P5: structural invariants (fuzz) | PASS | 600 (7233 assertions) |

## Known failures / unresolved issues
- None found by this milestone. The pre-existing, carried-over unresolved
  items from M09/M10/M11 (whether to further reconcile `mgp_mers.jl`'s
  demography, whether `mgp_filter.jl`'s generic stubs should ever be filled
  in) remain out of this milestone's scope and untouched.
- The randomized sweeps here (400-600 trials per property, `nmax` capped
  at 10-20 per deme for tractability) are broad but not exhaustive — as
  with every prior milestone's randomized testing (M08/M09's Gate-5
  sweeps), a property could in principle be violated outside the sampled
  range. No such violation was found in >14000 total assertions across
  five distinct properties and both models.

## Git state
- branch: `atpabuser-devel`
- commit: none (nothing committed or pushed by any milestone so far, per
  instructions)
- Files changed by this milestone:
  - `test/kli_properties_test.jl` (new)
  - `handoffs/M12_property_testing.md` (new, this file)
  - `test/runtests.jl` (modified — one more `include` line)
- Cumulative uncommitted state (M00-M12): all files listed in M09/M10's
  handoffs, plus `handoffs/M11_model_crossvalidation.md` (M11) and this
  milestone's three files above.

## Resume instructions
1. Advisor reviews `test/kli_properties_test.jl`'s five testsets and this
   handoff's "Findings" section (a clean pass, not a weakened-generator
   result).
2. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` before M13;
   must remain exactly 18759/18759 (this milestone's new baseline).

## Next milestone
M13 — per the master project plan, Indexed Event Families: `@eventfamily`
syntax for arbitrary-dimension models (e.g. a model with `D` demes indexed
generically, rather than SEIR's/MERS's fixed 2-deme structures with
hand-named events like `transmission_cc`/`transmission_hh`/`transmission_hc`/
`transmission_ch`). This is flagged explicitly here, per this dispatch's
own instruction, as a **bigger architectural undertaking than M11/M12**:
it requires extending `mgp_macro.jl`'s `@mgp`/`@event` DSL itself (not just
consuming the existing `Event`/`MGPModel` IR the way M02-M12 have done),
deciding how `Event.from`/`.into`/`.r` generalize to an indexed family of
events parameterized over deme pairs or tuples (rather than one `Event`
struct per concrete named event), and re-verifying that
`production_slots`/`enumerate_saturations`/`kli_binomial_ratio`/
`full_transitions`/`reduce_event_indicator` (all of which are already
demonstrated in M02-M12 to be generic over `event.r`/`event.from`/
`event.type` with no per-event lookup tables) compose correctly with
however the indexed-family DSL expands into concrete `Event`s. Given this
scope, M13 should probably get its own full dispatch rather than being
combined with anything else, consistent with how M08 (SEIR end-to-end) and
M09 (MERS end-to-end) each received their own dedicated dispatch for
comparably-sized pieces of work.

## Context note
The most important thing to preserve if this conversation were compacted:
M12 added `test/kli_properties_test.jl`, five broad randomized property
tests (400-600 trials each, ~14489 total `@test`s) over SEIR/MERS's
existing, already-verified KLI compiler core
(`mgp_phi.jl`/`mgp_transitions.jl`/`mgp_reduce.jl`), generalizing M02-M09's
targeted-instance and isolated-trial patterns into named, broad-coverage
properties: population support (`ℓ<=n`), full-to-reduced Rational-exact
conservation, proposal-support non-negativity, safe degradation on
structurally-infeasible inputs, and structural fuzzing of every
`KLITransition`/`ReducedTransition` subtype's invariants. **All five
properties passed with zero findings** — no bug was found or fixed, unlike
M02 (which found a real `kli_binomial_ratio` negative-`s` bug during its
own test-writing). No dependency was added (no QuickCheck-style library
exists in this repo; hand-written `MersenneTwister`-seeded randomized loops
were used, matching the pattern `kli_seir_compiled_test.jl`/
`kli_mers_compiled_test.jl` already established). M11
(`handoffs/M11_model_crossvalidation.md`, executed in the same dispatch)
is a short labeling milestone concluding, honestly, that **no model in this
repository has "exact oracle" status** — SEIR/MERS both sit at "handwritten
implementation comparison" (Gate 5, exact numerical matches against
independently-hand-coded `seir_naive.jl`/`mers_naive.jl`), and BDEI/BDSS
are both "not yet independently verified" because neither model exists
anywhere in the codebase. Baseline is now **18759/18759** (`4270` M11
baseline `+ 14489` new property tests). M13 (Indexed Event Families,
`@eventfamily` syntax) is flagged as a substantially bigger architectural
undertaking warranting its own dedicated dispatch.
