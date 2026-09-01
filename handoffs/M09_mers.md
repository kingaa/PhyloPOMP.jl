# Handoff — M09: End-to-End MERS Compiler

## Status
COMPLETE

## Objective
Apply M08's two general SEIR findings (decay-leftover, Q_u regular/singular
gating) to MERS's two-deme (camel/human), 12-event table, build
`mgp_mers_filter.jl`, and run the master project's Gate 5 (numerical
equivalence) against `mers_naive.jl` — MERS specifically exercises multiple
demes, cross-deme birth, mixed-deme forks, and destructive sampling.

## What changed
- Added `src/examples/mgp_mers_filter.jl`: `mers_compiled_event_rates!`,
  `mers_compiled_regular_part!` (structural mirrors of `NaiveMERS.event_rates!`/
  `regular_part!`, `mers_naive.jl:132-244` — identical control flow and
  RNG-call order, every likelihood term sourced from M01-M07's IR), and
  `mers_compiled_filter_pomp` (mirrors `NaiveMERS.filter_pomp`,
  `mers_naive.jl:252-324`, reusing `NaiveMERS.singular_part!` verbatim, same
  scoping decision M08 made for SEIR).
- Added `test/kli_mers_compiled_test.jl` (34 new `@test`s): TCC/THH
  weighted-aggregate algebraic check, THC/TCH identity+cross check, the M07
  decay instance reproduced through `compiled_decay`, a 400-trial isolated
  `mers_compiled_regular_part!` vs. `NaiveMERS.regular_part!` bit-exact
  check, and a 15-combo Gate-5 end-to-end sweep.
- `src/examples/Examples.jl`: one line (`include("mgp_mers_filter.jl")`).
- `test/runtests.jl`: one line (`include("kli_mers_compiled_test.jl")`).
- **No existing file modified**: `mers_naive.jl`, `seir_naive.jl`, all 8
  hand-coded filter modules, `mgp_filter.jl`'s stubs, `mgp_seir_filter.jl`,
  and `mgp_mers.jl` are all byte-for-byte unchanged.

### Naming note (deviation from the task's suggested names)
The task's prose suggested naming the new functions `compiled_decay`/
`compiled_event_rates!`/`compiled_regular_part!`/`compiled_filter_pomp` (M08's
names). Since `Examples.jl` `include`s every file flat into the single
`module PhyloPOMP` namespace (confirmed: `src/PhyloPOMP.jl:53`), and
`mgp_seir_filter.jl` already defines functions of those exact names with
identical (untyped, same-arity) signatures, reusing them verbatim would
**silently overwrite** the SEIR functions (`Method overwriting is not
permitted during Module precompilation` — confirmed by trying it first).
Renamed to `mers_compiled_event_rates!`/`mers_compiled_regular_part!`/
`mers_compiled_filter_pomp`. `compiled_decay` itself is **not** redefined at
all: it is already fully generic over `model.events` (verified — see
Finding 1), so `mers_compiled_event_rates!` calls the SAME, unmodified
`compiled_decay` from `mgp_seir_filter.jl` directly.

## Finding 1 (M08's decay-leftover pattern): confirmed to generalize directly, unchanged
`mgp_decay.jl`'s `total_decay`/`decay_contribution` and `mgp_seir_filter.jl`'s
`compiled_decay` were already fully generic over `model.events` before this
milestone touched anything (`compiled_decay` loops `for event in
model.events; event.type == DEATH || continue; ...`, no SEIR-only
assumption). No MERS-specific decay function was written; `mgp_mers_filter.jl`
calls `compiled_decay` directly. Verified (not assumed) by reproducing M07's
own hand-derived MERS instance end-to-end through the actual call path:
γ_c=1/2, γ_h=2/5, χ_c=1/10, χ_h=3/10, I_C=5>ℓ_C=2, I_H=3=ℓ_H=3 (boundary):
```
total_decay(MERS,...)                    = 13/5   (tex-literal λ)
leftover(removal_c) = γ_c·I_C·1{I_C>ℓ_C} − γ_c·(I_C−ℓ_C) = 5/2 − 3/2 = 1
leftover(removal_h) = 0  (I_H==ℓ_H, "above" gate false, both terms 0)
compiled_decay(MERS,...) = 13/5 + 1 + 0 = 18/5
```
— exactly M07/M08's cited `mers_naive.jl` total. Reproduced programmatically
in `test/kli_mers_compiled_test.jl`'s "M07 decay instance" testset (input
arithmetic in exact `Rational{Int}`; `compiled_decay`'s own internal
`Float64` conversion means the final comparison is `isapprox`, not `==`,
documented inline).

SAMPLE events (`sampling_c`/`sampling_h`) confirmed to need **no** leftover
correction: `mers_naive.jl:159` adds `chi_c*Ic + chi_h*Ih` unconditionally,
with no companion reduced-rate regular proposal anywhere for these marks
(SAMPLE events are singular, `event.regular==false`; `mgp_filter.jl`'s
regular driver already zeroes their `alpha`/`pi`, so there is no "regular
proposal used a reduced rate" gap to compensate — `total_decay`'s SAMPLE
branch already **is** the naive term, term for term).

## Finding 2 (M08's Q_u regular/singular gating): confirmed for THC/TCH, a genuinely new third case found for TCC/THH
**THC** (`transmission_hc`, r=(1,1), from=Camel) and **TCH**
(`transmission_ch`, r=(1,1), from=Human) are structurally identical to
SEIR's `infection` — confirmed directly against `mers_naive.jl`'s code
(k==3/5: identity only, no draw; k==4/6: exactly one lineage draw, realizes
CrossDeme only, never InlineSameDeme). The `boost(Φ,π)`/divide-out patterns
from M08 apply unchanged.

**TCC** (`transmission_cc`, r=(2,0)) and **THH** (`transmission_hh`, r=(0,2))
are genuinely different — the task explicitly flagged this as unverified,
and it turned out to matter. `mers_naive.jl`'s k==1/2 branches perform **no**
lineage draw at all and their inline term is `1 − C(ℓ_d,2)/C(n_d,2)`. This
does **not** equal `reduce_event_indicator`'s collapsed `:noop` group
(verified numerically to disagree: ℓ_C=2,I_C=5 gives `Φ_noop=3/5` from
`reduce_event_indicator`, but naive's raw term is `9/10`). The resolution:
because a TCC/THH event never distinguishes *which* of the ℓ_d tracked
lineages is silently "the same-deme parent" (the resulting coloring is
unchanged either way — `IdentityTransition` and `InlineSameDemeTransition`
are both `y'=y`), the aggregate probability must be weighted by the number
of ways each saturation's occupied slots could be filled by the ℓ_d tracked
lineages: `C(ℓ_d,0)·Φ_identity + C(ℓ_d,1)·Φ_inline = Φ_identity + ℓ_d·Φ_inline`.
Verified exactly equal to `mers_naive.jl`'s raw formula at four distinct
`(ℓ,n)` instances, plus a Chu-Vandermonde sum-to-1 cross-check
(`Φ_id + ℓ·Φ_inline + C(ℓ,2)·Φ_fork == 1` exactly), confirming the weighting
generally rather than by coincidence at one instance (`test/kli_mers_compiled_test.jl`
"TCC/THH weighted aggregate"). `ForkTransition` is never realized regularly
by any of the four BIRTH marks (only in `NaiveMERS.singular_part!`).

## Finding 3 (new, not anticipated by the task): `mgp_mers.jl`/`mers_naive.jl` hazard mismatch for `death_c`/`death_h`
`mgp_mers.jl:22-23` declares `death_c`/`death_h` with rate
`(S_c>0 ? B_c : 0.0)`/`(S_h>0 ? B_h : 0.0)` (a presence-gated constant), but
`mers_naive.jl:149-150` actually computes `alpha[11]=Bc*Sc/Nc`,
`alpha[12]=Bh*Sh/Nh` (population-proportional) — genuinely different
formulas, agreeing only in the degenerate case `Sc==Nc`. `birth_c`/`birth_h`
(`mgp_mers.jl:20-21`) **do** match `mers_naive.jl:147-148` exactly — only
the two DEATH-side NEUTRAL marks disagree. `mgp_mers.jl` is shared
model-layer infrastructure from an earlier milestone, not modified here
(out of scope, and modifying it risks other consumers). `mers_compiled_event_rates!`
uses `mers_naive.jl`'s own formula directly for `alpha[11]`/`alpha[12]`
(required for bit-exact Gate-5 matching — using `MERS`'s hazard would
desync the regular-step RNG stream the moment `Sc≠Nc`). Flagged for advisor
attention / a later milestone to reconcile `mgp_mers.jl` itself.

## NEUTRAL demography events (`birth_c`/`birth_h`/`death_c`/`death_h`)
Confirmed simple, as expected: `mers_naive.jl:223-230`'s k∈{9,...,12}
branches are plain population increments/decrements (`Sc+=1`, `Sh+=1`,
`Sc-=1`, `Sh-=1`) with **no** `ll` term at all — no Φ_u, no decay
complexity, matching the task's expectation exactly. `mers_compiled_regular_part!`
reproduces this directly (see Finding 3 above for the one hazard-value
wrinkle, which is about the *rate* used for timing, not about any
coloring/likelihood machinery for these marks).

## Event table (all 8 non-NEUTRAL marks + the 4 NEUTRAL marks)

| Mark | r | from | Regular saturations realized (naive) | Filter term |
|---|---|---|---|---|
| TCC (`transmission_cc`) | (2,0) | Camel | Identity + InlineSameDeme (**weighted**: `Φ_id+ℓ_C·Φ_inline`); Fork never | New Finding 2 formula |
| THH (`transmission_hh`) | (0,2) | Human | Identity + InlineSameDeme (**weighted**: `Φ_id+ℓ_H·Φ_inline`); Fork never | New Finding 2 formula |
| THC (`transmission_hc`) | (1,1) | Camel | Identity (untracked parent) OR CrossDeme (tracked parent, drawn); InlineSameDeme/Fork never | `boost`, same as SEIR infection |
| TCH (`transmission_ch`) | (1,1) | Human | Identity (untracked parent) OR CrossDeme (tracked parent, drawn); InlineSameDeme/Fork never | `boost`, same as SEIR infection |
| RC (`removal_c`) | (0,0) | Camel | DEATH — decay/rate-governed, not saturation | direct, decay-leftover (Finding 1) |
| RH (`removal_h`) | (0,0) | Human | DEATH — decay/rate-governed, not saturation | direct, decay-leftover (Finding 1) |
| SC (`sampling_c`) | (0,0) | Camel | SAMPLE — singular only, never regular | `total_decay`'s unconditional χ_C·I_C term |
| SH (`sampling_h`) | (0,0) | Human | SAMPLE — singular only, never regular | `total_decay`'s unconditional χ_H·I_H term |
| birth_c/birth_h/death_c/death_h | (0,0) | — | NEUTRAL — no coloring op at all | plain population update, no `ll` term |

Operator classification, φ_u, and reduction for TCC/THH/THC/TCH were already
fully derived and hand-verified at M02-M04 (`test/kli_phi_test.jl`,
`test/kli_full_transitions_test.jl`, `test/kli_reduce_test.jl`) — this
milestone's own contribution is the regular/singular **placement** (which
saturations the naive *proposal* actually visits at regular time, Finding 2
above), re-derived and re-confirmed here specifically for MERS, not
re-cited blindly.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before (M08 baseline): **4236/4236 passed**.
- After this milestone: **4270/4270 passed** (`4236 + 34`), ~48s wall time,
  0 failures, 0 errors. Delta is exactly the new `"Compiled MERS filter
  (M09)"` testset.

## The 5 verification gates

| Gate | Result | Evidence |
|---|---|---|
| 1. Structural (Δ, α, r, event semantics, wiring) | PASS (pre-existing) | `test/mgpaudit_test.jl` ("mgpaudit(MERS) completeness", 74 tests) |
| 2. Full KLI per-event transitions match trusted derivation | PASS (M03, confirmed) | `test/kli_full_transitions_test.jl` "TCC"/"THH"/"THC"/"TCH" (re-exercised live by `mers_compiled_regular_part!`'s `full_transitions` calls, cross-checked bit-exact) |
| 3. Reduced KLI (Φ_u) matches | PASS (M04, confirmed) | `test/kli_reduce_test.jl` "TCC"/"THH"/"THC"/"TCH"; this milestone additionally found `reduce_event_indicator`'s collapsed Φ is the WRONG quantity for TCC/THH's regular weighting (Finding 2) — a genuinely new result, not just a re-cite |
| 4. Filter structure + decay resolution | PASS | Finding 1 above, `test/kli_mers_compiled_test.jl` "M07 decay instance" |
| 5. Numerical equivalence | PASS | See below |

### Gate 5 detail
Same methodology as M08 (task option (b)): seed-matched `Np=1` `pfilter`
runs, `mers_compiled_filter_pomp` vs. a test-local `naive_oracle_filter_pomp`
wrapper (structurally identical to `NaiveMERS.filter_pomp` but parameterized
on `gen` instead of hardcoding the fixed empirical `mers_tree` — required
because `NaiveMERS.filter_pomp` itself takes no genealogy argument at all;
the wrapper calls `NaiveMERS.singular_part!`/`regular_part!` verbatim,
unmodified). Genealogies simulated with `demeset=NaiveMERS.Demes,
samplemap=[NaiveMERS.Camel,NaiveMERS.Human]` so the resulting deme values
are the exact enum instances `mers_naive.jl`'s code compares against.

- **Isolated unit check** (permanent test): 400 randomized `(Sc,Ic,Sh,Ih,
  ℓ_C,ℓ_H,params,dt,seed)` trials, mixed event realizations across all 12
  marks — **0 mismatches**.
- **Permanent Gate-5 sweep** (`test/kli_mers_compiled_test.jl`): 15
  randomized parameter draws (β_cc,β_hh∈[1,4]; **β_hc,β_ch∈[0.3,2.3], never
  zero**, so THC/TCH cross-deme marks are exercised in every combo;
  γ_c,γ_h∈[0.5,2.0]; χ_c,χ_h∈[0.2,0.6]) × a genealogy simulated per combo
  (target sample count 1-4, pop_c/pop_h∈{20,30,40}) × 100 seeds each. All
  15/15 combos successfully built their target genealogy, yielding **692
  finite log-likelihood comparisons, 692 exact matches**, worst
  `|Δll|≈7.1e-15` (floating-point noise, tolerance `atol=1e-6,rtol=1e-8`,
  same justification as M08).
- **Development-time broader sweep** (not committed, run and recorded for
  the advisor): 25 combos × 300 seeds, target sample count 1-5,
  pop_c/pop_h∈{20,30,40,60} — 23/25 combos built a genealogy within the
  retry budget (2 timed out, not a filter failure), **2903/2903 finite
  log-likelihood pairs matched exactly**, worst `|Δll|≈1.42e-14`. Confirms
  no silent gap across a wider parameter/sample-size range, including
  cross-species (camel↔human) samples in most combos.

## Known failures / unresolved issues
- Finding 3 (`mgp_mers.jl`'s `death_c`/`death_h` hazard formula does not
  match `mers_naive.jl`'s actual code) is a genuine, previously-unflagged
  gap — not fixed here (out of scope, shared model-layer file), worked
  around in `mers_compiled_event_rates!` by using `mers_naive.jl`'s formula
  directly. Flagged for advisor attention / a future milestone.
- Same SEIR `χ`/gap pattern noted in M07/M08 does not recur for MERS's own
  `chi_c`/`chi_h` (MERS's `sampling_c`/`sampling_h` events genuinely exist
  in `mgp_mers.jl`, unlike SEIR's phantom second sample type) — no analogous
  gap found here.
- `mgp_filter.jl`'s generic stubs remain unfilled, same as M08.
- Gate 5's sweep (15 combos × 100 seeds permanent, 25 × 300 during
  development) is broad but not exhaustive; no failures were found in
  >3500 total finite-likelihood comparisons across both runs.

## Git state
- branch: `atpabuser-devel`
- commit: none (nothing committed by any milestone so far, per instructions)
- new/modified files beyond M00-M08's:
  - `src/examples/mgp_mers_filter.jl` (new)
  - `test/kli_mers_compiled_test.jl` (new)
  - `handoffs/M09_mers.md` (new, this file)
  - `src/examples/Examples.jl` (modified — one more `include` line)
  - `test/runtests.jl` (modified — one more `include` line)

## Resume instructions
1. Advisor reviews `src/examples/mgp_mers_filter.jl`'s header (three
   Findings) and Finding 3 in particular (a genuine `mgp_mers.jl` bug this
   milestone worked around rather than fixed).
2. Decide whether to fix `mgp_mers.jl:22-23`'s `death_c`/`death_h` rates to
   match `mers_naive.jl`'s `Bc*Sc/Nc`/`Bh*Sh/Nh` (a small, isolated change,
   but touches shared model-layer infrastructure other milestones may
   depend on — verify nothing else consumes `MERS.events`'s `death_c`/
   `death_h` hazard before changing it).
3. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any
   further change; must not regress below 4270/4270.

## Next milestone
M10 — per the master plan, the stress test is "a progression/migration
recoloring [that] is not automatically collapsible with identity". **BDEI
does not exist anywhere in this repository** — confirmed by an explicit
grep (`bdei`, case-insensitive) across every `.jl`/`.tex`/`.md` file: zero
matches. `src/examples/` contains only SEIR and MERS (`mgp.jl`'s `@mgp SEIR`
block, `mgp_mers.jl`'s `@mgp MERS` block); no BDEI model, naive filter, or
test file exists at all, at any milestone.

Flagging for the advisor, not deciding unilaterally: SEIR's `progression`
event (`mgp.jl:196`, `move=swap(E => I)`, a pure MIGRATION with r=(0,1)) is
*already* exactly the "recoloring not collapsible with identity" stress test
the master plan describes for BDEI — `test/kli_reduce_test.jl`'s "SEIR
progression" testset already documents this (`reduce_event_indicator`
returns a single, non-collapsed `IdentityTransition` group with no
`InlineSameDemeTransition`/`Fork` possible at all, since r_E=0 rules out any
same-deme slot). It is plausible M10's actual job is "confirm SEIR's
`progression` already covers this, formally write up why" rather than
defining a new `@mgp BDEI` model from scratch — but a genuine BDEI model
(Birth-Death-Exposed-Infectious, with an explicit E→I progression AND its
own separate birth/death/sampling structure distinct from SEIR/MERS) may
still be the master plan's actual intent, since "BDEI" as a named
phylodynamic model conventionally has different demographic assumptions
(no S compartment, or exponential-growth demography) than SEIR — this
distinction matters for what M10 should build and was not resolved here, by
design (out of this milestone's scope; the task instructed "don't decide
this yourself, just flag it clearly for the advisor").

## Context note
The most important thing to preserve if this conversation were compacted:
M08's decay-leftover pattern generalizes to MERS UNCHANGED (`compiled_decay`
is already fully generic, reused verbatim, confirmed via M07's own 18/5
instance). M08's Q_u regular/singular-gating pattern (use unreduced
`IdentityTransition`/`CrossDemeTransition` directly, not `reduce_event_indicator`'s
collapsed groups) holds for MERS's cross-deme marks THC/TCH exactly as it
did for SEIR's `infection` — but TCC/THH (MERS's WITHIN-deme fork marks,
r=(2,0)/(0,2)) are a genuinely NEW third case M08 could not have predicted:
their regular-time "no visible fork" term is a `C(ℓ,s)`-WEIGHTED aggregate
(`Φ_identity + ℓ_d·Φ_inline`), not `reduce_event_indicator`'s plain
(unweighted) `Φ_noop` sum — verified to disagree numerically (3/5 vs. 9/10 at
the M07 instance) and re-derived from first principles (Vandermonde/
combinatorial weighting, since TCC/THH's regular proposal never performs an
explicit lineage draw, unlike THC/TCH's cross branch). A separate, previously
unflagged bug was also found: `mgp_mers.jl`'s `death_c`/`death_h` NEUTRAL
event hazards do not match `mers_naive.jl`'s actual code (`(S_c>0?B_c:0.0)`
vs. `Bc*Sc/Nc`) — worked around (not fixed) in `mers_compiled_event_rates!`
by using `mers_naive.jl`'s own formula directly. Gate 5 passes: 692/692
finite comparisons matched exactly in the permanent test (worst
`|Δll|≈7.1e-15`), 2903/2903 in a broader development sweep (worst
`|Δll|≈1.42e-14`). Baseline is now 4270/4270 (`4236` M08 baseline `+ 34` new
tests). BDEI does not exist in this repository at all — M10 needs advisor
input on whether it should formalize SEIR's `progression` as already
covering the master plan's stress test, or build a genuine new BDEI `@mgp`
model.
