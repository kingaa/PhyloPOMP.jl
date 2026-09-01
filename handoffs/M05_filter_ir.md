# Handoff — M05: Filter IR (RegularFlow / SingularFlow / DecayTerm)

## Status
COMPLETE

## Objective
Turn the M03/M04 empirical finding — `seir_naive.jl`'s regular-event proposal
never realizes `Fork` (nor `InlineSameDeme`) outside `singular_part!` — into a
principled, KLI-grounded Filter IR classification of M04's
`Vector{ReducedTransition}`: `RegularFlow` (continuous-time regular-event
flow), `SingularFlow` (fixed/observed data-event time only), `DecayTerm` (a
PROVENANCE-ONLY placeholder for KLI-incompatible-at-regular-time mass, not
the computed decay/driver weight — that is later-milestone work). No decay
WEIGHT formula (Eq. 47/B2), no proposal logic (`π_u`), no RC/RH/SC/SH
handling — all explicitly out of scope per the milestone spec.

## What changed
- Added `src/examples/mgp_filter_ir.jl`: abstract `FilterTerm`; concrete
  `RegularFlow`, `SingularFlow`, `DecayTerm` (each carrying `event::Event`,
  `reduced::ReducedTransition` — provenance back through M04's `Φ`/grouping
  to M03's `Vector{KLITransition}` — and a redundant `Φ::Rational{Int}`
  mirror); `DecayTerm` additionally carries `reason::Symbol`
  (`:fork_unobserved_at_regular_time`, the only reason this milestone
  constructs). `FilterSpec(event, regular, singular, decay)` container.
  `classify_filter_terms(event, reduced) -> FilterSpec` (the classification
  pass, with a provenance-integrity guard: every `KLITransition` inside every
  `ReducedTransition` must reference the same `event` argument, or it
  throws). `filter_spec(event, ℓ, n; Q=1) -> FilterSpec` convenience
  composition (`reduced_transitions |> classify_filter_terms`).
- Wired into `src/examples/Examples.jl`, immediately after `mgp_reduce.jl`
  and before `mgp_filter.jl` (which remains untouched).
- Added `test/kli_filter_ir_test.jl` (63 new `@test`s, registered in
  `test/runtests.jl` right after `kli_reduce_test.jl`): SEIR
  `infection`/`progression`, MERS `transmission_cc`/`transmission_hc`
  classification checks (re-derived via the actual function call chain, not
  just trusted from the M04 handoff table), a provenance-chain test walking
  `FilterTerm -> ReducedTransition -> KLITransition -> Event`, a mismatched-
  provenance-throws test, and a synthetic-`Event` test for the
  `!event.regular` branch (see "Gap" section below).

## Files added
- `src/examples/mgp_filter_ir.jl`
- `test/kli_filter_ir_test.jl`
- `handoffs/M05_filter_ir.md` (this file)

## Files modified
- `src/examples/Examples.jl` — one line added: `include("mgp_filter_ir.jl")`.
- `test/runtests.jl` — one line added: `include("kli_filter_ir_test.jl")`.
- `mgp_filter.jl` and all 8 hand-coded filter modules untouched (read-only
  oracles, cited by file:line, never edited).

## Classification table

| Event | `event.regular` | `ReducedTransition.kind` | Bucket | Oracle citation |
|---|---|---|---|---|
| SEIR `infection` (BIRTH, r=(1,1)) | true | `:noop` | `RegularFlow` | `seir_naive.jl:146-149` (k==1, untracked-parent branch — no `swap!`/`fork!` call) |
| SEIR `infection` | true | `:cross` | `RegularFlow` | `seir_naive.jl:150-156` (k==2, tracked-parent branch, unconditional `swap!(cols,Infec,Expos,b)`) |
| SEIR `infection` | true | `:fork` | `DecayTerm` (placeholder) | `seir_naive.jl:88,90` — `fork!` for infection-type branch points appears ONLY inside `singular_part!`'s `n.type==Node` branch, never in `regular_part!` |
| SEIR `progression` (MIGRATION, r=(0,1)) | true | `:noop` | `RegularFlow` | `seir_naive.jl:157-160` (k==3, untracked-E, no `swap!`) |
| SEIR `progression` | true | `:cross` | `RegularFlow` | `seir_naive.jl:161-167` (k==4, tracked-E, `swap!(cols,Expos,Infec,b)`) |
| SEIR `progression` | true | `:fork` | N/A — structurally impossible (M04 finding: single slot, always at target deme, `sum(s)∈{0,1}` only) | — |
| MERS `transmission_cc` (BIRTH, r=(2,0)) | true | `:noop` | `RegularFlow` | `mers_naive.jl:187-190` (`regular_part!` k==1: `Sc-=1; Ic+=1; ll += log(1-(ellc*(ellc-1)/Ic/(Ic-1)))` — no `fork!`/`swap!` call) |
| MERS `transmission_cc` | true | `:fork` | `DecayTerm` (placeholder) | `mers_naive.jl:71-78` (`singular_part!`'s Node branch, `n.lineage ∈ cols[Camel]`, `k==1` "camel-camel" sub-case) — `fork!(cols, Camel, ..., (Camel,Camel), children)` appears ONLY in `singular_part!` |
| MERS `transmission_hc` (THC, BIRTH, r=(1,1)) | true | `:noop` | `RegularFlow` | `mers_naive.jl:195-198` (`regular_part!` k==3: no cols mutation) |
| MERS `transmission_hc` | true | `:cross` | `RegularFlow` | `mers_naive.jl:199-205` (`regular_part!` k==4: `swap!(cols, Camel, Human, b)`) |
| MERS `transmission_hc` | true | `:fork` | `DecayTerm` (placeholder) | `mers_naive.jl`'s Node branch, `k==2`/`k==3` "camel-human" sub-case — `fork!(cols, Camel, ..., (Camel,Human)/(Human,Camel), children)`, singular-only |
| (synthetic) `regular=false` BIRTH event | false | `:noop`/`:cross`/`:fork` (all) | `SingularFlow` | No real SEIR/MERS example exists (see "Gap" below); rule justified by `mgp_filter.jl:166-167`'s `regular_step!` already zeroing `alpha`/`pi` for `!ev.regular` events at every non-fixed time |

`RC`/`RH`/`SC`/`SH` (DEATH/SAMPLE, r=(0,0)) are entirely out of scope: M03's
`full_transitions` already throws `ArgumentError` on these event types, so
`reduced_transitions` produces no input for `classify_filter_terms` to
classify for them — no parallel path was invented, per the milestone's
explicit instruction.

## Mathematical grounding for excluding `:fork` from `RegularFlow`
`mers_filter_suite.tex` explains this is principled, not an implementation
quirk of the naive proposal:
- Lines 760-762 (TCC/THH): "The missing mass corresponds to unobserved
  within-camel branch points and is accounted for through the regular birth
  inflow/outflow balance (§8), not through λ."
- Lines 841-844 ("No separate branch-point hazard"): "The within-deme unseen
  branch-point loss arises automatically from Eqs 7-8 inflow against the
  TCC/THH share of the (full-birth) outflow ... Adding either to λ would
  double-count."
- Lines 846-865 (driver box): births are kept FULL in the reduced regular
  exit rate (`α_TCC+α_THH+α_THC+α_TCH`, no ℓ-dependent reduction), while the
  INFLOW terms `I_7..I_12` (Eqs 7-12, lines 735-789) are built ONLY from the
  non-fork (noop/cross) reduced cases. The imbalance between the full
  outflow and the non-fork-only inflow is where the fork/branch-point mass
  goes — it must not be double-counted as `RegularFlow`.

## `DecayTerm`'s placeholder: what it carries vs. what remains
`DecayTerm` carries ONLY:
- `event` (the source `Event`)
- `reduced::ReducedTransition` (the `:fork` group, with full M03/M04
  provenance intact)
- `Φ::Rational{Int}` — copied VERBATIM from `reduced.Φ` (the raw sum of
  `phi_u` over the collapsed full transitions)
- `reason::Symbol` (`:fork_unobserved_at_regular_time`)

It explicitly does **NOT** carry KLI's actual `λ(t,x,y)` (Eq. 47/B2), which:
(a) is a rate-integrated quantity involving the population hazard `α_u`, not
just `Φ_u`; (b) per the tex (lines 841-844), is a *different mechanism* from
this fork/branch-point mass (`λ`'s own formula is sampling hazard +
sub-threshold removal, lines 833-836 — RC/RH/SC/SH territory, which this
milestone's classifier never even sees). Working out the exact mechanism by
which this fork mass is resolved (inflow/outflow imbalance vs. a decay-like
subtraction vs. some other Eq. 47/Appendix-B-B2 construction) is explicitly
deferred to the later driver/boost/decay milestone (M07 per the master
project plan) — this file only records that the mass exists, where it came
from (full provenance chain intact), and that it must not be silently folded
into `RegularFlow`.

## Gap found (as instructed, not fabricated)
No BIRTH/MIGRATION event with `event.regular == false` currently exists in
either shipped model. Checked exhaustively (asserted in
`test/kli_filter_ir_test.jl`'s `"synthetic singular BIRTH event"` testset):
every `regular == false` event in SEIR (`sampling`) and MERS
(`sampling_c`, `sampling_h`) is `SAMPLE`-type, which M03's `full_transitions`
already rejects by type — so no `ReducedTransition` for a singular
BIRTH/MIGRATION event can ever reach `classify_filter_terms` via
`reduced_transitions` in the current model set. The `!event.regular` branch
of `classify_filter_terms` is validated only against a synthetic, hand-built
`Event` (same structural fields as SEIR's `infection`, `regular=false`) —
documented in `mgp_filter_ir.jl`'s and the test file's comments as a known
gap, not silently assumed away. The classification rule chosen for that
branch (every reduced transition of a singular event → `SingularFlow`,
regardless of kind) is justified by `mgp_filter.jl:166-167`'s
`regular_step!`, which already routes a `!ev.regular` event's contribution
to zero on every open interval — so there is no decay bucket needed for a
singular event's own outcomes (whichever one the fixed observed tree
dictates at that instant is exactly consistent by construction). This
should be re-examined once a real singular BIRTH/MIGRATION event (or the
driver/decay milestone's treatment of SAMPLE-type singular events) exists to
check it against.

## Architecture decisions
- `classify_filter_terms` takes `event::Event` explicitly (not inferred from
  `reduced`) and verifies every `ReducedTransition`'s provenance
  (`transitions[i].event`) matches it — a defensive integrity check against
  a caller accidentally classifying one event's reduced transitions under a
  different `event` argument. Throws `ArgumentError` on mismatch (matches
  M03's "fail loud" convention for structurally-wrong calls).
- `RegularFlow`/`SingularFlow`/`DecayTerm` all carry a redundant `Φ` field
  (mirroring `reduced.Φ`) purely for ergonomic access without unpacking
  `.reduced.Φ` — same pattern M04 used for `ReducedTransition.kind`
  mirroring `key[1]`.
- `FilterSpec` does not attempt to also store `ℓ`/`n`/`Q` — those live only
  implicitly via `reduced.transitions[i].s`/`.phi`; adding them would
  duplicate information already reachable through the provenance chain.
- No new struct equality convention introduced; like `ReducedTransition`,
  `FilterTerm` subtypes have identity-based `==` (carry a `Vector`-bearing
  `ReducedTransition` field), so tests compare fields explicitly, not `==`
  on whole structs — consistent with M04's documented convention.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before (M04 baseline): **3920/3920 passed**.
- After this milestone: **3983/3983 passed** (`3920 + 63`), ~42s wall time,
  0 failures, 0 errors. The delta is exactly the new `Filter IR
  classification (M05)` testset in `test/kli_filter_ir_test.jl`. No
  pre-existing testset's pass count changed.
- One test-authoring bug found and fixed during this milestone: an initial
  assertion `spec.regular[1].Φ + spec.decay[1].Φ == 1` for MERS
  `transmission_cc` was WRONG — `Σ_s φ_u(s)` is not in general a probability
  distribution over saturations `s` (it is a sum of binomial-ratio
  compatibility factors, not a normalized categorical kernel); at the tested
  instance (`ℓ_C=2, I_C=5`), `Φ(noop)+Φ(fork) = 3/5+1/10 = 7/10 ≠ 1`. Fixed
  by removing the incorrect identity and documenting why in a comment
  instead.

## Verification status
| Layer | Status | Oracle |
|---|---|---|
| `classify_filter_terms` regular-event classification (`:noop`/`:cross` → `RegularFlow`, `:fork` → `DecayTerm`) | Verified (new) | `test/kli_filter_ir_test.jl`, all model-specific testsets |
| SEIR `infection`: 2 `RegularFlow` (8/15, 1/10), 1 decay term (1/30) | Verified (new); numbers re-derived via the live function chain, not trusted from M04's handoff table | `"SEIR infection"` |
| SEIR `progression`: 2/2 `RegularFlow`, empty decay bucket | Verified (new) | `"SEIR progression"` |
| MERS `transmission_cc`: 1 `RegularFlow` (3/5), 1 decay term (1/10) | Verified (new) | `"MERS transmission_cc"` |
| MERS `transmission_hc`: 2 `RegularFlow`, 1 decay term | Verified (new) | `"MERS transmission_hc"` |
| Provenance chain (`FilterTerm → ReducedTransition → KLITransition → Event`) unbroken | Verified (new) | `"provenance"` testset |
| Mismatched-provenance guard throws | Verified (new) | `"provenance"` testset, `@test_throws` |
| `!event.regular` branch | Verified only against a synthetic `Event` — no real SEIR/MERS example exists | `"synthetic singular BIRTH event"` testset; gap documented above |
| `DecayTerm`'s actual weight (Eq. 47/B2) | NOT implemented — explicitly out of scope, deferred further (later driver/boost/decay milestone) | N/A |
| RC/RH/SC/SH classification | NOT implemented — no `ReducedTransition`s ever reach this pipeline for them (M03's scope boundary) | N/A |
| Proposal logic (`π_u`) | NOT implemented — explicitly out of scope | N/A |

## Known failures / unresolved issues
- None outstanding in the delivered code (the one test-authoring bug found
  during this milestone's own test-writing was fixed within the milestone,
  see "Tests run").
- Gap carried forward (not a bug, a scope limitation documented above): no
  real singular BIRTH/MIGRATION event exists in SEIR/MERS today, so the
  `!event.regular` branch of `classify_filter_terms` is validated only
  against a synthetic `Event`. Flagging for the advisor: if a future model
  (or a reinterpretation of SEIR/MERS's `sampling`/`sampling_c`/`sampling_h`
  events once M03's SAMPLE-type scope boundary is revisited) introduces a
  genuine singular BIRTH/MIGRATION event, the `SingularFlow`-for-everything
  rule chosen here should be re-checked against that concrete case.
- Judgment call inherited from M04's handoff, now resolved by this
  milestone: M04 asked whether the `(:noop,...)` `ReducedTransition`'s `Φ`
  (which includes the InlineSameDeme contribution) is fully
  `RegularFlow`-eligible. This milestone's answer: YES — `InlineSameDeme`
  collapses into the SAME reduced `(:noop, d)` class as `Identity` (M04's
  collapsing rule, justified by `swap!`'s literal same-deme no-op
  behavior), so on the REDUCED color-only state there is no way to
  distinguish "untracked parent, no lineage effect" from "tracked parent,
  stays representing the same deme" — both are the identical `y'=y` outcome
  on the pruned genealogy. Nothing about `seir_naive.jl`'s specific choice
  to always realize CrossDeme (never InlineSameDeme) for a tracked-parent
  regular event contradicts this: `seir_naive.jl` never needs to
  distinguish "InlineSameDeme was chosen" from "Identity was chosen" at
  runtime BECAUSE they are literally the same reduced outcome — the
  `Φ(noop)` value correctly aggregates both, and `RegularFlow`'s `Φ` uses
  that aggregate, not a hypothetical split.

## Git state
- branch: `atpabuser-devel`
- commit if created: none
- uncommitted files (cumulative, M00-M05; nothing committed or pushed by any
  milestone so far, per instructions):
  - `docs/compiler/architecture_before.md`, `docs/compiler/compiler_roadmap.md` (M00)
  - `handoffs/M00_reconnaissance.md` .. `handoffs/M04_reduce_m.md` (M00-M04)
  - `src/examples/mgp_audit.jl`, `test/population_ir_test.jl` (M01)
  - `src/examples/mgp_phi.jl`, `test/kli_phi_test.jl` (M02)
  - `src/examples/mgp_transitions.jl`, `test/kli_full_transitions_test.jl` (M03)
  - `src/examples/mgp_reduce.jl`, `test/kli_reduce_test.jl` (M04)
  - `handoffs/M05_filter_ir.md` (new, this file)
  - `src/examples/mgp_filter_ir.jl` (new)
  - `test/kli_filter_ir_test.jl` (new)
  - `src/examples/Examples.jl` (modified — one more `include` line added)
  - `test/runtests.jl` (modified — one more `include` line added)

## Resume instructions
1. Advisor reviews `src/examples/mgp_filter_ir.jl` and
   `test/kli_filter_ir_test.jl`, in particular: (a) whether the
   `DecayTerm`-is-provenance-only design (no computed weight) is the right
   shape for M07 to build on, (b) the `mers_filter_suite.tex` line citations
   (760-762, 841-844, 846-865) actually support "fork excluded from regular
   flow" as a general principle (not just an SEIR/MERS coincidence), and (c)
   the synthetic-`Event` gap for the `!event.regular` branch.
2. M06 ("Provenance and `@mgpaudit`") should treat `FilterSpec` (this
   milestone's output) as the top of the provenance chain to walk: an
   `explain(term)` function following `FilterTerm -> ReducedTransition ->
   KLITransition -> Event` (all now in place, verified unbroken by this
   milestone's `"provenance"` testset), and an `@mgpaudit MODEL` macro
   printing the whole derivation for a model.
3. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any M06
   change and confirm the count does not regress below 3983/3983 (this
   milestone's new baseline).

## Next milestone
M06 — Provenance and `@mgpaudit`: a thin macro wrapper over ordinary
functions, likely `explain(term)` walking the full provenance chain built by
M01 (`audit_model`) through M05 (`FilterSpec`), and an `@mgpaudit MODEL`
macro showing the whole derivation for a model in one readable printout.
`FilterSpec`'s `regular`/`singular`/`decay` vectors, and each `FilterTerm`'s
`event`/`reduced`/`Φ` fields, are the concrete data M06 should walk — no
further Filter IR restructuring should be needed; M06 is a presentation
layer over what M01-M05 already built.

## Context note
The most important thing to preserve if this conversation were compacted:
`src/examples/mgp_filter_ir.jl`'s `classify_filter_terms(event, reduced) ->
FilterSpec` is now the canonical, tested classification from M04's
`Vector{ReducedTransition}` to the Filter IR's three buckets. The rule:
for a REGULAR BIRTH/MIGRATION event, `:noop`/`:cross` → `RegularFlow`,
`:fork` → `DecayTerm` (a PROVENANCE-ONLY placeholder, NOT the computed Eq.
47/B2 decay weight — that is later-milestone work, M07 per the master
plan); for a SINGULAR (`event.regular==false`) BIRTH/MIGRATION event, every
reduced transition → `SingularFlow`. This is verified against
`seir_naive.jl`'s and `mers_naive.jl`'s actual `regular_part!`/
`singular_part!` code (file:line citations in the classification table
above), and grounded in `mers_filter_suite.tex`'s explanation (lines
760-762, 841-844, 846-865) of why fork/branch-point mass is excluded from
regular flow (accounted for via the birth inflow/outflow balance, NOT the
same mechanism as the tex's own explicit `λ(t,x,y)` formula). A real
singular BIRTH/MIGRATION event does not currently exist in SEIR/MERS — this
is a documented gap, validated only via a synthetic `Event` in the test
suite, not silently assumed. `RegularFlow`/`SingularFlow`/`DecayTerm` all
carry unbroken provenance back through `ReducedTransition` to M03's
`Vector{KLITransition}` and the source `Event`, verified by a dedicated
provenance testset. Baseline is now 3983/3983.
