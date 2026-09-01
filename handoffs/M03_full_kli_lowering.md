# Handoff — M03: Full KLI Compatibility Lowering

## Status
COMPLETE

## Objective
Turn a saturation `s ∈ S_u(ℓ)` (M02's `enumerate_saturations`) plus its
`φ_u(s)` (M02's `kli_binomial_ratio`) into the classified, KLI-meaningful
full coloring transition `Y → Y'` it licenses, per `mers_filter_suite.tex`'s
Step D (lines 459-467): `IdentityTransition` (no slot occupied),
`InlineSameDemeTransition` (one slot occupied, ancestral deme == post-event
deme — a distinct full-KLI transition from Identity even though `y'=y`),
`CrossDemeTransition` (one slot occupied, ancestral deme != post-event deme
— `swap!`), `ForkTransition` (two-or-more slots occupied, common ancestor —
`fork!`). `ChopTransition` is defined as a documented placeholder type only
(DEATH/SAMPLE/NEUTRAL compatibility is decay/rate-driven, out of scope). No
`Q_u` gating (regular-vs-singular reachability), no reduction/marginalization
over `m` (Φ_u, M04), and no wiring into `mgp_filter.jl`'s stubs was
attempted — all explicitly out of scope per the milestone spec.

## What changed
- Added `src/examples/mgp_transitions.jl`: the `KLITransition` abstract type
  and its four derivable concrete subtypes (`IdentityTransition`,
  `InlineSameDemeTransition`, `CrossDemeTransition`, `ForkTransition`) plus
  the stub `ChopTransition`, and `full_transitions(event, ℓ, n; Q=1)`, which
  calls M02's `enumerate_saturations`/`kli_binomial_ratio` and classifies
  each returned saturation generically from `sum(s)` and a comparison of
  the occupied slot's deme(s) against `event.from` (the ancestral deme for
  every slot, per the tex's Step A, lines 438-441) — no per-model/per-event
  literal saturation→classification table.
- Wired into the include chain in `src/examples/Examples.jl`, immediately
  after `mgp_phi.jl` and before `mgp_filter.jl` (which remains untouched).
- Added `test/kli_full_transitions_test.jl` (78 new `@test`s, registered in
  `test/runtests.jl` right after `kli_phi_test.jl`) with term-by-term checks
  against `mers_filter_suite.tex`'s TCC/THH/THC/TCH derivations and Complete
  Reference Table, hand-derived + `seir_naive.jl`-cross-checked SEIR
  `infection`/`progression` tables, boundary cases, and an explicit
  scope-boundary test (`full_transitions` throws `ArgumentError` on every
  DEATH/SAMPLE/NEUTRAL event in both models).

## Files added
- `src/examples/mgp_transitions.jl`
- `test/kli_full_transitions_test.jl`
- `handoffs/M03_full_kli_lowering.md` (this file)

## Files modified
- `src/examples/Examples.jl` — one line added: `include("mgp_transitions.jl")`.
- `test/runtests.jl` — one line added: `include("kli_full_transitions_test.jl")`.
- No existing struct, macro, filter, or likelihood-bearing file touched.
  `Event`, `MGPModel`, `mgp_macro.jl`, `mgp_mers.jl`, `mgp_phi.jl`,
  `mgp_audit.jl`, `mgp_filter.jl`, and all eight hand-coded filter modules
  are byte-for-byte unchanged (`seir_naive.jl` and `mgp.jl`/`mgp_mers.jl`
  were read extensively as oracles, never edited).

## Mathematical decisions

### 1. Migration semantics (resolved, not assumed)
The milestone's hypothesis — that `progression`'s `r=(0,1)` (single slot at
the TARGET deme I, no independent E-slot) means `s_I=1 → CrossDeme`
(`swap!(y,E,I,b)`) and `s_I=0 → Identity` — is **confirmed**, verified
directly against the actual, tested, hand-coded `seir_naive.jl`, not
assumed:
- `src/examples/seir_naive.jl:161-167` (`regular_part!`, `k==4`, the
  tracked-E-lineage branch): `ll += log(ellE)` (the `-log(q)`
  lineage-selection correction, `q=1/ellE`), then
  `ellE, ellI = swap!(cols,Expos,Infec,b)` at line 164 — moves lineage `b`
  from deme `Expos` (E) to deme `Infec` (I) — exactly `swap!(y, event.from,
  target_deme, b)` with `event.from=E`, `target_deme=I`. Then
  `ll -= log(I)` (line 167, `I` post-increment) is exactly `log(φ_u(s=1))`
  once the `-log(pi[4])=-log(ellE/E)` term already subtracted at line 145
  and the `+log(ellE)` correction at line 162 are algebraically combined
  (they cancel `ellE`'s contribution entirely, per KLI's Boost `Ψ_u =
  φ_u/π_u` decomposition).
- `src/examples/seir_naive.jl:157-160` (`k==3`, untracked-E branch): no
  `swap!` call, coloring unchanged — the `Identity` case — and
  `ll += log(1-ellI/I)` (line 160, `I` post-increment, `ellI` unchanged
  since untracked) is exactly `log((n_I-ℓ_I)/n_I)`, matching `φ_u(s=0)`.
- **Correction to the milestone's own hypothesized closed form**: the
  milestone spec suggested the CrossDeme case might have `φ_u = ℓ/I`. This
  is **wrong**. Both the general M02 formula (`r=(0,1)` gives
  `φ(s=1)=C(n-ℓ,0)/C(n,1)=1/n`, with **no** `ℓ`-dependence) and the
  `seir_naive.jl` line-156/167 algebra above agree exactly:
  **`φ_u(s=1) = 1/I` (post-event `I`), not `ℓ/I`.** Asserted in
  `test/kli_full_transitions_test.jl`, testset `"SEIR progression"`.
- `event.from`/`event.into` for `progression`: `mgp.jl:96-97`
  (`Event(:progression, ..., [0, 1], MIGRATION, 1, [2], ...)`) — `from=1`
  (E), production `r=[0,1]` (the slot lives at I, `into=[2]`), confirming
  the "production at target, not source" encoding is correct as designed
  by `_decode_move`'s `:swap` branch (`mgp_macro.jl:87-93`).

### 2. TCC/THH/THC/TCH term-by-term match against `mers_filter_suite.tex`
Every saturation's `φ_u` value AND operator classification, generated by
the generic `full_transitions` (no per-event literal table), was checked
against the tex exactly, at these line ranges:

| Event | tex derivation | tex table rows | Checked |
|---|---|---|---|
| TCC (r=(2,0)) | 472-480 | 554-556 | 3/3 saturations, `test/kli_full_transitions_test.jl` `"TCC"` |
| THH (r=(0,2)) | 482-485 | 557-559 | 3/3 saturations, `"THH"` |
| THC (r=(1,1)) | 487-498 | 560-563 | 4/4 saturations, `"THC"` — including the `(0,1)`-CrossDeme vs `(1,0)`-InlineSameDeme distinction the spec calls out as most error-prone |
| TCH (r=(1,1)) | 500-509 | 564-567 | 4/4 saturations, mirrored roles, `"TCH"` |
| Step A (ancestral deme = parent's deme) | 438-441 | — | confirmed generic (`event.from` used uniformly) reproduces every one of the above |
| Step D (operator classification rule) | 459-467 | — | `full_transitions`'s `sum(s)`/ancestral-match classification reproduces it exactly |

All four MERS events matched exactly, term by term — no discrepancy found.
`ForkTransition`'s `ancestral_deme`/`slot_demes` fields were also checked
against the tex's `κ^{bb'}_{d_anc,d_b,d_b'}` subscripts (e.g. THC's `(1,1)`
gives `ancestral_deme=1(C)`, `slot_demes={1,2}={C,H}`, matching
`κ^{bb'}_{CHC}`).

### 3. Why RC/RH/SC/SH/NEUTRAL are out of scope
Per `mers_filter_suite.tex:511-539`: removal (RC/RH) has `r=(0,0)` forced,
`φ=1` (the empty product, Eq. 9), and its actual compatibility is a
**decay/outflow** condition (`I_d-1 ≥ ℓ_d` for the outflow; the sub-threshold
case `I_d ≤ ℓ_d` instead contributes to `λ` — tex lines 511-517, deferred to
`§8`/the driver-decay milestone). Sampling (SC/SH) similarly has `r=(0,0)`,
`φ=1`, and splits by whether the event time is an observed data-event time
(singular chop `χ^{a∅}_C`) or not (contributes to `λ` only) — tex lines
519-539. Neither case is a saturation/`φ_u` derivation problem — the
`enumerate_saturations`/`kli_binomial_ratio` machinery is degenerate for
them — so `full_transitions` rejects all `DEATH`/`SAMPLE`/`NEUTRAL` events
by **type**, and `ChopTransition` is defined only as a placeholder so the
IR vocabulary is complete.
**One nuance found and documented (not in the original milestone framing)**:
`r_u=(0,...,0)` is NOT universal for `SAMPLE` events — SEIR's singular
`sampling` event (built from `move=sample(I)`, `mgp.jl:105-106`, decoded via
`mgp_macro.jl`'s `:sample` branch) has `r=(0,1)`, a genuine production slot
for the inserted leaf, unlike MERS's `sample_remove`-based
`sampling_c`/`sampling_h` (`r=(0,0)`). Both are still correctly rejected by
`full_transitions` because the scope boundary is `event.type`, not whether
`r` happens to be all-zero — verified explicitly in
`test/kli_full_transitions_test.jl`'s `"out-of-scope event types"` testset,
which asserts both the differing `r` values and the uniform rejection.

## Architecture decisions
- `full_transitions(event::Event, ℓ, n; Q=1) -> Vector{KLITransition}`, one
  arity only (unlike M02's raw-vector/Event two-arity design) — the
  classification step genuinely needs `event.from`/`event.type` (not just
  `r`), so there is no meaningful raw-vector overload; callers needing a
  bare-vector interface can still call M02's functions directly.
- Every non-`ChopTransition` subtype carries `event`, `s`, `phi` as common
  fields (by convention, not a shared field via a mutual supertype
  requirement) plus type-specific structural fields
  (`deme`/`ancestral_deme`+`target_deme`/`ancestral_deme`+`slot_demes`) —
  deme indices only, no lineage IDs (explicitly out of scope per the
  milestone spec; lineage-ID bookkeeping is a runtime/execution concern).
- Classification order note: `sum(s) >= 2` is treated uniformly as `Fork`
  regardless of exact count (only 0, 1, or 2 occur in SEIR/MERS today,
  since max `r_u` total is 2) — this is forward-compatible with a
  hypothetical `r_u` summing to 3+ in a future model, though the resulting
  "branch point" semantics for 3+ occupied slots is not validated by any
  current test (no such event exists in SEIR/MERS).
- `full_transitions` throws `ArgumentError` (not silently returning `[]`)
  for out-of-scope event types and for `event.from == 0` — matches M02's
  "fail loud on invalid input, don't silently degrade" pattern for
  structurally-wrong calls (as opposed to M02's own `s`-out-of-range
  defensive-zero convention, which is about numerically-degenerate-but-
  well-typed input).

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before (M02 baseline): **3777/3777 passed**.
- First run of this milestone's new suite: **3 failures + 1 error** out of
  the new tests, all self-inflicted test-authoring bugs, not bugs in
  `full_transitions` itself:
  1. THH test reused M02's `ℓ_H=1` concrete instance, which was valid for
     directly exercising the raw `kli_binomial_ratio` formula (M02's test
     doesn't call `enumerate_saturations`) but is combinatorially too small
     for `full_transitions` (which enumerates via `min(r,ℓ)`, so
     `min(r_H=2,ℓ_H=1)=1` only ever yields `s_H∈{0,1}`, never the `s_H=2`
     fork case) — fixed by switching to `ℓ_H=2, n_H=6` and re-deriving the
     three φ values by hand.
  2. The "out-of-scope event types" test wrongly assumed `r=(0,0)` for
     every `SAMPLE`/`DEATH`/`NEUTRAL` event; SEIR's `sampling` actually has
     `r=(0,1)` (see "Mathematical decisions" §3) — fixed by asserting the
     correct expected `r` per event instead of a blanket `zeros(...)`.
- After both fixes: **3855/3855 passed** (`3777 + 78`), ~41s wall time, 0
  failures, 0 errors. The delta is exactly the new `Full KLI compatibility
  lowering (M03)` testset in `test/kli_full_transitions_test.jl`. No
  pre-existing testset's pass count changed.

## Verification status
| Layer | Status | Oracle |
|---|---|---|
| `full_transitions` generic classification (no per-event lookup table) | Verified (new) | `test/kli_full_transitions_test.jl`, all testsets |
| TCC (3/3 saturations, φ + operator) | Verified (new) | `mers_filter_suite.tex:472-480,554-556` |
| THH (3/3 saturations, φ + operator) | Verified (new) | `mers_filter_suite.tex:482-485,557-559` |
| THC (4/4 saturations, φ + operator, InlineSameDeme-vs-CrossDeme direction) | Verified (new) | `mers_filter_suite.tex:487-498,560-563` |
| TCH (4/4 saturations, φ + operator, mirrored roles) | Verified (new) | `mers_filter_suite.tex:500-509,564-567` |
| SEIR `infection` (4/4 saturations, φ hand-derived) | Verified (new) | general M02 binomial formula + structural analogy to THC; CrossDeme case independently confirmed against `seir_naive.jl:150-156` |
| SEIR `progression` (2/2 saturations, migration-semantics resolution) | Verified (new); one hypothesis-correction | `seir_naive.jl:157-167`, exact algebraic match |
| Boundary: `ℓ=0` forces Identity-only | Verified (new) | `test/kli_full_transitions_test.jl` `"boundary ell=0"` |
| Boundary: `ℓ=n` forces `φ=0` on untracked-outcome saturations | Verified (new) | `test/kli_full_transitions_test.jl` `"boundary ell=n"` (progression + THC) |
| Scope: DEATH/SAMPLE/NEUTRAL rejected (type-based, not r-based) | Verified (new); one nuance found (SEIR `sampling` r≠0) | `test/kli_full_transitions_test.jl` `"out-of-scope event types"` |
| `ChopTransition` derivation | NOT implemented — documented stub, explicitly out of scope | N/A, deferred to driver/decay milestone |
| `Q_u` regular-vs-singular gating (e.g. Fork only reachable at singular events) | NOT implemented — explicitly out of scope | N/A, deferred (see "Known failures" below) |
| Reduction over `m` (Φ_u) | NOT implemented — explicitly out of scope | N/A, deferred to M04 |

## Known failures / unresolved issues
- None outstanding in the delivered code — both test-authoring bugs found
  during this milestone's own test-writing were fixed within the milestone
  (see "Tests run"), not deferred.
- **Judgment call flagged for the advisor — `seir_naive.jl`'s regular-event
  proposal only ever visits a SUBSET of the full compatible outcome set,
  and this is NOT a disagreement with `full_transitions`'s φ_u math.**
  During the migration-semantics/infection cross-check, I found that
  `seir_naive.jl`'s `regular_part!` NEVER exercises the
  `InlineSameDemeTransition` (`s=(0,1)`, tracked I-parent lineage stays
  represented in I) or `ForkTransition` (`s=(1,1)`) cases for a REGULAR
  `infection` event — its `k==2` branch (tracked parent) unconditionally
  performs the CrossDeme swap (`s=(1,0)`, `seir_naive.jl:153`), never the
  InlineSameDeme alternative, and Fork is only ever constructed in
  `singular_part!` (observed branch nodes, `seir_naive.jl:88,90`), never in
  `regular_part!`. I worked out why this is consistent, not a bug: between
  fixed singular (observed) event times, the pruned/tracked lineage set
  `cols` cannot gain new members (new lineage IDs only enter via `plant!`
  at the root or as children of `chop!`/`fork!` calls inside
  `singular_part!`, i.e. only at the FIXED times the input tree dictates) —
  so a REGULAR event's `Fork` saturation (`s=(1,1)`, requiring a NEW branch
  point) has `Q_u=0` always: it would imply an unobserved coalescence of
  the pruned tree at a time the fixed input tree says nothing branches.
  This is exactly the `Q_u`/regular-vs-singular GATING that M02's handoff
  flagged as M03's deferred work and that this milestone's spec explicitly
  keeps out of scope (`full_transitions` returns the full, ungated
  structure; gating is later-milestone/driver territory) — so I did not
  attempt to encode this gating in `full_transitions` itself, only
  documented the reasoning here and in `mgp_transitions.jl`'s module
  docstring. Whether `seir_naive.jl`'s specific choice to ALWAYS realize
  CrossDeme (rather than sometimes InlineSameDeme) for a tracked-parent
  regular infection is itself a validated/unbiased importance-sampling
  scheme is a filter-correctness question outside this milestone's scope
  (mgp_filter.jl's stubs / the "Validation gates" section) — flagging for
  the advisor in case it matters for M04+ design, but it did not block or
  contradict any φ_u value asserted here (the CrossDeme φ_u value itself
  was independently confirmed exactly, algebra shown above).
- No ambiguity found in `mers_filter_suite.tex` itself that blocked this
  milestone — Steps A-D and the TCC/THH/THC/TCH derivations were
  unambiguous and reproduced exactly by the generic classifier.

## Git state
- branch: `atpabuser-devel`
- commit if created: none
- uncommitted files (cumulative, M00-M03; nothing committed or pushed by
  any milestone so far, per instructions):
  - `docs/compiler/architecture_before.md`, `docs/compiler/compiler_roadmap.md` (M00)
  - `handoffs/M00_reconnaissance.md`, `handoffs/M01_population_ir.md`, `handoffs/M02_kli_phi.md` (M00-M02)
  - `src/examples/mgp_audit.jl`, `test/population_ir_test.jl` (M01)
  - `src/examples/mgp_phi.jl`, `test/kli_phi_test.jl` (M02)
  - `handoffs/M03_full_kli_lowering.md` (new, this file)
  - `src/examples/mgp_transitions.jl` (new)
  - `test/kli_full_transitions_test.jl` (new)
  - `src/examples/Examples.jl` (modified — one `include` line added)
  - `test/runtests.jl` (modified — one `include` line added)

## Resume instructions
1. Advisor reviews `src/examples/mgp_transitions.jl` and
   `test/kli_full_transitions_test.jl`, in particular: (a) the migration-
   semantics resolution and its correction of the milestone's own `ℓ/I`
   hypothesis to `1/I` (see "Mathematical decisions" §1), (b) the
   `seir_naive.jl` regular-vs-singular `Q_u`-gating observation flagged
   above (not implemented, only documented), and (c) the SEIR `sampling`
   `r≠(0,0)` nuance (§3).
2. M04 ("Explicit Reduction over `m`") should treat `full_transitions`'s
   output as its input: group the `KLITransition`s that collapse to the
   same reduced `d → d'` transition on the color-only state (in particular,
   `IdentityTransition` and `InlineSameDemeTransition` both have `d'=d`) and
   compute `Φ_u(d,d') = Σ φ_u` over the full transitions reducing to that
   `d→d'`, per KLI's `m`-marginalization (mentioned in `mers_filter_suite.tex`
   around lines 430-432, "When several collapse under the reduction, `Φ_u`
   is their sum ... by the Chu–Vandermonde identity") and referenced in
   `mgp_filter.jl`'s `apply_move!` docstring (the "Chu–Vandermonde collapse
   (Eq. 8)" validation gate). `ForkTransition`s reducing together (e.g. a
   symmetric pair) may also need this treatment — work out from the tex
   whether Fork ever collapses with anything else, or is always its own
   reduced class.
3. Before M04 starts, resolve whether the `Q_u` regular-vs-singular gating
   observation above (Fork/InlineSameDeme reachability depending on whether
   `event.regular`) belongs in M04's scope or is still deferred further (to
   the driver/proposal milestone) — the milestone spec's own scoping
   suggests the latter, but this hasn't been explicitly re-confirmed since
   the observation was made during M03, after the spec was written.
4. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any
   M04 change and confirm the count does not regress below 3855/3855 (this
   milestone's new baseline).

## Next milestone
M04 — Explicit Reduction over the event indicator `m`: group `full_transitions`'s
output by reduced `d → d'` transition and compute `Φ_u(d,d') = Σ φ_u` over
the full transitions that collapse to it (KLI's `m`-marginalization, closed
form via Chu–Vandermonde — `mgp_filter.jl`'s `apply_move!` docstring cites
this as "Eq. 8", and `mers_filter_suite.tex` lines 430-432 as the general
statement). `IdentityTransition` and `InlineSameDemeTransition` are the
concrete case this milestone's IR was deliberately built to keep distinct
so M04 has the right building blocks to sum over.

## Context note
The most important thing to preserve if this conversation were compacted:
`src/examples/mgp_transitions.jl`'s `full_transitions(event, ℓ, n; Q=1)` is
now the canonical, tested, generic classifier from a saturation `s` (M02's
`enumerate_saturations`) to a full KLI coloring transition
(`IdentityTransition`/`InlineSameDemeTransition`/`CrossDemeTransition`/
`ForkTransition`), verified term-by-term against `mers_filter_suite.tex`'s
TCC/THH/THC/TCH derivations and independently cross-checked against the
tested, hand-coded `seir_naive.jl` for SEIR's `infection`/`progression`.
The migration-semantics question (whether `progression`'s `r=(0,1)`
encoding correctly implies a single CrossDeme/Identity choice, no
independent E-slot) is RESOLVED and CONFIRMED against `seir_naive.jl`, with
one correction to the milestone's own hypothesized closed form: the
CrossDeme case is `φ_u=1/I`, not `ℓ/I`. `Q_u`'s regular-vs-singular gating
(e.g. why `seir_naive.jl` never realizes `Fork` outside `singular_part!`)
was investigated, understood, and documented, but deliberately NOT
implemented — `full_transitions` returns the full ungated compatibility
structure per the milestone's explicit scope, and gating is left to a later
milestone. `ChopTransition` is a real type but a documented no-op stub;
`DEATH`/`SAMPLE`/`NEUTRAL` events are rejected by `full_transitions` based
on `event.type`, not on whether `r_u` happens to be zero (SEIR's `sampling`
event is a nonzero-`r` counterexample that still must be, and is, rejected).
