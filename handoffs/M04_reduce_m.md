# Handoff — M04: Explicit Reduction over the Event Indicator `m`

## Status
COMPLETE

## Objective
Turn M03's full (uncollapsed) `Vector{KLITransition}` for a single
event/state pair into the REDUCED, color-only transition table KLI actually
filters on: group the full transitions by the outcome they produce on the
`Coloring{D,N}` state (no event-indicator `m` component), and sum the exact
`Rational{Int}` `phi_u` compatibility factors within each group to obtain
`Phi_u`, per `mers_filter_suite.tex`'s "Coloring convention" paragraph and
`Phi_u`/boost definition (lines 417-432). Purely mechanical grouping — no
`Q_u` gating, no Filter IR classification (regular/singular/decay — M05), no
wiring into `mgp_filter.jl`'s stubs.

## What changed
- Added `src/examples/mgp_reduce.jl`: `ReducedTransition` (struct: `key`,
  `kind`, `Φ::Rational{Int}`, `transitions::Vector{KLITransition}`
  provenance), `reduced_key(t::KLITransition) -> Tuple` (dispatched per
  concrete subtype), `reduce_event_indicator(transitions) ->
  Vector{ReducedTransition}` (the grouping/summing pass), and the
  convenience composition `reduced_transitions(event, ℓ, n; Q=1) ->
  Vector{ReducedTransition}` (`full_transitions |> reduce_event_indicator`).
- Wired into `src/examples/Examples.jl`, immediately after
  `mgp_transitions.jl` and before `mgp_filter.jl` (which remains untouched).
- Added `test/kli_reduce_test.jl` (65 new `@test`s, registered in
  `test/runtests.jl` right after `kli_full_transitions_test.jl`): hand-summed
  `Rational{Int}` equality assertions for MERS TCC/THH/THC/TCH and SEIR
  `infection`/`progression`, all reusing M03's own already-verified
  `(ℓ,n)`/`phi` numbers unchanged; a provenance test (collapsed groups
  contain exactly the right transition *types*, not just counts); a
  structural confirmation that `progression` never produces
  `InlineSameDemeTransition`/`ForkTransition`; and generic partition-
  invariant checks (`Σ Φ == Σ φ`, every input transition appears in exactly
  one output group).

## Files added
- `src/examples/mgp_reduce.jl`
- `test/kli_reduce_test.jl`
- `handoffs/M04_reduce_m.md` (this file)

## Files modified
- `src/examples/Examples.jl` — one line added: `include("mgp_reduce.jl")`.
- `test/runtests.jl` — one line added: `include("kli_reduce_test.jl")`.
- `mgp_phi.jl` and `mgp_transitions.jl` untouched (consumed, not modified,
  per the milestone's explicit constraint). No filter module touched.

## Mathematical decisions

### 1. The collapsing rule, derived from `swap!`'s actual code, not assumed
`src/coloring.jl:36-45`'s `swap!(y, i, j, b)` is `delete!(y[i], b); push!(y[j],
b)`. When `i == j` (`InlineSameDemeTransition`'s case, since its slot's
post-event deme always equals `event.from` per M03's Step-A-derived
classification), this is `delete!` then `push!` of the *same* lineage id into
the *same* deme's `BitSet` — a genuine no-op on `y.cols`, confirmed by
reading the code rather than assuming the tex's informal "`σ^b_{CC}y=y`" is
literally what the data structure does. This is the basis for collapsing
`IdentityTransition` (no slot occupied) and `InlineSameDemeTransition` (one
slot occupied, ancestral==post-event deme) into a single `(:noop, d)` reduced
class, `d = event.from` (verified identical for both — `InlineSameDemeTransition.deme`
is by construction always `== event.from`, since M03's classifier only
produces that subtype when the slot's deme equals `event.from`).

### 2. `CrossDemeTransition` and `ForkTransition` never collapse with noop or each other
A `CrossDemeTransition` moves a tracked lineage between two *different* demes
on the reduced `Coloring` (`swap!` with `i != j` is not a no-op — it changes
which `BitSet` contains `b`), a genuinely different reduced outcome from "no
change". Keyed by `(:cross, ancestral_deme, target_deme)`, so distinct
deme-pairs never collapse into each other either. A `ForkTransition` *adds* a
lineage to the coloring (`fork!` pushes into two target `BitSet`s after
deleting the ancestor) — again not reducible to "no change" or a "move".
Keyed generically by `(:fork, ancestral_deme, sorted_slot_demes_tuple)`, so a
hypothetical future event with two *different* Fork saturations landing on
the identical `(ancestral_deme, slot_demes multiset)` would still correctly
collapse — this does not occur in any current SEIR/MERS event (verified: TCC
and THH each have exactly one Fork-classified saturation, since their only
2-slot saturation is `s=r_u` in the single nonzero deme; THC/TCH's only
2-slot saturation is `(1,1)`), so the generic key was exercised only as a
singleton group in every test here, per the milestone's own framing
("implement the grouping key generically... this does not occur for any
current SEIR/MERS event").

### 3. `SEIR progression` is a genuine "nothing collapses" contrast case
`progression`'s production vector is `r=(0,1)` — a single slot that always
lives at the TARGET deme (I), never the ancestral deme (E) (M03's migration-
semantics finding, `handoffs/M03_full_kli_lowering.md` §1). Consequently
`sum(s) ∈ {0,1}` only (never `≥2`, so `ForkTransition` is structurally
impossible — there is only one slot, full stop) and the one occupied slot's
deme is always `I ≠ E = event.from`, so `InlineSameDemeTransition` is also
structurally impossible (M03's classifier only emits it when the occupied
slot's deme equals `event.from`). `test/kli_reduce_test.jl`'s `"SEIR
progression"` testset asserts this structurally (`all(t isa
Union{IdentityTransition,CrossDemeTransition} for t in ts)`), not just
numerically — confirming the "2 full → 2 reduced, no collapse" case rather
than assuming it.

## Architecture decisions
- `reduce_event_indicator` takes `AbstractVector{<:KLITransition}` (not
  restricted to the concrete `Vector{KLITransition}` M03 happens to return),
  so it composes with any iterable of `KLITransition`s a caller might
  construct.
- `ReducedTransition.key` is an opaque generic `Tuple`, distinguishable
  structurally by its first element (`:noop`/`:cross`/`:fork`) and remaining
  fields (deme indices, or a sorted `Tuple` of slot demes for Fork) — chosen
  as a `Tuple` (not a custom key struct) because it needs no behavior beyond
  `==`/`hash`, both of which `Tuple` gets for free, including for the nested
  `Tuple` used inside Fork's key (avoided a `Vector` there specifically to
  keep the key immutable/hashable without caveats).
- `ReducedTransition.kind` is a redundant-but-convenient `Symbol` mirror of
  `key[1]`, added because `test/kli_reduce_test.jl` (and, per the spec,
  later audit tooling) frequently wants to dispatch/group by kind without
  unpacking the key tuple.
- Provenance (`ReducedTransition.transitions`) is a `Vector{KLITransition}`
  in first-encountered input order — deliberately NOT discarded, per the
  milestone's explicit instruction that later audit/explain milestones need
  it.
- Default Julia struct `==` is identity-based (`===`) for structs carrying
  `Vector` fields (confirmed by direct test — see "Tests run" below) — this
  is why `test/kli_reduce_test.jl` compares `ReducedTransition`s via
  `Set((rt.kind, rt.Φ) for rt in ...)` rather than `==` on the structs
  themselves; noted here so a future milestone doesn't assume
  `ReducedTransition` supports value equality out of the box.
- `reduce_event_indicator` is `AbstractVector{<:KLITransition}` in, plain
  `Vector{ReducedTransition}` out — grouping order is "first key
  encountered in the input", not sorted/canonicalized, since no test or
  downstream consumer needs a canonical order (callers needing a specific
  kind look it up, e.g. via `Dict(rt.kind => rt for rt in rts)` as the tests
  do).

## Full → reduced → Φ tables (all six required marks)

**TCC** (`r=(2,0)`, `I_C=5, ℓ_C=2` — instance from `full_transitions(tcc,
[2,0],[5,3])`):
```
TCC: s=0 (Identity, φ=3/10) + s=1 (InlineSameDeme, deme C, φ=3/10) → noop-in-C → Φ = 3/10+3/10 = 3/5
TCC: s=2 (Fork, C+C, φ=1/10)                                       → fork-in-C  → Φ = 1/10
```

**THH** (`r=(0,2)`, `I_H=6, ℓ_H=2` — instance from `full_transitions(thh,
[0,2],[3,6])`):
```
THH: s=0 (Identity, φ=2/5) + s=1 (InlineSameDeme, deme H, φ=4/15) → noop-in-H → Φ = 2/5+4/15 = 2/3
THH: s=2 (Fork, H+H, φ=1/15)                                      → fork-in-H  → Φ = 1/15
```

**THC** (`r=(1,1)`, `I_C=5,ℓ_C=2,I_H=4,ℓ_H=1` — instance from
`full_transitions(thc, [2,1],[5,4])`):
```
THC: s=(0,0) (Identity, φ=9/20) + s=(1,0) (InlineSameDeme, deme C, φ=3/20) → noop-in-C     → Φ = 9/20+3/20 = 3/5
THC: s=(0,1) (CrossDeme C→H, φ=3/20)                                       → cross C→H     → Φ = 3/20
THC: s=(1,1) (Fork, C+H, φ=1/20)                                           → fork (C anc.) → Φ = 1/20
```
(4 full → 3 reduced; the (1,0) InlineSameDeme collapses with the (0,0)
Identity, NOT with the (0,1) CrossDeme — confirmed by assertion, not
assumed.)

**TCH** (mirror of THC, `I_H=5,ℓ_H=2,I_C=4,ℓ_C=1` via the same
`(2,1),(5,4)` call with roles swapped):
```
TCH: s=(0,0) (Identity, φ=9/20) + s=(0,1) (InlineSameDeme, deme H, φ=3/20) → noop-in-H     → Φ = 9/20+3/20 = 3/5
TCH: s=(1,0) (CrossDeme H→C, φ=3/20)                                       → cross H→C     → Φ = 3/20
TCH: s=(1,1) (Fork, H+C, φ=1/20)                                           → fork (H anc.) → Φ = 1/20
```

**SEIR `infection`** (`r=(1,1)`, parent deme I, `n_E=6,ℓ_E=2,n_I=5,ℓ_I=2` —
instance from `full_transitions(infection, [2,2],[6,5])`):
```
infection: s=(0,0) (Identity, φ=2/5) + s=(0,1) (InlineSameDeme, deme I, φ=2/15) → noop-in-I    → Φ = 2/5+2/15 = 8/15
infection: s=(1,0) (CrossDeme I→E, φ=1/10)                                      → cross I→E    → Φ = 1/10
infection: s=(1,1) (Fork, E+I, φ=1/30)                                          → fork (I anc.) → Φ = 1/30
```

**SEIR `progression`** (`r=(0,1)`, MIGRATION, `n_I=5,ℓ_I=2` — instance from
`full_transitions(progression, [3,2],[9,5])`):
```
progression: s=0 (Identity, φ=3/5) → noop-in-E  → Φ = 3/5   (singleton, no collapse)
progression: s=1 (CrossDeme E→I, φ=1/5) → cross E→I → Φ = 1/5   (singleton, no collapse)
```
2 full → 2 reduced, unchanged — the contrast case: no `InlineSameDemeTransition`
or `ForkTransition` is even structurally possible for this event (single
slot, always at the target deme).

## Chu–Vandermonde closed-form note (optional cross-check, done)
For TCC, `Φ(noop) = φ(0)+φ(1) = C(I_C-ℓ_C,2)/C(I_C,2) + 2(I_C-ℓ_C)/[I_C(I_C-1)]`.
Combining over the common denominator `I_C(I_C-1)`:
`Φ(noop) = [(I_C-ℓ_C)(I_C-ℓ_C-1) + 2(I_C-ℓ_C)] / [I_C(I_C-1)]
         = (I_C-ℓ_C)(I_C-ℓ_C+1) / [I_C(I_C-1)]`.
At `I_C=5, ℓ_C=2`: `(3)(4)/(5·4) = 12/20 = 3/5` — matches the direct sum
exactly. This is the Chu–Vandermonde-flavored simplification the tex alludes
to (summing binomial-ratio terms over adjacent `s` telescopes into a single
binomial-difference form); confirmed algebraically for TCC only, not
re-derived for the other five marks (the direct-summation match is the
acceptance gate; this closed form is a bonus cross-check, not re-verified
elsewhere to avoid excessive symbolic time as the milestone spec allows).

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before (M03 baseline): **3855/3855 passed**.
- After this milestone: **3920/3920 passed** (`3855 + 65`), ~42s wall time,
  0 failures, 0 errors. The delta is exactly the new `Explicit reduction
  over m (M04)` testset in `test/kli_reduce_test.jl`. No pre-existing
  testset's pass count changed.
- One test-authoring finding during this milestone (documented in
  "Architecture decisions" above, not a bug in `mgp_reduce.jl`): confirmed
  by a standalone Julia snippet that default struct `==` for a type with a
  `Vector` field is identity-based (`===`), not field-wise — adjusted the
  two convenience-wrapper-agreement assertions in `test/kli_reduce_test.jl`
  to compare `(kind, Φ)` sets instead of raw `==` on `ReducedTransition`
  vectors.

## Verification status
| Layer | Status | Oracle |
|---|---|---|
| `reduce_event_indicator` generic grouping/summing (no per-event lookup table) | Verified (new) | `test/kli_reduce_test.jl`, all testsets |
| TCC: 3 full → 2 reduced, `Φ(noop)=φ(0)+φ(1)`, `Φ(fork)=φ(2)` alone | Verified (new) | hand-summed from M03's own verified `phi` values; `test/kli_reduce_test.jl` `"TCC"` |
| THH: mirror of TCC | Verified (new) | `"THH"` |
| THC: 4 full → 3 reduced, critical (1,0)-collapses-with-(0,0)-not-(0,1) distinction | Verified (new) | `"THC"`, including a `key`-inequality assertion between the noop and cross groups |
| TCH: mirror of THC | Verified (new) | `"TCH"` |
| SEIR `infection`: 4 full → 3 reduced, same shape as THC | Verified (new) | `"SEIR infection"` |
| SEIR `progression`: 2 full → 2 reduced, structurally confirmed no-collapse | Verified (new) | `"SEIR progression"`, includes a structural (not just numeric) assertion that Inline/Fork are impossible |
| Provenance (collapsed groups list the right transition *types*) | Verified (new) | `"TCC"`'s `noop_types == Set([IdentityTransition, InlineSameDemeTransition])` assertion |
| Partition invariants (`Σ Φ == Σ φ`, every input in exactly one group) | Verified (new) | `"invariants"` testset, `objectid`-based set equality |
| Chu–Vandermonde closed form (TCC only, optional) | Verified by hand algebra (not a `@test`) | this file, "Chu–Vandermonde closed-form note" |
| `Q_u` regular-vs-singular gating | NOT implemented — explicitly out of scope, deferred further (M05/driver) | N/A |
| Filter IR classification (regular/singular/decay) | NOT implemented — explicitly out of scope | N/A, M05 |

## Known failures / unresolved issues
- None outstanding. No test-authoring bug survived past the equality-fallback
  finding above (fixed within this milestone).
- Judgment call carried over, unresolved, from M03 (not re-litigated here,
  since M04's spec explicitly kept `Q_u` gating out of scope): whether
  `seir_naive.jl`'s regular-event proposal choosing CrossDeme unconditionally
  (never realizing the InlineSameDeme/Fork alternatives this milestone's
  `Φ_u` sums include) reflects a validated importance-sampling design is
  still a filter-correctness question for a later milestone, not this one.
  Flagging again only because M05 (Filter IR) is exactly where this becomes
  load-bearing: if `Φ_u`'s InlineSameDeme contribution corresponds to a
  saturation `seir_naive.jl` never proposes at regular events, M05 needs to
  decide whether that piece of `Φ_u` belongs to `RegularFlow` at all, or
  only ever contributes at singular (observed) events.

## Git state
- branch: `atpabuser-devel`
- commit if created: none
- uncommitted files (cumulative, M00-M04; nothing committed or pushed by any
  milestone so far, per instructions):
  - `docs/compiler/architecture_before.md`, `docs/compiler/compiler_roadmap.md` (M00)
  - `handoffs/M00_reconnaissance.md` .. `handoffs/M03_full_kli_lowering.md` (M00-M03)
  - `src/examples/mgp_audit.jl`, `test/population_ir_test.jl` (M01)
  - `src/examples/mgp_phi.jl`, `test/kli_phi_test.jl` (M02)
  - `src/examples/mgp_transitions.jl`, `test/kli_full_transitions_test.jl` (M03)
  - `handoffs/M04_reduce_m.md` (new, this file)
  - `src/examples/mgp_reduce.jl` (new)
  - `test/kli_reduce_test.jl` (new)
  - `src/examples/Examples.jl` (modified — one more `include` line added)
  - `test/runtests.jl` (modified — one more `include` line added)

## Resume instructions
1. Advisor reviews `src/examples/mgp_reduce.jl` and
   `test/kli_reduce_test.jl`, in particular: (a) the collapsing-rule
   justification traced to `swap!`'s actual `BitSet` semantics (not just the
   tex's informal notation), (b) the THC critical-distinction test (`(1,0)`
   collapses with `(0,0)`, not `(0,1)`), and (c) the generic Fork key design
   (never exercised as a real multi-member group by current SEIR/MERS, per
   design).
2. M05 ("Filter IR — RegularFlow/SingularFlow/DecayTerm classification")
   should treat `reduce_event_indicator`'s `Vector{ReducedTransition}` as its
   input: each `ReducedTransition` becomes a candidate Filter IR term, and
   M05 must decide, per term, whether it is reachable as `RegularFlow`
   (continuous-time, unobserved-interval proposal) or only as `SingularFlow`
   (at a fixed/observed data-event time) or contributes to `DecayTerm` — this
   is exactly the `Q_u` regular-vs-singular gating M03 found and both M03
   and M04 deliberately left unimplemented. In particular: M03 found that
   `seir_naive.jl`'s `regular_part!` never realizes `Fork` outside
   `singular_part!` (see M03's handoff "Known failures"); M04 additionally
   found `seir_naive.jl`'s regular step also never realizes
   `InlineSameDemeTransition` for `infection` (always chooses CrossDeme) —
   M05 needs to work out, from the tex or from filter-correctness principles,
   whether that means the `(:noop,...)` `ReducedTransition`'s `Φ` (which
   includes the InlineSameDeme contribution, per M04's collapsing) is fully
   `RegularFlow`-eligible or whether some sub-decomposition survives into
   M05's classification.
3. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any M05
   change and confirm the count does not regress below 3920/3920 (this
   milestone's new baseline).

## Next milestone
M05 — Filter IR: wrap each `ReducedTransition` (this milestone's output) into
a classified term — `RegularFlow` (continuous, unobserved-interval,
`Q_u`-reachable), `SingularFlow` (fixed/observed data-event time only), or
`DecayTerm` (the outflow/rate-driven pieces M03 deferred entirely, e.g.
removal's `I_d-1≥ℓ_d` condition and sampling's background-rate contribution,
`mers_filter_suite.tex:511-539`) — with provenance threaded through from
`ReducedTransition.transitions` (and, transitively, M03's `KLITransition`
provenance). M03 already found, and M04 corroborated, that `Fork` and (per
M04's new finding) `InlineSameDeme`-derived reduced transitions are not
realized by `seir_naive.jl`'s regular-event proposal — M05 needs to formalize
why (via `Q_u`) rather than leave it as an observation in a handoff doc.

## Context note
The most important thing to preserve if this conversation were compacted:
`src/examples/mgp_reduce.jl`'s `reduce_event_indicator(transitions) ->
Vector{ReducedTransition}` is now the canonical, tested, generic reduction
pass from M03's full `KLITransition`s to KLI's reduced, color-only `Φ_u`
table. The collapsing rule is: `IdentityTransition` and
`InlineSameDemeTransition` ALWAYS collapse into one `(:noop, event.from)`
group (verified via `swap!`'s literal same-deme-no-op code, not just the
tex's informal notation); `CrossDemeTransition` and `ForkTransition` NEVER
collapse with the noop group or (for distinct deme-pairs/slot-deme-multisets)
with each other. All six required marks (TCC/THH/THC/TCH/SEIR
infection/progression) were hand-verified against M03's own already-verified
`phi` numbers — no new `phi` values were computed, only their grouped sums.
`ReducedTransition` keeps full provenance (`transitions::Vector{KLITransition}`)
deliberately, for later audit tooling. Default Julia struct `==` on
`ReducedTransition` is identity-based (contains a `Vector` field), NOT
field-wise — any future code comparing `ReducedTransition`s for value
equality must compare fields explicitly (e.g. `(kind, Φ)` tuples), not rely
on `==`. `Q_u` regular-vs-singular gating remains deliberately unimplemented,
now with a second concrete data point (`InlineSameDeme` also never realized
by `seir_naive.jl`'s regular step, not just `Fork`) documented for M05.
