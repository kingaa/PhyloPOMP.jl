# Handoff — M01: Population IR

## Status
COMPLETE

## Objective
Establish a clean, mechanically-inspectable Population IR on top of the
existing `Event`/`MGPModel` scaffold, without changing any likelihood
math: confirm the IR is structurally complete for both `SEIR` and `MERS`,
and add a callable `audit_model`/`validate_model` pair so the static model
representation can be inspected and semantically checked outside of
`@mgp`'s macroexpansion-time `error()` calls.

## What changed
- Evaluated whether `Event` should become a typed-subtype hierarchy
  (`BirthEvent`, `MigrationEvent`, ...) for future multiple-dispatch in
  M02/M03. Decided **no** — left `Event`/`EventType` structurally
  unchanged. See "Architecture decisions" below for the reasoning.
- Added `test/population_ir_test.jl`: walks every event of `SEIR` and
  `MERS` and asserts Δ_u, α_u, r_u, and W_u (from/into) are all present
  and well-typed; evaluates every hazard closure at a real state/parameter
  point and checks the result is finite and non-negative; cross-checks
  `audit_model(SEIR)` field-by-field against `SEIR_REFERENCE`; checks
  `audit_model(MERS)`'s deme-index invariants; checks `validate_model`
  returns no issues for `SEIR`, `MERS`, and `SEIR_REFERENCE`, and *does*
  report issues for a deliberately broken hand-built `MGPModel`; smoke
  tests the `show` output. Registered in `test/runtests.jl` immediately
  after `seir_macro_equivalence.jl`.
- Added `src/examples/mgp_audit.jl`:
  - `audit_model(model::MGPModel) -> ModelAudit` — ordinary function (not
    a macro), returns a structured `ModelAudit`/`EventAudit` value
    exposing compartments, demes, and per-event semantic type,
    regular/singular, observed, Δ, hazard-presence, r, and from/into
    wiring (resolved to both raw indices and deme names). `Base.show`
    (`MIME"text/plain"`) renders a human-readable table; verified at the
    REPL against both `SEIR` and `MERS`.
  - `validate_model(model::MGPModel) -> Vector{String}` — separately
    callable, non-throwing semantic validation covering the checks M00
    found folded into `@mgp`'s macroexpansion (duplicate names, deme ⊆
    compartments) plus checks that did not exist anywhere before this
    milestone (every event's `from`/`into` indices valid against
    `model.demes`, `r` length matching `length(model.demes)`, hazard
    callability, Δ referencing declared compartments, per-`EventType`
    wiring sanity). Implemented as the Step 5 stretch goal — turned out
    clean and additive, so it's included as COMPLETE, not partial.
- Wired `mgp_audit.jl` into the include chain in
  `src/examples/Examples.jl`, between `mgp.jl` (defines `Event`/
  `MGPModel`) and `mgp_filter.jl` (untouched).

## Files added
- `src/examples/mgp_audit.jl`
- `test/population_ir_test.jl`
- `handoffs/M01_population_ir.md` (this file)

## Files modified
- `src/examples/Examples.jl` — one line added: `include("mgp_audit.jl")`.
- `test/runtests.jl` — one line added: `include("population_ir_test.jl")`.
- No existing struct, macro, filter, or likelihood-bearing file touched.
  `Event`, `MGPModel`, `mgp_macro.jl`, `mgp_mers.jl`, `mgp_filter.jl`, and
  all eight hand-coded filter modules are byte-for-byte unchanged.

## Mathematical decisions
- None. This milestone is purely structural/inspective — no φ_u, Φ_u,
  Q_u, or saturation logic was written, per the milestone's explicit
  constraint. The only "convention chosen" is naming: `EventAudit`/
  `ModelAudit` field names (`delta`, `has_hazard`, `r`, `from`/`into`)
  intentionally echo the Δ_u/α_u/r_u/W_u vocabulary from `Event`'s
  docstring (`src/examples/mgp.jl:34-58`) and
  `docs/compiler/compiler_roadmap.md` item 2, so later milestones/docs can
  refer to both consistently.

## Architecture decisions
- **Typed events (Step 1): kept the single `Event` struct + `EventType`
  enum, did not refactor into a `BirthEvent <: MGPEvent` subtype
  hierarchy.** Reasoning:
  - Every place that will eventually need to dispatch on event semantics
    is a small, closed set: `kli_select`, `kli_decay`, `apply_move!`,
    `singular_update!` (all stubs in `mgp_filter.jl`) plus the already-
    working `kli_hazard`/`apply_pop` (which are semantics-agnostic and
    would gain nothing from subtyping). That's at most 4 functions that
    will ever need an `if ev.type == ...`/`@match` over 5 cases — roughly
    20 branches total across the whole future filter core, not a
    combinatorial dispatch problem.
  - A subtype hierarchy would require either (a) duplicating the 8 common
    fields (`name`, `Δ`, `hazard`, `r`, `from`, `into`, `regular`,
    `observed`) across 5 concrete structs, or (b) an abstract supertype
    with accessor functions replicating what a `NamedTuple`-backed single
    struct already gives for free — added ceremony with no dispatch
    benefit until a function actually needs *different fields* per event
    type (none currently do; `from`/`into`/`r` are already uniformly
    present and just semantically unused for e.g. `NEUTRAL`/`DEATH`).
  - It would also directly threaten the one hard invariant this milestone
    was told not to break: `test/seir_macro_equivalence.jl`'s exact
    struct-field equality between the `@mgp`-generated `SEIR.events` and
    the hand-written `SEIR_REFERENCE` literal table
    (`src/examples/mgp.jl:91-107`) — changing `Event`'s type would force
    rewriting that oracle, which is exactly the kind of invasive,
    unforced change the milestone spec warned against.
  - Concrete evidence against redesign: 2 models, 5 event types, no
    function today branches on `ev.type` at all (`kli_hazard`/`apply_pop`
    are structurally generic), and the codebase is small enough (per M00)
    that a future genuine need for dispatch can be revisited in M02/M03
    once `kli_select`/`apply_move!` are actually being filled in and it's
    clear whether their per-type logic is naturally expressed as
    `if`/`elseif` (cheap) or wants real multiple dispatch (then worth the
    refactor, with a real reason instead of a speculative one).
  - Documented, not silently decided: this reasoning lives here and should
    be the first thing M02 reconsiders if `kli_select`/`apply_move!`
    implementations start feeling awkward under the enum-tag approach.
- `audit_model`/`validate_model` were placed in a new file
  (`src/examples/mgp_audit.jl`) rather than appended to `mgp.jl`, to keep
  `mgp.jl` scoped to "model layer is DATA" (per its own header comment)
  and keep the audit/validation *tooling* separable — this also sets up
  cleanly for a later `@mgpaudit` macro (mentioned in the milestone spec
  as landing in M06) to be a thin formatting wrapper around
  `audit_model`, per the project's macros-are-syntax-not-architecture
  rule.
- `validate_model` returns `Vector{String}` rather than a richer `Issue`
  struct (with severity/field/event provenance) — kept intentionally
  simple for M01. If M02+ wants machine-actionable validation results
  (e.g. to gate `compile_filter`), revisit this as a small `struct Issue`
  with a `severity`/`message`/`event` triple; not needed yet since nothing
  currently consumes `validate_model`'s output programmatically.
- Did not attempt to split `parse_model` from `lower_population` (the
  design fork M00 flagged) — out of scope for this milestone's stated
  task list, which treated the DSL/lowering fusion as already-acceptable
  and focused this milestone on the *output* (`Event`/`MGPModel`)
  instead. Left as an open fork for the advisor per M00's handoff.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before (M00 baseline, confirmed in `handoffs/M00_reconnaissance.md`):
  **3248/3248 passed**.
- After (this milestone's changes applied): **3596/3596 passed**, ~39s
  wall time, 0 failures, 0 errors. The delta (+348) is exactly the new
  `Population IR structural audit` testset in
  `test/population_ir_test.jl` (visible as its own line in the per-testset
  summary: `Population IR structural audit | 348 348 0.5s`). No
  pre-existing testset's pass count changed.

## Verification status
| Layer | Status | Oracle |
|---|---|---|
| Structure — Δ_u/α_u/r_u/W_u present & well-typed for every SEIR/MERS event | Verified (new) | `test/population_ir_test.jl` |
| Structure — hazards evaluate to finite, non-negative reals | Verified (new) | `test/population_ir_test.jl` (real x/θ per model) |
| `audit_model(SEIR)` vs. `SEIR_REFERENCE` | Verified (new) | `test/population_ir_test.jl` field-by-field comparison |
| `audit_model(MERS)` deme-index invariants | Verified (new) | `test/population_ir_test.jl`; also eyeballed against `mgp_mers.jl:3-24` at the REPL (see below) |
| `validate_model` clean on SEIR/MERS/SEIR_REFERENCE | Verified (new) | `test/population_ir_test.jl` |
| `validate_model` catches broken models | Verified (new) | `test/population_ir_test.jl` deliberately-broken `MGPModel` fixture |
| `show`/REPL readability | Manually verified | `julia --project=. -e 'using PhyloPOMP; show(stdout, MIME("text/plain"), audit_model(PhyloPOMP.MERS))'` — output eyeballed to match `mgp_mers.jl:3-24` event-by-event (self-loop r=(2,0)/(0,2) for same-species transmission, r=(1,1) for cross-species, chop/sample_remove wiring, `none`-move demographic events) |
| Existing DSL→Event lowering (`seir_macro_equivalence.jl`) | Unaffected — still passing | `test/seir_macro_equivalence.jl` (unchanged, still 3006/3006) |
| Filter math / KLI weighting | Not touched this milestone | N/A — explicitly out of scope |

## Known failures / unresolved issues
- None introduced by this milestone. Everything from M00's "Known
  failures / unresolved issues" list remains open and unaddressed here
  (by design — out of scope): the `parse_model`/`lower_population` split
  question, the proposal-taxonomy formalization question, and the
  `enumerate_saturations`/`derive_full_compatibility`/`compute_phi`/
  `reduce_event_indicator` mathematical gap. See
  `handoffs/M00_reconnaissance.md` for the full list; M01 did not
  resolve or touch any of it.
- `validate_model`'s per-`EventType` wiring checks are deliberately
  minimal (only "does a from/into exist where the type requires one") —
  they do not check e.g. that a `BIRTH` event's `into` is non-empty, or
  that `r` is consistent with `into` counts. Flagging as a possible
  follow-up if M02 wants stronger IR-level guarantees before building
  `enumerate_saturations`/`compute_phi` against it, but not filed as a
  blocking gap for this milestone.

## Git state
- branch: `atpabuser-devel`
- commit if created: none
- uncommitted files:
  - `docs/compiler/architecture_before.md` (from M00, still uncommitted)
  - `docs/compiler/compiler_roadmap.md` (from M00, still uncommitted)
  - `handoffs/M00_reconnaissance.md` (from M00, still uncommitted)
  - `handoffs/M01_population_ir.md` (new, this file)
  - `src/examples/mgp_audit.jl` (new)
  - `test/population_ir_test.jl` (new)
  - `src/examples/Examples.jl` (modified — one `include` line added)
  - `test/runtests.jl` (modified — one `include` line added)
  - No file was committed or pushed by this milestone, per instructions.

## Resume instructions
1. Advisor reviews `src/examples/mgp_audit.jl` and
   `test/population_ir_test.jl`, and in particular the Step 1
   typed-events decision above — confirm agreement before M02 starts
   building `kli_select`/`apply_move!` logic against the existing
   enum-tag `Event` (rather than assuming subtypes will appear later).
2. M02 ("Generic Production Slots, Saturations, and φ") should treat
   `audit_model`'s `EventAudit`/`ModelAudit` structures as a convenient
   read-only inspection tool while implementing
   `enumerate_saturations`/`compute_phi` — not as something that itself
   needs to grow KLI math; keep the audit layer static-IR-only.
3. Before M02 starts hand-deriving `φ_u`/`Φ_u` generically, resolve M00's
   open fork on how to bring `StructuredMGPs.pdf`'s formulas into the
   repo as tracked, citable source material (excerpt vs. full PDF) — this
   milestone did not touch that question.
4. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any
   M02 change and confirm the count does not regress below 3596/3596
   (this milestone's new baseline).

## Next milestone
M02 — Generic Production Slots, Saturations, and φ

## Context note
The most important thing to preserve if this conversation were compacted:
`src/examples/mgp_audit.jl`'s `audit_model`/`validate_model` are now the
canonical, tested way to inspect and semantically check an `MGPModel`
without touching the macro or the filter stubs — use them (don't
reinvent) when M02+ needs to iterate over a model's events
programmatically. The Step 1 decision to *not* introduce typed `Event`
subtypes was deliberate and reasoned (see "Architecture decisions"), not
an oversight — if M02's `kli_select`/`apply_move!` implementations end up
wanting real multiple dispatch on event semantics, that is the first place
to revisit this call, with concrete evidence from actually writing that
code, rather than reopening it speculatively.
