# Handoff — M02: KLI φ_u

## Status
COMPLETE

## Objective
Implement, as ordinary testable Julia functions (not macros), the two
genuinely-missing generic pieces flagged by M00/M01 as the compiler's
mathematical core: `enumerate_saturations` (the saturation space `S_u(ℓ)`
for a production vector and lineage-count state) and `kli_binomial_ratio`
(the KLI production-slot compatibility ratio `φ_u`), computed from the
Population IR's `Event.r` field alone, with `Q_u` left as an external
parameter (Q_u derivation is M03's job). No hand-coded filter module was
touched; no Q_u derivation, coloring-operator logic, or marginalization
over `m` (Φ_u) was implemented — those are explicitly out of scope per the
milestone spec.

## What changed
- Added `production_slots(event)`, `enumerate_saturations(r, ℓ)` /
  `enumerate_saturations(event, ℓ)`, and `kli_binomial_ratio(r, s, ℓ, n;
  Q=1)` / `kli_binomial_ratio(event, s, ℓ, n; Q=1)` to a new file
  `src/examples/mgp_phi.jl`, wired into the include chain in
  `src/examples/Examples.jl` immediately after `mgp_audit.jl` (before
  `mgp_filter.jl`, which remains untouched).
- `production_slots` finally reads `Event.r` — M00 found it populated by
  the DSL but never consumed by any downstream code; this is the first
  consumer.
- `enumerate_saturations` is a genuine nested/product iteration over each
  deme's feasible range `0:min(r_d, ℓ_d)` (via `Iterators.product`), not a
  per-model/per-event literal list — verified directly against `Event`
  objects pulled at runtime from the real `SEIR`/`MERS` `MGPModel`s (not
  hand-typed duplicates of their fields), satisfying the M02 acceptance
  gate's "no hard-coded saturation list" requirement.
- `kli_binomial_ratio` computes the exact formula
  `φ_u(s) = Q · ∏_d C(n_d-ℓ_d, r_d-s_d) / C(n_d, r_d)` using a custom
  `safe_binomial` wrapper (see "Mathematical decisions") and returns an
  exact `Rational{Int}`, never `Float64`.
- Added `test/kli_phi_test.jl` (181 new `@test`s, registered in
  `test/runtests.jl` right after `population_ir_test.jl`) with hand-derived
  reference values for every case the milestone spec required, plus
  cross-checks against real `SEIR`/`MERS` events and `mers_filter_suite.tex`.

## Files added
- `src/examples/mgp_phi.jl`
- `test/kli_phi_test.jl`
- `handoffs/M02_kli_phi.md` (this file)

## Files modified
- `src/examples/Examples.jl` — one line added: `include("mgp_phi.jl")`.
- `test/runtests.jl` — one line added: `include("kli_phi_test.jl")`.
- No existing struct, macro, filter, or likelihood-bearing file touched.
  `Event`, `MGPModel`, `mgp_macro.jl`, `mgp_mers.jl`, `mgp_filter.jl`,
  `mgp_audit.jl`, and all eight hand-coded filter modules are byte-for-byte
  unchanged.

## Mathematical decisions
- **φ_u formula**: implemented exactly as given in the M02 task spec
  (already vetted against the KLI paper by the project supervisor), and
  cross-checked term-by-term against `mers_filter_suite.tex`'s Step C
  (lines 450-457: `φ_u = ∏_d C(I_d-ℓ_d, r_d-s_d)/C(I_d, r_d)`, where the
  tex's `I_d` is exactly this spec's post-event `n_d`).
- **`C(a,b)` zero convention**: `Base.binomial(a,b)` in Julia 1.12 already
  returns `0` for `b < 0` and for `0 <= a < b` — verified with `@assert`s
  at module load time (`binomial(5,-1)==0`, `binomial(5,6)==0`) rather than
  assumed. It does **not** return `0` for `a < 0` — it falls back to the
  generalized/Pascal binomial-coefficient extension, e.g.
  `binomial(-1,1) == -1`, confirmed by `@assert binomial(-1,1) == -1` in
  `mgp_phi.jl`. Since `a` here is always `n_d - ℓ_d` and the model
  invariant is `ℓ_d <= n_d`, `a < 0` should not arise from valid input —
  but I wrote a `safe_binomial(a,b)` wrapper that explicitly zeroes
  `a < 0 || b < 0 || b > a` anyway, so the function degrades gracefully
  instead of silently returning a wrong negative "count" if that invariant
  is ever violated by a future caller.
- **`n_d < r_d` (denominator `C(n_d,r_d) = 0`)**: rather than dividing by
  zero, `kli_binomial_ratio` returns exactly `0` for the whole product in
  that case. Documented as a defensive convention, not a mathematical
  derivation — the *model-level* invariant is that `n_d >= r_d` always
  holds post-event (production already happened), so this branch should
  not be hit by valid M03+ callers, but the function does not crash if it
  is. Tested explicitly (`test/kli_phi_test.jl`, "boundary n<r").
- **Negative `s_d` (found during test-writing, not assumed in advance)**:
  I initially assumed `safe_binomial`'s `b<0` rule alone would zero out any
  `s` outside `[0, min(r,ℓ)]`, including `s_d < 0`. This is **false**:
  `s_d < 0` makes the numerator's `b = r_d - s_d` *larger* than `r_d`, not
  negative, so e.g. `kli_binomial_ratio([1],[-1],[5],[10])` originally
  returned `1//1`, not `0`. Caught by my own test suite (see "Tests run"
  below — this was a real, not merely hypothetical, bug found before the
  handoff was written). Fixed by adding an explicit `any(s_d < 0)` guard at
  the top of `kli_binomial_ratio` that returns `0` immediately. `s_d > r_d`
  remains correctly caught by the pre-existing `b<0` rule (numerator
  argument `r_d - s_d < 0`). This is documented in the function's
  docstring as a deliberate, tested design choice, not left implicit.
- **`Q_u` handling**: accepted as an external `Real` keyword, default `1`,
  multiplied into the `Rational{Int}` product. No derivation attempted —
  explicitly M03's job per the milestone constraints.

## Architecture decisions
- Two-arity design for both `enumerate_saturations` and
  `kli_binomial_ratio`: a low-level form taking raw `r`/`s`/`ℓ`/`n` vectors
  (matches the spec's exact signatures, and is what all the hand-computed
  reference-value tests exercise directly, without needing to construct a
  dummy `Event`), plus a convenience overload taking an `Event` (which
  simply calls `production_slots(event)` and delegates). This lets the
  same tests both assert exact hand-derived numbers *and* cross-check the
  generic function against real, unmodified `SEIR`/`MERS` `Event` objects
  pulled by name from `model.events` at runtime.
- Return type of `kli_binomial_ratio` is always `Rational{Int}` (never
  `Float64`), per the spec's "exact arithmetic, not floating point"
  instruction — makes every boundary case (`0`, `1`, `1/20`, etc.) an
  exact `==` comparison in tests, not an approximate one.
- Placed in a new file (`mgp_phi.jl`) between `mgp_audit.jl` and
  `mgp_filter.jl` in the include chain, per the milestone spec — keeps
  `mgp_filter.jl`'s stubs (`kli_select`/`kli_decay`/`apply_move!`/
  `singular_update!`) as the landing spot for M03+ code that will *consume*
  `enumerate_saturations`/`kli_binomial_ratio`, without touching those
  stubs yet.

## Tests run
- `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl`
- Before (M01 baseline, confirmed in `handoffs/M01_population_ir.md`):
  **3596/3596 passed**.
- First run of this milestone's new suite: **1 failure** out of 181 new
  tests — the negative-`s` defensiveness assertion described above under
  "Mathematical decisions." Fixed by adding the explicit guard to
  `kli_binomial_ratio`, not by weakening the test.
- After fix: **3777/3777 passed**, ~40s wall time, 0 failures, 0 errors.
  The delta (+181) is exactly the new `KLI φ_u (M02)` testset in
  `test/kli_phi_test.jl`. No pre-existing testset's pass count changed.

## Verification status
| Layer | Status | Oracle |
|---|---|---|
| `production_slots` reads `Event.r` correctly | Verified (new) | `test/kli_phi_test.jl`, cross-checked against `SEIR`/`MERS` `Event.r` values printed at the REPL |
| Saturation enumeration — bounds always `0 <= s_d <= min(r_d,ℓ_d)` | Verified (new) | `test/kli_phi_test.jl` "enumeration bounds" testset, 14 `(r,ℓ)` combinations, including a cardinality check (`length(S) == ∏(min(r_d,ℓ_d)+1)`) |
| Saturation enumeration — `ℓ_d=0` / `r_d=0` collapse | Verified (new) | `test/kli_phi_test.jl` "boundary ℓ=0" |
| φ_u exact values — r=(0),(1),(2),(1,1) | Verified (new), hand-computed | `test/kli_phi_test.jl`; see "Report" numbers below |
| φ_u boundary: ℓ=n (fully saturated) | Verified (new), hand-computed | `test/kli_phi_test.jl` "boundary ℓ=n" |
| φ_u boundary: n<r (denominator 0) | Verified (new), documented convention | `test/kli_phi_test.jl` "boundary n<r" |
| φ_u defensiveness on out-of-range `s` | Verified (new); one real bug found+fixed during this milestone | `test/kli_phi_test.jl` "defensive φ_u" |
| Cross-check vs `mers_filter_suite.tex` (TCC, THH, THC, TCH, RC/RH/SC/SH) | Verified (new) | see table below, exact line citations |
| Julia `binomial` boundary behavior | Verified (new), not assumed | `@assert`s at the top of `mgp_phi.jl` |
| Q_u / marginalization over m (Φ_u) | Not implemented — explicitly out of scope | N/A, deferred to M03/M04 |

### Cross-check vs `mers_filter_suite.tex`

| MERS event (real `Event` from `PhyloPOMP.MERS.events`) | tex derivation lines | tex table rows | Match |
|---|---|---|---|
| `transmission_cc` (TCC, r=(2,0)) | 472-480 | 554-556 | φ(0,0)=3/10, φ(1,0)=3/10, φ(2,0)=1/10 — exact match |
| `transmission_hh` (THH, r=(0,2)) | 482-485 | 557-559 | φ(0,0)=1/2, φ(0,1)=1/2, φ(0,2)=1/6 — exact match |
| `transmission_hc` (THC, r=(1,1)) | 487-498 | 560-563 | φ(0,0)=9/20, φ(0,1)=3/20, φ(1,0)=3/20, φ(1,1)=1/20 — exact match |
| `transmission_ch` (TCH, r=(1,1)) | 500-509 | 564-567 | same 4 values as THC (mathematically identical product for r=(1,1); confirms the generic function reproduces both independently-derived tex sections) — exact match |
| `removal_c`/`removal_h`/`sampling_c`/`sampling_h` (r=(0,0)) | 511-539 | 568-571 | φ=1 (empty product), s forced to (0,0) — exact match |
| Step B (saturation enumeration rule) | 443-449 | — | `enumerate_saturations` reproduces `s_d ∈ {0,...,min(r_d,ℓ_d)}` exactly |
| Step C (binomial ratio formula) | 450-457 | — | `kli_binomial_ratio` is a direct transcription, verified numerically above |

All cross-checks used `n=(5,4)` / `ℓ=(2,1)` (or the single-deme analogues
`n=5,ℓ=2` / `n=4,ℓ=1`) as concrete instances, computed by hand from the
tex's closed-form expressions and from the general binomial formula
independently, then asserted equal in code.

## Known failures / unresolved issues
- None outstanding. The one bug found (negative-`s` non-defensiveness) was
  fixed within this milestone, not deferred.
- **Judgment call flagged for the advisor**: the milestone spec left it as
  "your choice, document it" whether `kli_binomial_ratio` should be
  defensive against an out-of-bounds `s`. I chose to make it fully
  defensive (`s_d > r_d` and `s_d < 0` both return `0`) rather than only
  relying on callers to pass saturations from `enumerate_saturations`.
  Rationale: cheap to guarantee, and it converts a class of "silently wrong
  number" bugs (like the one I found) into "well-defined zero" instead.
  If M03 wants `kli_binomial_ratio` to be a hot inner loop, this adds two
  cheap `Int` comparisons per call — negligible, but worth noting since
  it's an extra branch not in the literal spec formula.
- **Not an ambiguity in the tex, but worth flagging**: `mers_filter_suite.tex`'s
  φ_u table (lines 552-577) displays every ratio multiplied by a support
  indicator `\Supp` (KLI Eq. 37, `𝟙_{I_d≥ℓ_d}`) that is separate from the
  binomial-ratio product itself. This milestone's `kli_binomial_ratio` does
  *not* implement `\Supp` — it's subsumed by the `n_d < r_d` defensive-zero
  behavior in some but not all cases (e.g. `\Supp` also gates the `r=(0,0)`
  rows RC/RH/SC/SH, where `φ_u=1` unconditionally per Step C, but the
  transition is only *compatible* — `Q_u` — when `I_d≥ℓ_d`; that's a
  `Q_u`/support concern, correctly out of scope for M02, not a φ_u concern).
  Flagging so M03 does not conflate "φ_u itself already encodes `\Supp`"
  with the actual situation ("φ_u is the bare binomial ratio; `\Supp`/`Q_u`
  are separate, M03 work").
- No ambiguity found in `mers_filter_suite.tex` itself that blocked this
  milestone — Steps B/C and the reference table were unambiguous and
  reproduced exactly.

## Git state
- branch: `atpabuser-devel`
- commit if created: none
- uncommitted files:
  - `docs/compiler/architecture_before.md` (from M00, still uncommitted)
  - `docs/compiler/compiler_roadmap.md` (from M00, still uncommitted)
  - `handoffs/M00_reconnaissance.md` (from M00, still uncommitted)
  - `handoffs/M01_population_ir.md` (from M01, still uncommitted)
  - `src/examples/mgp_audit.jl` (from M01, still uncommitted)
  - `test/population_ir_test.jl` (from M01, still uncommitted)
  - `handoffs/M02_kli_phi.md` (new, this file)
  - `src/examples/mgp_phi.jl` (new)
  - `test/kli_phi_test.jl` (new)
  - `src/examples/Examples.jl` (modified — one `include` line added)
  - `test/runtests.jl` (modified — one `include` line added)
  - No file was committed or pushed by this milestone, per instructions.

## Resume instructions
1. Advisor reviews `src/examples/mgp_phi.jl` and `test/kli_phi_test.jl`,
   in particular the two documented judgment calls above (the `n_d < r_d`
   defensive-zero convention and the negative-`s` guard), and confirms
   agreement before M03 builds on top of them.
2. M03 ("Full KLI Compatibility Lowering") should treat
   `production_slots`/`enumerate_saturations`/`kli_binomial_ratio` as
   building blocks to call, not reinvent: it needs to derive `Q_u(y,y')`
   (the compatibility indicator, currently an external `Q` argument here)
   from actual genealogy/coloring state, and the coloring-operator
   transitions themselves (chop χ / swap σ / fork κ, `y → y'`) that a given
   saturation `s` corresponds to — `enumerate_saturations` gives you the
   `s` values; M03 must decide, per `s`, which coloring transition(s) it
   licenses (see `mers_filter_suite.tex` Step D, lines 459-467, and the
   Complete Reference Table's `Q_u(y,y')` column, lines 552-577, as the
   worked-example oracle for that mapping).
3. Before M03 starts, resolve the `\Supp` flag above: decide whether the
   support indicator `𝟙_{I_d≥ℓ_d}` becomes part of `Q_u` (my working
   assumption) or needs a separate generic function.
4. Run `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` after any
   M03 change and confirm the count does not regress below 3777/3777
   (this milestone's new baseline).

## Next milestone
M03 — Full KLI Compatibility Lowering

## Context note
The most important thing to preserve if this conversation were compacted:
`src/examples/mgp_phi.jl`'s `production_slots`/`enumerate_saturations`/
`kli_binomial_ratio` are now the canonical, tested, generic implementations
of KLI Eq.9's production-slot binomial ratio and its saturation space —
verified via exact `Rational{Int}` arithmetic against hand-derived numbers
and against `mers_filter_suite.tex`'s TCC/THH/THC/TCH/RC/RH/SC/SH
derivations, using real `Event` objects from the unmodified `SEIR`/`MERS`
models (not hard-coded per-event saturation lists). `Q_u` is deliberately
left as an external, un-derived parameter (default 1) — that derivation,
plus the actual coloring-operator transition logic (`y → y'`), is M03's
job and should consume these three functions rather than re-deriving
saturation enumeration or the binomial ratio from scratch. One real bug was
found and fixed during this milestone's own test-writing (negative `s` not
originally zeroed by `kli_binomial_ratio`) — a concrete example of why the
milestone's "write explicit reference-value tests, not just property
tests" instruction mattered in practice.
