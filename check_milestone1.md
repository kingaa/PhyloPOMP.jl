# Review brief: Milestone 1 — generic forward genealogy simulator

**Audience:** a reviewing agent that has *not* seen the implementation
session. This document is self-contained; it does not assume access to the
conversation that produced it.

**Branch:** `atpabuser-devel`. **Status:** implemented, tested, full suite
green. Not yet committed at the time of writing.

## 1. Task and confirmed gap

PhyloPOMP.jl (the Julia phylodynamic-inference package) had no forward
genealogy simulator — only particle-filter kernels that *consume* an
already-known genealogy (`src/examples/{seir,mers}_{naive,soft,guided,hard}.jl`).
Confirmed before writing any code:
- `TODO.md:5` listed "genealogical simulations" as unimplemented.
- All example trees (`src/examples/seir_trees.jl`, `mers_tree.jl`) are
  static Newick string literals baked into source, not simulator output.
- No YAML dependency, no model-config mechanism anywhere in the repo.

R's `phylopomp` (kingaa/phylopomp, the sibling R package) has such a
simulator: a Doob-Gillespie SSA kernel (`popul_proc_t`/`master_t`/
`genealogy_t`/`ball_t`/`inventory_t` in its `src/*.h`) that is generic C++,
parameterized per-model via YAML → compile-time C++ codegen
(`yaml/add_model.R`). The decision taken here (confirmed with the user) was
to build a **Julia-native** equivalent rather than ingest R's YAML files
directly (R's YAML `rate:`/`jump:` fields contain literal C++ snippets;
ingesting them would mean writing a C++-subset-to-Julia translator, out of
scope). Julia doesn't need R's per-model compile step at all: rate/jump
closures are already ordinary JIT-specialized Julia functions, so one
generic engine can drive any declaratively-specified model.

This repo already had the right declarative vocabulary for that: `@mgp`/
`Event`/`MGPModel` in `src/examples/mgp.jl` — a Julia macro DSL producing an
auditable event table (name, population stoichiometry `Δ`, hazard closure,
one of five `EventType`s, source/target deme indices), previously used only
by an unfinished particle-filter scaffold (`src/examples/mgp_filter.jl`,
whose own header says its weight-bearing functions are "UNVERIFIED STUBS").
Milestone 1 adds a **second, independent consumer** of `MGPModel` — a
forward simulator — and does not touch or depend on `mgp_filter.jl` at all.

## 2. Files touched

| File | Change |
|---|---|
| `Project.toml` | Added `Random` as an explicit dependency (stdlib; needed for `AbstractRNG` threading). |
| `src/rcateg.jl` | Added an optional `rng::AbstractRNG` keyword (default `Random.default_rng()`) to all four `rcateg` methods. Purely additive — every pre-existing positional call site is unaffected; verified via the full test suite (see §5). |
| `src/simulate.jl` | **New.** The generic simulator engine (detailed in §3). |
| `src/PhyloPOMP.jl` | Added `include("simulate.jl")`, placed *after* `include("examples/Examples.jl")` since the engine depends on `MGPModel`/`Event`/`SEIR`, which are defined there. |
| `test/seir_simulate.jl` | **New.** The acceptance test (detailed in §4). |
| `test/runtests.jl` | Added `include("seir_simulate.jl")` between `seir_macro_equivalence.jl` and `seir_naive.jl`. |

No existing file's *behavior* changed. `src/examples/mgp_filter.jl`,
`src/coloring.jl`, and every hand-coded filter kernel are untouched.

## 3. Design of `src/simulate.jl`

### 3.1 Why not reuse `Coloring` (`src/coloring.jl`)?

`Coloring`'s `BitSet`s hold *filter-tracked lineage indices* (size ~ℓ, the
branches a particle filter colors against an already-fixed tree). Forward
simulation needs one live entry per **currently-extant individual in the
simulated population** (size ~ population count), a different growth
regime — a `BitSet` allocates up to its max element, which would grow
unboundedly over a simulation's event count if lineage IDs were minted the
way `Coloring` expects. A new type, `SimInventory` (`Vector{Vector{Name}}`,
one vector per deme, holding genealogy-node *names*), is used instead. It's
deliberately pattern-compatible with `Coloring`'s vocabulary (a migration
relabels deme membership with no new node, exactly matching what `swap!`
does) without inheriting its filter-specific invariants.

### 3.2 The core loop (`simulate`, bottom of the file)

Standard Doob-Gillespie direct method:
1. Compute `hazard(x,θ)` for every event in `model.events`.
2. Total rate `≤ 0` → absorbing state, stop.
3. Draw `Δt ~ Exp(total)`; if it would cross `tmax`, stop (censor).
4. Otherwise pick an event via `rcateg` (weighted by hazards, now
   `rng`-threaded), apply its `Δ` to the population state (`apply_delta`),
   and apply its genealogy effect (`apply_event!`).
5. **Self-consistency assertion** after every event: for every deme `d`,
   `length(inventory[d]) == getproperty(x, model.demes[d])`. This holds by
   construction (forward simulation tracks *every* extant individual,
   unlike a filter which tracks a sparse subset) and is cheap to check —
   included specifically so an engine bug fails loudly and immediately
   rather than producing a silently-wrong tree.

### 3.3 `apply_event!` — where genealogy topology is built

This is the part most worth scrutinizing; it's the one place the SSA's
event stream turns into tree structure. Per `EventType`:

- **BIRTH**: pick a uniformly random live lineage in deme `ev.from`
  (`rand(rng, inv[d])`). Its current node gets **one child continuing the
  same lineage** (same deme) plus **one new child per deme in `ev.into`**.
  E.g. SEIR's `infection` event (`from=I, into=[E]`) bifurcates an
  infectious lineage into "stays infectious" + "newly exposed" — the
  infector doesn't vanish, so this must NOT be modeled as the parent lineage
  moving to the child deme; it must gain a sibling. Verified this is what
  the code does (`c1` continues in `d`, `cj` for each `j ∈ ev.into` is new).
- **MIGRATION**: relabels the chosen lineage's deme membership in the
  inventory; **no new node** — a migration is not itself an observable
  genealogy feature in this model class (matches `Coloring`'s `swap!`).
- **DEATH**: the lineage is removed from the inventory; its current node is
  left as an unsampled dead-end tip (0 children, `type==Node`), to be
  removed by `prune!`.
- **SAMPLE**: the current node is marked `type = Sample`, and (if
  `samplemap` given — Milestone 2 concern, `nothing` for SEIR) its `deme`
  field is set. Whether the lineage continues is read off `ev.Δ`: **if `Δ`
  decrements the sampled deme's own compartment, the sample was destructive**
  (lineage removed, matches `mers_naive.jl`'s `@assert length(n.children)==0`
  for MERS's `sample_remove`-only sampling); otherwise it's non-destructive
  and gets a single continuation child (a serially-sampled pass-through
  node, matching `NaiveSEIR.singular_part!`'s `n.type==Sample` handling of a
  node with exactly one child). SEIR's `sampling` event has `Δ=()`, so this
  path always takes the non-destructive branch in Milestone 1's test.
- **NEUTRAL**: no lineage involved at all (`ev.from==0`, e.g. SEIR's
  `waning`); only `Δ` (already applied by the caller) has any effect.

**A subtlety worth checking carefully:** the destructive/non-destructive
distinction is derived generically from `Δ`, not hardcoded per model. This
was necessary because `Event` itself doesn't otherwise distinguish
`move=sample(...)` (non-destructive) from `move=sample_remove(...)`
(destructive) at the type level — both produce `EventType.SAMPLE`
(`src/examples/mgp_macro.jl`'s `_decode_move`); only their `Δ` differs in
practice. Reviewer: check this inference is actually correct for every
model's sampling event, not just SEIR's.

### 3.4 `prune!` — the gap identified in `repair!`

`repair!` (`src/genealogy.jl:108`) only drops nodes that are *both*
parentless and childless (an isolated dead root). It does **not** prune
unsampled-extant tips (lineages alive at `tmax`, never sampled) or extinct
dead-ends, and does not collapse resulting degree-1 internal nodes. Without
this, a freshly-simulated tree would contain dangling, non-`Sample` tips a
filter never expects.

`prune!` runs as two fixed-point passes over a `Dict{Name,GenealNode}`
working copy (chosen over in-place `Vector` mutation to avoid index-shift
bugs while deleting):
1. Repeatedly remove `type==Node` tips with 0 children (never touches
   `Root` or `Sample` types) until none remain.
2. Repeatedly collapse `type==Node` nodes left with exactly 1 child
   (reattach grandchild directly to grandparent — correct without summing
   branch lengths manually, since `slate` is an absolute time and the
   grandchild's own `slate` already reflects it).
3. A `Root` node that ends up childless (its entire founding lineage was
   never sampled) is deliberately left in place — `repair!`, called right
   after `prune!` in `simulate`, already drops parentless-and-childless
   nodes, so this composes correctly without duplicating that logic.

**Reviewer: this is hand-rolled tree surgery with no independent
implementation to compare against** (R's equivalent, `genealogy.h`'s
prune/obscure, wasn't ported, only its *purpose* was matched). The
strongest evidence it's correct is indirect (§4's headline test would very
likely fail if it produced a structurally malformed tree — a filter
computing rates over children/parent pointers is sensitive to exactly this
kind of error), but a direct code read of `prune!` is worth the time.

## 4. The acceptance test (`test/seir_simulate.jl`)

The headline check, and the reason to trust this more than a round-trip
test: **an independently-written, already-tested filter accepts a tree the
new engine built from scratch, at the *same* parameters, and assigns it
finite likelihood.**

```julia
θ = (β=4.0, σ=1.0, γ=1.0, ω=1.0, ψ=0.10, χ=0.0, N=100.0)
x0 = (S=99, E=0, I=1, R=0)
g = simulate(PhyloPOMP.SEIR, θ; x0=x0, graft=[0,1], tmax=20.0, rng=...)
# ... (retry until nsample(g) ≥ 5, deterministic given the seed)
p = NaiveSEIR.filter_pomp(g; β=4.0, σ=1.0, γ=1.0, ω=1.0, ψ=0.10, χ=0.0,
                           pop=100, S0=0.99, E0=0.0, I0=0.01, R0=0.0)
pf = pfilter(p, Np=1000)
@test isfinite(logLik(pf))
```

Why this discriminates and a Newick/CBLV round trip alone would not: a
round trip only proves `newick`/`parse_newick` (or `cblv`/`parse_cblv`) are
mutual inverses on whatever tree comes out — it would faithfully serialize
a *structurally wrong* tree just as well as a correct one. The filter
check instead requires actual topological/timing correctness: `NaiveSEIR`'s
`singular_part!` has hard assertions on children counts by node type
(`Root` exactly 1, internal `Node` exactly 2, `Sample` 0 or 1) that a
malformed tree would trip, and the guided-coloring math would not converge
to a finite likelihood if the branching structure didn't reflect a
genuine, compatible realization of the same SEIR rate structure the filter
assumes.

Full test list (21 assertions): argument validation (`x0`/`graft`
mismatch, wrong-length `graft` → `ArgumentError`), non-degenerate
simulation (single root, tip/sample-type consistency, time bounds),
Newick round trip, CBLV round trip, then the headline check above.

## 5. Validation performed (exact commands and results)

```
julia --project=. -e 'using PhyloPOMP; println("OK")'
  → OK

julia --project=. test/seir_simulate.jl (standalone driver)
  → 21/21 passed

RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl
  → 3281/3281 passed (full suite, all pre-existing suites unaffected)

julia --project=. test/runtests.jl  (heavy benchmarks on, default)
  → 3281/3281 passed, ~3 min
```

## 6. Explicitly out of scope / known limitations

- **Piped continuation** (R's `simulate(x, time=...)` resuming a
  serialized SSA clock) — not implemented. Noted in `simulate`'s docstring.
  Would need the sim state to carry the next-scheduled-but-unfired
  `(next, event)` if added later.
- **R-YAML ingestion** — deliberately deferred; see §1.
- **Performance** — `prune!`'s fixed-point `Dict` scans are O(n²)
  worst-case; fine at test scale (tens of nodes), not benchmarked at scale.
- **Multi-lineage grafts** (`graft` with an entry `>1`, producing a
  multi-root forest) are supported by the code and exercised nowhere in
  the test suite beyond argument-validation checks — only the single-root
  (`graft=[0,1]`) case is exercised end-to-end. The codebase elsewhere
  clearly anticipates forests (`output/seir_forest/`,
  `scripts/seir_cblv_forests.jl`), so this isn't speculative generality,
  but it is untested; the `newick(g)`/`parse_newick` multi-string path in
  particular is unverified for simulator output.

## 7. Suggested reviewer focus (highest to lowest value)

1. `apply_event!`'s BIRTH/SAMPLE branches (§3.3) — read against the docs
   in the code and confirm the deme-assignment logic for MERS specifically
   (see `check_milestone2.md`).
2. `prune!` (§3.4) — hand-trace a small example with a birth whose
   off-branch subtree is entirely unsampled, confirm the collapsed result
   is a valid tree (parent/children pointers consistent, no orphaned
   entries).
3. The self-consistency assertion in the main loop (§3.2) — confirm it
   would actually fire for a plausible class of bugs (e.g., temporarily
   comment out one `add!`/`remove!` call in `apply_event!` and confirm the
   test suite fails loudly rather than silently).
4. Multi-root/forest path (§6) — if time permits, construct a `graft`
   vector with an entry `>1` and confirm `newick`/`parse_newick` round-trip
   through the `Vector{String}` path.
