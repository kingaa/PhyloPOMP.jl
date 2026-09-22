# Handoff — M14: Sync with Aaron's guided-filter rewrite (September 2026)

## Status
COMPLETE (2026-09-18). Full test suite: see "Verification" below.

## Objective
Aaron King's repo (`../PhyloPOMP.jl_aarons`, HEAD `44edaf6`) rewrote
`src/examples/mers_guided.jl` (his commit `25039d0`, "guided filter for MERS
problem"), restored `mers_naive.jl` (`44edaf6`), removed `@inbounds` from
`guide.jl` (`a97b501`), and moved to the registered
PartiallyObservedMarkovProcesses (POMP) v0.10. Commit `93ebd17` on our branch
had already copied his `mers_guided.jl`, `mers_naive.jl`, `mers_tree.jl` and
the two tests in, but nothing around them was updated. This milestone brings
the rest of `atpabuser-devel` into line: the Soft/Hard MERS kernels, the
dependency set, the compiler layer, the forward simulator, and the documents.

## Relay to Aaron: two likelihood errors in the new `mers_guided.jl`

Both are fixed in our copy (`src/examples/mers_guided.jl`); Aaron's repo does
not have the fixes. Every other MERS kernel in the repo (`mers_naive.jl`, the
old `mers_funs.jl`, the compiled `mgp_mers_filter.jl`) already had the correct
terms, and the compiler IR computes both generically.

### Bug 1 — within-host fork boost is missing a factor 2 (constant offset)
`singular_branch!` charged `-log(I_c*(I_c-1))` on a camel→camel fork (and the
human analogue). The correct factor is `1/binomial(I_c,2) = 2/(I_c(I_c-1))`,
i.e. `-log(I_c*(I_c-1)/2)`, as in `mers_naive.jl` and as the IR gives
(`full_transitions(transmission_cc, ℓ, n)` → `ForkTransition.phi` =
`1//1, 1//3, 1//10, 1//28` at `I_c = 2, 3, 5, 8`). Effect: a `log 2` penalty
on every within-host coalescence in the sampled coloring. When the data force
the host assignment at every internal node (e.g. an all-camel tree with
β_hc = β_ch = 0) this is a pure offset of `(#internal nodes)·log 2`. On a
mixed tree the likelihood sums over colorings, and the bug halves the
within-host component (∝ β_cc, β_hh) while leaving the cross-host components
(∝ β_hc, β_ch) intact, so it tilts the surface toward cross-host explanations
and biases the camel/human transmission split. Milder than bug 2, but real on
the empirical MERS tree, where many internal nodes are host-ambiguous.

End-to-end check (4 camel tips, 3 internal nodes, β_hc = β_ch = 0, 40 reps of
`pfilter(Np=4000)`): Aaron's version `-27.02 ± 0.31`, patched `-24.98 ± 0.30`;
gap `2.04 ± 0.43` versus the predicted `3·log 2 = 2.08`.

### Bug 2 — regular within-host births charge no "no visible fork" factor (biases inference)
`regular_transmission_cc!`/`regular_transmission_hh!` returned a zero
log-weight. A within-deme birth in regular time (no observed node) must be
weighted by the probability that it did **not** join two of the `ℓ_d` tracked
lineages, `1 - binomial(ℓ_d,2)/binomial(I_d,2)`. Independent anchors:
`mers_naive.jl` charges exactly this on its k==1/k==2 branches; the
Chu–Vandermonde testset in `test/kli_mers_compiled_test.jl` shows the no-fork
mass is `Φ_identity + ℓ·Φ_inline = 1 - C(ℓ,2)·Φ_fork`, strictly below 1 when
`ℓ_d ≥ 2`; and Aaron's kernel proposes cc/hh births at the full population rate
with proposal probability 1 while applying neither the `×Φ_noop` jump factor
nor the equivalent `α(1-Φ_noop)` decay term, so that mass is never charged.
`seir_guided.jl` has no analogue because SEIR infection is cross-deme (I→E);
its identity/cross boosts check out against `choose_branch`'s normalization.

Numerical check (4 camel tips with long coexisting branches, β_cc = 3,
N_c = 40, I_c0 = 0.25, 30 reps of `pfilter(Np=4000)`):

| kernel | logLik |
|---|---|
| pre-rewrite GuidedMERS (`93ebd17~1`, has the factor) | −28.63 ± 0.04 |
| Aaron's rewrite, bug 1 fixed only | −22.58 ± 0.30 |
| Aaron's rewrite, both fixes (now on disk) | −28.58 ± 0.05 |

The effect scales with the number of within-host births while `ℓ_d ≥ 2`, so
it biases estimates of β_cc, β_hh and the initial conditions strongly, and it
does so even when every node's host assignment is forced by the data.

### Minor — `check()` assertion strings
`GuidedMERS.check`/`GuidedSEIR.check` interpolate `n.time` in their `@assert`
messages, but `GenealNode` has `slate`, not `time`. A failing assertion throws
a field-access error instead of the intended message. Left as-is (upstream
cosmetic); see the simulator agent's confirmation below.

## What changed on `atpabuser-devel`

### Load break fixed
`93ebd17` replaced `mers_tree.jl` with Aaron's, which defines `mers_newick`
only. `mers_soft.jl:30` and `mers_hard.jl:27` still called
`first(mers_trees)`, so the package failed to precompile
(`UndefVarError: mers_trees`). Both now parse `mers_newick`.

### Dependencies (Project.toml, test/Project.toml)
- POMP compat `0.8` → `0.10` (registered release; the `[sources]` devel pins
  are gone from both project files). POMP 0.10 adds `mif`, `traces`,
  `geometric_cooling` (used by Aaron's `test/mers_guided.jl`), drops the
  `paramsymbs` export and `flexmap.jl`, and changes `logdmeasure`'s `y` from a
  3-D array to a vector. Nothing under `src/`, `test/`, or `scripts/` used the
  removed pieces.
- `Distributions` added as a dependency (parity with upstream; the test uses
  `LogNormal`). `BenchmarkTools`, `Crayons`, `Tally` moved out of the main
  `[deps]` (test-only). `Random` kept (used by `rcateg` and `simulate.jl`).
- Version `0.0.13-3`, matching upstream. `.github/workflows/CI.yml` now tests
  Julia 1.13 as upstream does.
- `Manifest.toml` (gitignored) was regenerated from upstream's.

### Synced verbatim from upstream
`src/guide.jl` (no `@inbounds`), and the `[`@demes`](@ref `@demes`)` docstring
spelling in `genealogy.jl`, `parse.jl`, `cblv.jl`, `coloring.jl`.
`mers_guided.jl` is Aaron's file plus the two fixes above (four lines).

### Soft/Hard MERS kernels refactored to the new architecture
`mers_funs.jl` is now the MERS analogue of `seir_funs.jl`, written in the
new `mers_guided.jl` style and `include`d by `SoftMERS` and `HardMERS` only:

- Copied verbatim from the (fixed) `mers_guided.jl`: `knowledge!`, `check`,
  `deme_occupancy`, `live_condition`, `singular_root!`, `terminal_sample!`,
  `singular_sample!`, `singular_branch!`, `singular_part!`, `mers_rinit`,
  `filter_pomp`. All three MERS kernels therefore share one singular part,
  one nested `state=(;S_c,I_c,S_h,I_h)` NamedTuple with a `live` flag and
  `init_state`, one `check(gen)`, and the same defaults (β_ch = β_hc = 1,
  χ_h = 1, B = 10; the old Soft/Hard default χ_h = 0 forced `-Inf` at the
  first human tip).
- Shared rate pieces in `seir_funs.jl`'s `on*/off*` keyword style so Soft and
  Hard use one `event_rates!` (12 slots): `transmission!` (cc, hh,
  hc-identity, hc-swap, ch-identity, ch-swap; leftover `rate - Σα` returned
  as decay), `removal!` (`pi = 1-ℓ/I`, gated on `I > ℓ`, leftover
  `rate - Σα`), `demography!`, `sampling`. Soft passes `onC = ellC`; Hard
  passes `onC = sum_relhaz(rh, node, cols, Camel, Human)`, `offC = I_c-ellC`,
  exactly as `seir_soft.jl`/`seir_hard.jl` do.
- `mers_soft.jl` / `mers_hard.jl` now contain only the module preamble and
  `regular_part!(cols, state, guide, n, t, tf; kwargs...)`, taking the
  interval from `rprocess` rather than from `guide[n].tbeg/tend`, with the
  `k > 0` guard. Soft's removal bookkeeping was unified with Hard's
  (`pi = 1-ℓ/I` charged by the generic `-log(pi[k])` line); the value-level
  equivalence to the old hand-charged convention is proved case by case in a
  comment at the top of `mers_soft.jl`.
- Public signatures unchanged: `SoftMERS.filter_pomp(gen, m; kwargs...)`,
  `HardMERS.filter_pomp(gen, m; kwargs...)`, `*.mers_tree`, `*.Demes`.

Regression against the pre-refactor kernels on the 3-tip `small_tree`
fixture (`B_c = B_h = 0` so the runs are like-for-like; Np = 1000, 10 reps):

| kernel | before | after (three seeds) |
|---|---|---|
| Soft | −111.88 ± 0.66 | −109.40 ± 0.88, −109.82 ± 1.06, −112.83 ± 0.49 |
| Hard | −106.64 ± 2.22 | −109.74 ± 1.36, −110.68 ± 0.63, −110.72 ± 1.59 |

The `logmeanexp` estimator is heavy-tailed on this fixture; all values sit
inside its spread. `test/mers_soft.jl` replaces the old "soft vs guided
differ under a shared seed" check (true of any two kernels, so vacuous) with
an estimate-level agreement test: Soft, Hard and Guided `logmeanexp` over 15
reps at Np = 2000 must agree pairwise within 3 combined SE (measured
z = 1.00, 0.58, 1.12), plus a finite-logLik check at the new defaults.

Independent higher-precision check (same fixture and guide, `B = 0`,
Np = 10000, 20 reps each, separate seeds), including the pre-rewrite guided
kernel from `93ebd17~1` as a fourth estimator:

| kernel | logLik |
|---|---|
| Soft (refactored) | −106.86 ± 1.19 |
| Hard (refactored) | −105.78 ± 2.85 |
| Guided (Aaron + two fixes) | −106.12 ± 0.65 |
| pre-rewrite Guided | −106.92 ± 0.62 |

All six pairwise z-scores are below 0.9. The four kernels estimate the same
likelihood; Guided has the smallest variance.


### Forward simulator (`src/simulate.jl`)
No code changes were needed; POMP 0.10's `simulate.jl` is unchanged from 0.8,
our `simulate(model::MGPModel, θ; ...)` method is disjoint from POMP's
`simulate(object; nsim, ...)`, and `Test.detect_ambiguities(PhyloPOMP)` is
empty. New discriminating tests feed simulated genealogies straight into the
guided filters (the filters that now carry the trusted math):

- `test/mers_simulate.jl`: the headline check validates through both
  `SoftMERS.filter_pomp` and `GuidedMERS.filter_pomp` (same realization,
  re-simulated with `demeset=GuidedMERS.Demes`,
  `samplemap=[GuidedMERS.Camel, GuidedMERS.Human]`, asserted identical via
  `newick`). A new testset asserts `GuidedMERS.check(g) === nothing`, the four
  degree/deme invariants, and a finite `logLik` at `Np=2000` with guide
  `fsmarkov(Camel=>0.5, Human=>0.5, (Camel,Human)=>0.05)`. 36 tests (was 20).
- `test/seir_simulate.jl`: new testset asserts `GuidedSEIR.check`, a finite
  `logLik` at `Np=1000`, and that the simulated tree contains an inline
  (non-destructive) sample, exercising `GuidedSEIR.inline_sample!`. 30 tests
  (was 21).

Review findings (read-only, cited by the agent against `src/simulate.jl`):

- Structure always satisfies the guided `check()`s: a BIRTH creates one node
  and opens both lineages on it (`simulate.jl:176-183`), so no node exceeds
  two children; `prune!` removes childless and collapses degree-1 `Node`s;
  the Root never enters the inventory. MERS samples are terminal because the
  `@mgp` sampling events decrement `I_c`/`I_h`. Verified on 60 genealogies
  per model.
- Deme metadata is exactly what `knowledge!`/`terminal_sample!` consume: only
  `Sample` nodes get `samplemap[d]`; a mismatched `samplemap` module fails
  loudly with a `MethodError`, never silently.
- Edge case: an extinct-unsampled run returns an **empty** `Genealogy`;
  `check` passes vacuously but `filter_pomp` throws `BoundsError` at
  `guide.jl:303` (`timezero(g::Guide) = g.nodes[1].tbeg`). Callers need a
  retry loop (both tests have one).
- Naming traps between the `@mgp` tables and the filters: SEIR's population
  is `N` in the table but `pop` in the filter; `@mgp SEIR` declares `χ` but has
  no destructive-sampling event, so the SEIR filters must be given `χ=0` for
  the two encodings to agree. MERS hazards match one-for-one (including the
  `B_c*S_c/N_c` death hazard). `GuidedMERS` has no non-destructive sampling
  path, so a future `sample(I_c)` event would be rejected by its `check`.
- `check()` assertion strings: confirmed on Julia 1.13 that a Sample/Node
  degree violation throws `FieldError: type GenealNode has no field 'time'`
  rather than the intended message (Root branches use `slate` and render
  correctly). Also `n.deme in [Human, Camel]` with a genuinely `missing` deme
  throws `TypeError: non-boolean (Missing)` instead of "unknown deme".
- Wart: `include("mgp_mers.jl")` appears in both `mgp_macro.jl` and `mgp.jl`,
  so `MERS` and `MERSDemes` are defined twice per load (identical, harmless).
  Not changed here.


### Tests
- `test/mers_guided.jl` (Aaron's): the ~170-iteration `mif` block on the
  548-node tree is now inside `if heavy ... end` (`RUN_HEAVY_TESTS`, default
  on) and asserts on its result; the `Np=1000` `pfilter` smoke check always
  runs. Measured once: the full file takes 4 min 41 s (the 100-iteration
  `mif` alone 211 s); with `RUN_HEAVY_TESTS=no` it takes 9 s.


## Compiler layer (`mgp_*.jl`): verdict
No rewrite needed now. Both compiled filters (`mgp_seir_filter.jl`,
`mgp_mers_filter.jl`) anchor on `NaiveSEIR.singular_part!` /
`NaiveMERS.singular_part!`, whose signatures did not change, and both Gate-5
bit-exact tests pass under POMP 0.10 (SEIR: 238/238 comparisons, worst
|Δll| = 1.4e-14; MERS: 718/718, worst 7.1e-15).

The forward-looking point is the more useful one. Aaron's new file structure
is already the compiler's structure:

| Aaron's `mers_guided.jl` piece | `@mgp` / IR counterpart |
|---|---|
| `transmission!` (4 hazards, zero decay) | the four `BIRTH` events' `hazard` closures |
| `removal!` (reduced proposal `rate·(1-ℓ/I)`, leftover `rate - Σα`) | `DEATH` events + `compiled_decay`'s leftover term (M08) |
| `demography!` | `NEUTRAL` events |
| `sampling` (adds to decay) | `SAMPLE` events in `total_decay` |
| `singular_root!` via `deme_occupancy` | `n(x)` over `model.demes` |
| `singular_branch!` fork weights | `ForkTransition.phi` (`1/(n_i n_j)` or `2/(n(n-1))`) |
| `regular_transmission_cc!` no-fork factor | `Φ_identity + ℓ·Φ_inline` (Chu–Vandermonde aggregate, M09 Finding 2) |
| `event_rates!` composition | `regular_step!` loop in `mgp_filter.jl` |

The two upstream bugs are exactly the two entries the IR computes generically
(bug 1 = `ForkTransition.phi`; bug 2 = the no-fork aggregate). A compiled
guided MERS filter driven from the `@mgp MERS` table would have produced both
terms correctly by construction. That is the concrete argument for finishing
the generic `kli_select`/`apply_move!`/`singular_update!` slots in
`mgp_filter.jl` with Aaron's modular signatures as the target shape: each of
his rate pieces is one `EventType` group, and each of his singular functions is
one `NodeType` dispatch. Not done in this milestone.

## Verification
- `julia --project -e 'using Pkg; Pkg.test()'` with `RUN_HEAVY_TESTS=no`:
  **passed** (all 30 files included by `test/runtests.jl`, POMP 0.10.0, Julia 1.13.0).
- Heavy blocks run individually: `test/mers_guided.jl` with the `mif` block
  (4 min 41 s, passes); `test/mers_soft.jl`/`test/mers_hard.jl` heavy
  benchmarks (pass; see the refactor agent's numbers above).
- Compiled-filter Gate 5 under POMP 0.10: SEIR 238/238 bit-exact, MERS
  718/718 bit-exact.
- Kernel agreement at Np = 10000: table above.
- Bug 1 / bug 2 numerical checks: tables above (scratch scripts, not kept).


## Files touched
Modified: `.github/workflows/CI.yml`, `Project.toml`, `test/Project.toml`,
`README.md` (upstream header), `src/guide.jl`, `src/genealogy.jl`,
`src/parse.jl`, `src/cblv.jl`, `src/coloring.jl` (upstream sync),
`src/examples/mers_guided.jl` (two fixes), `src/examples/mers_funs.jl`,
`src/examples/mers_soft.jl`, `src/examples/mers_hard.jl` (refactor),
`src/examples/mers_filter_suite.tex` (changelog), `test/mers_guided.jl`
(heavy gating), `test/mers_soft.jl` (agreement test), `test/mers_simulate.jl`,
`test/seir_simulate.jl` (guided-filter validation).
Untracked additions: `handoffs/M14_upstream_guided_sync.md` (this file),
`docs/compiler/guided_kernel_assessment.{tex,pdf}` (addendum added; the file
was already untracked), `mers_guided_part1.md` (note prepended).
Not touched: `src/examples/Examples.jl`, `test/runtests.jl`, all `mgp_*.jl`,
`mers_naive.jl`, `seir_*.jl`, historical milestone notes.
Nothing was committed.

## Open items
1. Relay bug 1, bug 2 and the `check()` message issue to Aaron (upstream).
2. Decide whether `simulate` should throw a clear error on an empty
   genealogy rather than letting `filter_pomp` fail at `guide.jl:303`.
3. Remove the duplicate `include("mgp_mers.jl")` (`mgp.jl` vs `mgp_macro.jl`).
4. Optional: compile a guided MERS kernel from the `@mgp MERS` table using
   Aaron's modular signatures as the target (see the compiler verdict).

