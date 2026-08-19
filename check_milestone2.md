# Review brief: Milestone 2 — generalizing the simulator to MERS

**Audience:** a reviewing agent that has *not* seen the implementation
session. Reads standalone, but assumes `check_milestone1.md` (Milestone 1's
review brief) has been read first — this document only covers the delta.

**Branch:** `atpabuser-devel`. **Status:** implemented, tested, full suite
green. Not yet committed at the time of writing.

## 1. Goal and what turned out to already exist

The plan for Milestone 2 was "extend `@mgp` so a new model can be declared
without hand-building `Event` structs, then repeat Milestone 1's
acceptance pattern against MERS." On inspection, **the declarative model
spec already existed**: `src/examples/mgp_mers.jl` already had a complete
`@mgp MERS begin ... end` block (12 events: 4 transmission routes, 2
removals, 2 destructive samples, 4 demography events for the susceptible
pools), exported as `PhyloPOMP.MERS`. So Milestone 2 required **no new
front-end/DSL code at all** — only:
1. A generalization to the simulator engine itself (§2), needed because
   MERS's genealogy convention differs from SEIR's in one specific way.
2. A different validating filter than originally planned (§3) — a
   pre-existing API constraint made `NaiveMERS.filter_pomp` unusable for
   this purpose.
3. The acceptance test itself (§4), including a real correctness lesson
   learned by hitting it (§5).

## 2. The one real design gap: deme is observed for MERS, latent for SEIR

Milestone 1 hardcoded `deme = missing` on every simulated node, matching
SEIR's convention (`Genealogy{Unstructured}`, compartment membership never
recorded on the tree). Confirmed by reading `mers_naive.jl`'s
`singular_part!`: its `Sample` branch checks `n.deme == Camel` / `n.deme ==
Human` directly, and the genealogy it filters (`mers_tree`) is parsed with
`demes=Demes` (a real two-instance deme enum), not `Unstructured`. This
makes sense epidemiologically: which host species a sample came from is
directly observable, even though which species a lineage passed through
earlier (internal/root nodes) is not.

**Fix**, entirely additive to `src/simulate.jl`, defaults preserve
Milestone-1 behavior exactly:
- `simulate(...; demeset::Module = Unstructured, ...)` — sets the
  *returned* `Genealogy`'s type parameter. MERS callers pass
  `demeset = SoftMERS.Demes` (the exact module instance the target filter
  uses — Julia dispatches on `Genealogy{D}` by object identity, so a
  locally-invented demeset with the same deme *names* would NOT be
  interchangeable with `SoftMERS.Demes`; the caller must supply the real
  one).
- `simulate(...; samplemap::Union{Nothing,AbstractVector} = nothing, ...)`
  — if given, `samplemap[d]` is the concrete enum instance recorded as
  `deme` on any node marked `Sample` when the fired event's MGPModel-deme
  index is `d`. `nothing` (default) leaves every node's deme `missing`,
  i.e. Milestone 1's SEIR behavior is `demeset=Unstructured,
  samplemap=nothing` and is bit-for-bit unchanged (reverified — see §4 of
  `check_milestone1.md`'s validation, re-run after this change with the
  same 21/21 result).
- `apply_event!` gained a `samplemap` parameter, used in exactly one place:
  the `SAMPLE` branch does `isnothing(samplemap) || (cur.deme =
  samplemap[d])` right after marking `cur.type = Sample`. Root/internal
  nodes are never touched — only sampling observes species identity,
  matching the hand-coded filters' convention exactly.

**No other part of the engine changed.** The Gillespie loop, `apply_delta`,
`prune!`, and every other branch of `apply_event!` are untouched and
model-agnostic — MERS's 12 events (4 `BIRTH`/transmission, 2 `DEATH`/
removal, 2 `SAMPLE`/destructive, 4 `NEUTRAL`/demography) all ran correctly
through the *existing* Milestone-1 code paths on first attempt (see §5) —
this is the strongest evidence for the plan's central claim that one
generic engine, not per-model code, is enough in Julia.

**Reviewer:** confirm `model.demes[event.from]` (used to look up `sym` for
the destructive/non-destructive Δ check, see `check_milestone1.md` §3.3)
still resolves correctly for all 4 of MERS's transmission events, whose
`from`/`into` are less symmetric than SEIR's single-deme-pair case (e.g.
`transmission_hc`: `from=I_c, into=[I_h]` — a camel infects a human, the
*camel* lineage continues, a *new* human lineage is created; verify this
is what `@mgp`'s `_decode_move` actually encodes by reading
`src/examples/mgp_mers.jl` against `src/examples/mgp_macro.jl`'s `fork`
case, not just trusting the test result).

## 3. Filter substitution: `SoftMERS`, not `NaiveMERS`

The plan said "re-run the acceptance pattern against `NaiveMERS.filter_pomp`."
This is not possible as written: `NaiveMERS.filter_pomp`
(`src/examples/mers_naive.jl:252`) takes **no genealogy argument at all** —
it always filters the fixed, static, empirical `mers_tree`
(`filter_pomp(; Beta_cc=4.0, ...)`, no `gen` parameter). This is confirmed,
pre-existing, upstream-synced behavior (see `handoff.md`'s "Follow-up round
4" section), not a bug to fix here.

`SoftMERS.filter_pomp(gen::Genealogy, m::FSMarkovProc; ...)`
(`src/examples/mers_funs.jl:369`, shared by the Soft/Guided/Hard MERS
kernels) **does** accept an arbitrary genealogy, exactly like
`NaiveSEIR.filter_pomp` did for Milestone 1. It additionally requires an
auxiliary guide process `m` (built via `fsmarkov(...)`, e.g.
`fsmarkov(Camel=>0.3, Human=>0.7, (Camel,Human)=>1)`) — this is
particle-filter machinery unrelated to the simulator; the test constructs
one following the exact pattern already used in `test/mers_soft.jl:28`.

**Reviewer:** confirm this substitution is a legitimate like-for-like swap
and not quietly testing something weaker. `SoftMERS` shares its
`singular_part!`/rate structure with `NaiveMERS` via `mers_funs.jl` (per
`mers_filter_suite.tex`'s "Implementation Map" section, cited in
`handoff.md`), and both encode the same MERS rate equations independently
of `@mgp MERS` — so the same "two independent encodings agree" argument
from Milestone 1 still holds. It is a slightly weaker independence claim
than SEIR's (Soft's proposal law differs from Naive's, but both consume
the same `β_cc/γ_c/χ_c/...` params against the same target rate structure);
worth a second opinion on whether that's an acceptable substitution or
whether `NaiveMERS.filter_pomp` should instead be extended to accept a
`gen` argument in a future round (out of scope here — see
`check_milestone1.md` §6's "untouched by design" principle, which was
extended to cover all hand-coded filter kernels, not just SEIR's).

## 4. The acceptance test (`test/mers_simulate.jl`)

```julia
θ = (β_cc=3.0, β_ch=0.5, β_hc=0.5, β_hh=3.0, γ_c=1.0, γ_h=1.0,
     χ_c=0.3, χ_h=0.3, B_c=0.0, B_h=0.0, N_c=20.0, N_h=20.0)
x0 = (S_c=19, I_c=1, S_h=20, I_h=0)
g = simulate(PhyloPOMP.MERS, θ; x0=x0, graft=[1,0], tmax=10.0, rng=...,
             demeset=SoftMERS.Demes, samplemap=[Camel,Human])
# retry (deterministic given the seed) until 4 ≤ nsample(g) ≤ 12 AND both
# species are represented among the samples
p = SoftMERS.filter_pomp(g, m; β_cc=3.0, ..., N_c=20, N_h=20)
seed!(20260814); pf = pfilter(p, Np=5000)
@test isfinite(logLik(pf))
```

Structural checks unique to MERS (beyond the SEIR-analogous set): every
`Sample` node has **exactly zero children** (MERS's sampling is always
`sample_remove`/destructive — no serial-sampling continuation path is
exercised here, unlike SEIR), every `Sample` node has a **non-missing**
deme, and every **non-**`Sample` node's deme is missing — directly checks
the §2 generalization did what it claims.

Full test list (20 assertions): argument validation (including a new
`samplemap` length-mismatch check), non-degenerate cross-species
simulation, the Sample/deme structural checks above, Newick round trip
(demes preserved through the trip, checked via `Set` comparison of sampled
demes before/after), CBLV round trip, then the headline filter check.

## 5. A real correctness-adjacent lesson: tree-size tuning, not an engine bug

First attempt at this test used `N_c=N_h=100`, `χ_c=χ_h=0.5`, and accepted
any tree with `nsample(g) ≥ 4` and both species present. Running under the
full suite (not standalone) intermittently failed the headline check with
`logLik = -Inf`. Diagnosis (not a suite-ordering/RNG-seeding bug in the
simulator — `simulate` is called with an explicit local `rng`, fully suite-
order-independent by construction):

- The accepted tree had **61 samples** — the "≥4" lower bound let a large
  stochastic outbreak through as soon as it also happened to clear it,
  with no upper bound to stop it. `pfilter(Np=5000)` collapses to `-Inf` on
  such a tree a non-negligible fraction of the time (~20%, measured
  directly: 4/20 fixed global seeds gave `-Inf` on the accepted 61-sample
  tree) — an `Np`-vs-tree-size adequacy issue in the *filter*, not evidence
  of a wrong tree.
- Separately, the test wasn't seeding the filter's RNG at all (`pfilter`
  draws from the global RNG, unlike `simulate`, which takes an explicit
  `rng`), so the outcome depended on however much global RNG state prior
  test files had consumed — reproducible standalone, flaky inside the full
  suite.

**Fix** (both applied): tuned `β`/`χ`/`N` down and added an upper bound
(`nsample(g) ≤ 12`) to the retry condition, verified via direct sampling
(500 independent realizations) that ~54% land in the target `[4,12]` ×
both-species window — comfortably within the 300-attempt retry budget;
added an explicit `seed!(20260814)` immediately before `pfilter`, matching
the established pattern in `test/mers_soft.jl:26` and elsewhere in this
suite (every other filter test in the repo seeds explicitly; this one
originally didn't). Reverified: full suite green 3× in a row, both with
and without `RUN_HEAVY_TESTS`.

**Reviewer:** this is worth an explicit second look — confirm the
diagnosis (tree-size-dependent filter collapse rate, not a hidden simulator
bug) rather than taking it on faith. The diagnostic approach used: fix the
simulated tree, sweep 20 different global `seed!` values through
`pfilter`, and observe the failure rate directly. Re-running that sweep
independently (or with a different tree/seed) would be the cheapest way to
confirm this wasn't papering over something else.

## 6. Validation performed (exact commands and results)

```
julia --project=. test/mers_simulate.jl (standalone driver)
  → 20/20 passed

RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl   (run 3×)
  → 3301/3301 passed, all 3 runs identical

julia --project=. test/runtests.jl  (heavy benchmarks on, default)
  → 3301/3301 passed, ~3 min
```

(3301 = the 3281 from Milestone 1's validation + 20 new MERS assertions.)

## 7. Out of scope / not attempted

- Extending `NaiveMERS.filter_pomp` to accept an arbitrary `gen` (§3) —
  flagged as a possible future improvement, not attempted here since it
  would modify hand-coded, already-tested production filter code, which
  both milestones' plans deliberately avoided touching.
- Milestone 3 (distributional cross-validation against local R
  `phylopomp-fork`) — not started.
