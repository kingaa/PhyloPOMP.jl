# Review brief: Milestone 3 — cross-language validation against R phylopomp

**Audience:** a reviewing agent that has *not* seen the implementation
session. Assumes `check_milestone1.md`/`check_milestone2.md` read first.

**Branch:** `atpabuser-devel`. **Status:** implemented, and — unlike
Milestones 1–2 — this one found and fixed a real bug in `src/simulate.jl`.
Not yet committed at the time of writing.

## 1. Goal

Milestones 1–2 validated the simulator against Julia's own hand-coded
filters (`NaiveSEIR`/`SoftMERS`). Milestone 3 is the "the filter agrees
with the simulator" check's independent counterpart: does the Julia
simulator's *output distribution* match an entirely separate
implementation of the same model — R's `phylopomp::runSEIR`? This is the
"cross-language gold-standard validation (never implemented)" item flagged
in `handoff.md`'s "Remaining work" section.

## 2. A real bug was found: sample events were timestamped wrong

**Symptom:** at matched parameters, `nsample` distributions and extinction
fractions matched R closely (KS p=0.52, not significant) from the start,
but first-sample-time was **systematically and significantly earlier** in
Julia (mean 1.41 vs R's 1.95 at N=400; KS D=0.24, p≈0 — unambiguous, not
sampling noise).

**Diagnosis path** (each step is independently reproducible):
1. Ruled out measurement error: parsed one of R's own saved trees through
   Julia's `parse_newick`, computed first-sample-time from it, and
   compared against R's own `lineages(x)` event log for the *same*
   realization (`event_type==1` = sample, per `src/node.h`'s documented
   coding) — exact match (`3.1613...` both ways). The extraction
   methodology was not the problem.
2. Ruled out the exponential-draw/hazard machinery: directly verified the
   mean waiting time to the first event of *any* kind, over 200,000 draws
   at `x0`, matched the closed-form `1/total_hazard(x0)` to 4 significant
   figures.
3. Found the actual bug by reading `src/simulate.jl`'s `apply_event!`
   SAMPLE branch: it did `cur.type = Sample` — **relabeling the existing
   node in place without updating its timestamp**. `cur.slate` was
   whatever time the *previous* event on that lineage fired, not the
   current sample event's actual time `t`. Every sample was recorded at
   the time of the lineage's prior event, not the sample itself —
   systematically too early, and non-obviously so, since it doesn't affect
   `nsample` (count is unaffected by which node gets which timestamp).
4. Confirmed by reading R's `sample()` (`src/master.h` in the R source):
   it calls `make_node()` to create a **fresh node at the current time**,
   rather than reusing the lineage's existing holder — exactly the
   behavior Julia was missing.

## 3. The fix, and a second bug it exposed

Naive fix (create a new, correctly-timestamped `Sample` node, then — if
non-destructive — a separate continuation child at the same instant)
**broke immediately**: `parse_newick` throws `AssertionError: dropping
zero-length branch collapses multiple samples` (`src/parse.jl`'s
`clip_zlb!`). This is pre-existing, load-bearing code — worth understanding
why before working around it, not patching over it.

Empirically isolated (see `check_milestone3.md`-adjacent scratch tests,
not committed — reproducible via the two hand-built `Genealogy` examples
described below) exactly what `clip_zlb!` will and won't collapse:
- A zero-length edge whose *surviving* parent is `Sample`-typed: refuses
  (the assert), because collapsing would silently overwrite the `Sample`
  tag with the child's type — a correctness guard, not a limitation to
  route around.
- A zero-length edge whose surviving parent is plain `Node`-typed, even
  with 2 children **each independently having their own further
  descendants**: collapses correctly, with no polytomy — verified directly
  by round-tripping a hand-built `Root(0)->Sample(1)->Node(3){Node(3)
  ->Sample(5), Node(3)->Sample(7)}` tree through `newick`/`parse_newick`
  and confirming the middle `Node(3)` ends up with exactly its 2 intended
  grandchildren, not a flattened merge.

This directly implies the fix's shape: **never create a zero-length edge
below a `Sample` node; zero-length edges below a plain `Node` are fine.**
Final design (`src/simulate.jl`'s `apply_event!`):
- **SAMPLE**: create a fresh, correctly-timestamped `Sample` node as
  `cur`'s child. If non-destructive, that new node **itself** becomes the
  new open tip in the inventory — no separate continuation node at all
  (the zero-length-below-Sample case above is why not).
- **BIRTH**: unchanged in the common case. New exception: if the chosen
  lineage's current node is itself `Sample`-typed (its previous event was
  a non-destructive sample, and this birth is its very next event), a
  plain `Node` is interposed first to hold the bifurcation, since a
  `Sample` node may have at most 1 child
  (`NaiveSEIR.singular_part!`'s `@assert length(n.children)<2`, and
  empirically confirmed no R-emitted `Sample` node ever has 2 children
  either). This intermediate's own zero-length edge to its 2 children is
  the *safe* case verified above.

## 4. Validation that the fix is actually correct (not just "different")

1. **Round-trip self-consistency**: over 314 non-extinct draws, direct
   in-memory `first-sample-time` vs the same statistic recomputed after a
   full `newick`/`parse_newick` round trip now agree to ~1e-6 (floating
   rounding only, from the writer's 6-significant-digit branch-length
   rounding) — before the fix, only 3/314 were even exactly equal, and the
   mean shifted by a full 0.085 across the round trip (a live symptom of
   the same bug, independent of any R comparison).
2. **No regression**: full suite (`RUN_HEAVY_TESTS=no julia --project=.
   test/runtests.jl`) still **3301/3301**, including Milestones 1–2's own
   simulate-then-filter acceptance tests, both of which are somewhat
   insensitive to fine intra-branch timing (they check finite likelihood
   and structural invariants, not exact times) and so did not themselves
   catch this bug — a fact worth flagging for future work (see §7).
3. **The actual cross-validation, before vs after** (N=400, same seeds
   both times):

   | statistic | before fix | after fix |
   |---|---|---|
   | first-sample-time KS D, p | 0.24, ≈0 (**significant**) | 0.076, 0.34 (not significant) |
   | first-sample-time mean (Julia / R) | 1.41 / 1.95 | 1.86 / 1.95 |

## 5. Final N=1000 cross-validation result

Command sequence (fully reproducible, fixed seeds both sides):
```
Rscript scripts/seir_crossvalidate.R 1000 seir_crossvalidate_r_trees.txt
julia --project=. scripts/seir_crossvalidate.jl 1000 seir_crossvalidate_r_trees.txt
```
Full report: `seir_crossvalidate_results.md` (repo root). Summary:

| statistic | Julia | R | test | result |
|---|---|---|---|---|
| extinction fraction | 0.243 | 0.237 | two-proportion z | z=0.31, not significant |
| nsample | mean 89.8 | mean 89.5 | KS | D=0.036, p=0.54 |
| first-sample-time (nonextinct) | mean 1.79 | mean 1.84 | KS | D=0.052, p=0.26 |

All three checks: no evidence of a distributional difference. Per the
extinction-spike caveat in the script's own header (repeated here since
it matters for interpreting this table): the KS tests are *conservative*
in the presence of a shared ~24% mass point at `nsample=0` on both sides
(the ECDFs agree exactly there, deflating D) — which is exactly why the
extinction fraction is checked as its own separate two-proportion test,
not folded into the KS read.

`seir_crossvalidate_r_trees.txt` (the raw N=1000 R-generated Newick
dump, ~6.5MB) was **not** committed — fully regenerable from the command
above with the fixed seed (`set.seed(20260813)` in the R script). It was
moved to the session scratchpad rather than left at the repo root.

## 6. Design/placement decisions

- **E0=0 parameter choice**: R's `rinit` (`yaml/seir.yml`) grafts
  `round(pop*E0)` lineages into Exposed and `round(pop*I0)` into
  Infectious. With `E0=0`, this reduces to exactly Julia's single-root
  `graft=[0,1]` case. Using R's own default `E0=0.05` would require
  exercising the multi-root/forest path flagged as untested in
  `check_milestone1.md` §6 — deliberately avoided here to keep this
  milestone isolated to the engine question it was designed to test, not
  compounded with a second unverified code path.
- **Placement**: `scripts/seir_crossvalidate.{R,jl}`, matching the
  existing `scripts/gen_seir_for_crosscheck.R`/`julia_cblv_crosscheck.jl`
  convention. Deliberately **not** wired into `test/runtests.jl`: it
  requires a working R + phylopomp installation, and a stochastic KS test
  has no place in CI.
- **No new dependency**: the two-sample KS test is hand-rolled in the
  Julia script (standard asymptotic Kolmogorov formula) rather than adding
  `HypothesisTests` to `Project.toml` for a one-off validation script.

## 7. Suggested reviewer focus

1. **The fix itself** (`src/simulate.jl`'s `apply_event!`) — re-derive
   independently that a `Sample`-typed node can legitimately end up as the
   `cur` for a later `BIRTH`, and that the interposed-`Node` handling is
   correct, rather than trusting this document's trace.
2. **Why Milestones 1–2's own tests didn't catch this** — their acceptance
   criterion (finite filter likelihood) is a real but coarse signal;
   consider whether a cheap, permanent regression test asserting
   `slate(sample_node) == firing_time` at the unit level (not just via
   full-scale KS comparison) is worth adding to `test/seir_simulate.jl`
   going forward, so this class of bug is caught at Milestone-1 speed next
   time rather than requiring a full R cross-validation run.
3. **The `clip_zlb!` boundary** (§3) — this is pre-existing code, not
   written this session; confirm the empirical characterization here
   (collapses fine below `Node`, refuses below `Sample`) actually matches
   its written logic rather than just this session's two test cases.
