# Review: Gillespie simulator and docstrings (2026-10-05)

Branch `simulator-phase2`, working tree as of this date. Nothing in `src/` or `test/` was edited.

## Page 1

**First action:** done 2026-10-05. Findings 1, 2 and 3 are applied (see below).

| Item | Verdict |
|---|---|
| Engine branches on model name | No. Only on `EventType` (5 values). Pass. |
| Models run through one engine | SEIR, SIR, SI2R, MERS, and the test models `Division`/`Twins`. MTBD also runs, declared ad hoc (finding 2). Pass. |
| Line-level correctness | No bugs found in `apply_event!`, `prune!`, `rcateg`, `repair!`. |
| Simulator tests | `simulate_checks`, `sir_simulate`, `seir_simulate`, `mers_simulate`, `mtbd`: 6505 of 6505 pass (69 s). |
| Engine checks the IR it is given | Yes, since fix 1 (2026-10-05). Before that, only `@mgp` checked it. |
| Tests drive the real loop | Partly. `simulate_checked` re-implements the loop (finding 3). |
| Docstrings | Rewrite: 13 src files (6 whole-file, 7 header-only) and 6 test files. About 525 lines to delete and about 2,650 to cut down to a few hundred. |

---

## Part 1: Simulator

### What passes

The design matches the "inspectable IR + generic engine" target. `@mgp` produces an `MGPModel`, a plain table of `Event`s. `_simulate_loop!` evaluates `ev.hazard`, picks one event with `rcateg`, adds the precomputed `δ` tuple, and calls `apply_event!`. `apply_event!` branches only on `ev.type`. No code path names a model.

Checked and correct:

- `BIRTH`: one new `Node` at time `t`. `r[j]` copies are opened in deme `j`. The drawn slot is reused when `r[from] ≥ 1`. `fork(A => B, B)` (the parent leaves its deme) works.
- `SAMPLE`: the `Sample` node becomes the open tip when `r[from] == 1`, so no zero-length edge sits below a `Sample`.
- `prune!`: the parent-first pass splices a whole chain of degree-1 `Node`s in one pass, because a spliced child's `parent` is rewritten before that child is visited.
- `rcateg`: the remainder past a zero last weight, and a zero draw on a zero first weight, both land on a positive-weight category.
- `repair!`: the (time, depth, name) order is a valid strict order. In simulated trees, same-time nodes occur only as founder roots.
- Per-step invariant: lineage count equals compartment count in every deme. A `Δ`/`r` disagreement on a deme is caught on its first firing.

### Findings, ranked

**1. The engine trusts the IR it is given.** All structural checks live in `@mgp` (`src/examples/mgp_macro.jl`). `validate_model` (`src/examples/mgp_audit.jl:171`) exists, but `simulate` never calls it. A hand-built `MGPModel` (for example `SEIR_REFERENCE`) skips every check. Probe results (`docs/review_simulator_2026-10-05_probe.jl`; its four bad `Event`s are ready-made `@test_throws` cases):

| Bad `Event` | What `simulate` does |
|---|---|
| `BIRTH` with `sum(r) == 1` | Accepted. Returns a tree. The birth leaves no trace. |
| `Δ = [:I=>-1, :Rr=>+1]` (typo for `:R`) | Accepted. `:Rr` is dropped: `I` goes 3 → 0 and `R` stays 0. |
| `SAMPLE` with `r[from] == 2` | Accepted. Treated as non-destructive. |
| `from = 2` in a 1-deme model | `BoundsError` from inside `SimInventory`, not an `ArgumentError`. |

Fix (about 30 min including tests):

1. At the top of `_simulate`, run `issues = validate_model(model)`. If `issues` is non-empty, throw `ArgumentError(join(issues, "\n"))`.
2. Add to `validate_model`: `BIRTH ⇒ sum(r) == 2`; `SAMPLE ⇒ r[from] ∈ {0, 1}`; `MIGRATION ⇒ length(into) == 1`; and move-vs-pop agreement. That last check is `_check_move_vs_pop`'s rule: `r[d] - (d == from)` equals the `Δ` entry for deme `d`.
3. Make `@mgp` call the same function, so the rules live in one place.

**Applied 2026-10-05 (steps 1 and 2).** `validate_model` now also checks, per type: `BIRTH` has 2 products; `MIGRATION` has one `into` and `r` is 1 there; `DEATH` has `r` all 0; `SAMPLE` has `r[from]` 0 or 1 and 0 elsewhere; `NEUTRAL` has `from = 0`. It also checks Δ against `r` and `from` on every deme, and that no compartment appears twice in Δ. `_simulate` throws `ArgumentError` listing every issue. All four probe models are now rejected. Tests were added in `test/population_ir_test.jl` (one per rule) and `test/sir_simulate.jl` (the four probes, plus `SEIR_REFERENCE` as a hand-built model that still runs).

Step 3 was not done. `@mgp SEIR` is expanded inside `mgp_macro.jl`, which loads before `mgp_audit.jl`, so `validate_model` does not exist yet at that point. The macro's errors are also phrased in its own syntax (`fork`, `pop`, `move`). The macro already enforces the same rules on everything it emits, so `validate_model` only has to cover hand-built models. Sharing one function would first need the IR moved out of `src/examples/` (next paragraph).

Related: the IR types (`Event`, `MGPModel`, `@mgp`, `validate_model`) are under `src/examples/`. `src/simulate.jl` depends on them, and `Examples.jl` is included before `simulate.jl`. The IR belongs in `src/`, with only the model declarations in `src/examples/`.

**2. MTBD is the missing model, and it declares cleanly.** The MTBD filters (`mtbd_naive.jl`, `mtbd_guided.jl`) have no `@mgp` declaration. The probe declares one with no engine change. It uses 10 events (it leaves out the two `r < 1` sampling events): `fork(I_1 => I_1, I_2)` with `pop=(I_2=+1)`, `sample_remove`, `chop`, and `swap`. It simulates a 69-sample tree (θ in the probe, `tmax = 6`). Recommended: add `@mgp MTBD` to `src/examples/`. Then add an `mtbd` case to `scripts/mgp_crossvalidate.{R,jl}` using R phylopomp's `runMTBD2`. That is the simulator check. MTBD's exact ODE likelihood checks the filters, not the simulator.

The cross-validation covers only `seir`, `seirchi` and `mers`. R phylopomp also exports `runSIR`, `runSI2R`, `runBDEI` and `runBDSS`. SIR and SI2R are already declared with `@mgp`, so they can be added with no new model code.

**Applied 2026-10-05.**
- `src/examples/mgp_mtbd.jl` declares `@mgp MTBD`. It has 12 events, the same as R's MTBD2, including `r1`/`r2` (the probability that a sample removes the host). It is exported and passes `validate_model`.
- `test/mtbd_simulate.jl` runs the shared invariants (one founder, and a forest, with `r1 = 0.7`). It also passes one simulated tree (`r1 = r2 = 1`) through `GuidedMTBD.filter_pomp` and gets a finite log-likelihood.
- `scripts/mgp_crossvalidate.{R,jl}` gained `sir`, `si2r` and `mtbd` cases, a `SEED` override, and per-deme sample counts for any deme set.

Results at N = 2000 per side:

| Model | Statistics | p < 0.05 |
|---|---|---|
| SIR | 11 | 0 |
| MTBD | 13 | 0 |
| SI2R | 14 | 2 (mean internal-node time 0.029, total branch length 0.017). With `SEED=7` and N = 6000: 0. |

The bundled `SI2R` sampled without removing the host (`move=sample`), as in `si2r_model.qmd`. R's `runSI2R` removes it (`sample_death`). `SI2R` now has a parameter `r`, the probability that sampling removes the host, as MTBD has: `r = 0` is the qmd model and `r = 1` is R's. The `si2r` case runs `SI2R` with `r = 1`. It reproduces the output of the earlier script-local copy exactly. `runBDEI` and `runBDSS` were not added.

Expressibility limits worth knowing for the acid test:
- Importation of an infected host from outside the model is rejected (`move=none pop=(I=+1)` fails `_check_move_vs_pop`). The engine could support it by opening a new `Root` mid-run.
- Compound events (KLI §3.2(f)) are not representable.

**3. The per-event test copies the engine loop.** `simulate_checked` (`test/simulate_checks.jl:19-90`) re-implements the Gillespie loop and uses `apply_delta`, the `Dict` path. So the real loop's parts are only tested end to end: the `δ`-tuple update, the hazard guards, the negativity check, and the inventory check. If `_simulate_loop!` drifts, the per-event checks keep passing against the copy. Fix: generalize the `rec` argument of `_simulate_loop!` into a per-event callback `(G, inv, ev, t, x)`. `simulate_trajectory` then becomes one such callback, and `simulate_checked` passes another. After that, `apply_delta` is used only in `test/sir_simulate.jl:162`. `mgp_filter.jl` defines a duplicate, `apply_pop`.

**Applied 2026-10-05.** `_simulate` takes an `onevent` function in place of `rec`. It is called once before the first event (`ev = nothing`) and after every event, with `(G, inv, ev, t, x)`, and draws nothing from `rng`. `simulate_trajectory` records through it. `simulate_checked` now runs `PhyloPOMP._simulate` with a checking function instead of its own copy of the loop. It also checks every new state against `apply_delta`, so the `Dict` path and the `δ`-tuple path are compared on every event.

**4. Each hazard call is dynamically dispatched.** `Event.hazard::Function` is an abstract field (`isconcretetype` returns `false`), so every hazard evaluation is a dynamic call. That is K calls per step. This is a cost, not a bug: the SIR scaling test runs in under 5 s. If it matters later, store the hazards as a tuple type parameter on `MGPModel`.

### Beyond the simulator

The filter layer does not pass the same acid test yet. In `mgp_filter.jl`, the generic `kli_select`, `kli_decay`, `apply_move!` and `singular_update!` all call `error()` and have no callers. `mgp_seir_filter.jl` and `mgp_mers_filter.jl` reuse `NaiveSEIR.singular_part!` and `NaiveMERS.singular_part!`. This is the `@proposal` half of the design in the pasted note.

---

## Part 2: Docstring sweep

Rubric applied by four reviewers (read-only):

- Keep what a caller needs: signature, arguments, return value, what throws, one-line invariant, short paper reference.
- Cut: history, milestone/Gate/handoff/scratch references, "VERIFIED"/"STATUS"/"confirmed", line numbers into other files or the tex, ALL-CAPS, banners, self-justification ("documented, not hidden", "cited inline", "deliberately"), comments that restate the code, derivation essays.

### Verdict by file

| Verdict | src | test |
|---|---|---|
| **rewrite** | `mgp_filter_ir.jl` (277 of 348 lines are comments), `mgp_decay.jl` (224/281), `mgp_mgpaudit.jl` (217/304), `mgp_mers_filter.jl` (269/511), `mgp_seir_filter.jl` (199/400), `mgp_transitions.jl` (170/253). Header essays only: `mgp_reduce.jl`, `mgp_phi.jl`, `mgp_proposal.jl`, `mgp_filter.jl`, `mers_funs.jl`, `mers_soft.jl`, `mers_hard.jl` | `kli_kingman_moran_test.jl`, `kli_full_transitions_test.jl`, `kli_properties_test.jl`, `kli_filter_ir_test.jl`, `kli_mers_compiled_test.jl`, `kli_proposal_test.jl` |
| **trim** | `mgp_explain.jl`, `mgp_audit.jl`, `mgp.jl`, `mgp_si2r.jl`, `mgp_sir.jl`, `mtbd_funs.jl`, `mtbd_naive.jl`, `mtbd_guided.jl`, `mers_naive.jl`, `mers_guided.jl`, `guide.jl`, `genealogy.jl`, `simulate.jl`, `rcateg.jl`, `parse.jl`, `fsmarkov.jl` | `kli_si2r_test.jl`, `kli_phi_test.jl`, `kli_reduce_test.jl`, `kli_decay_test.jl`, `kli_seir_compiled_test.jl`, `mgpaudit_test.jl`, `population_ir_test.jl`, `mers_simulate.jl`, `mers_soft.jl`, `mtbd.jl`, `seir_macro_equivalence.jl`, `seir_simulate.jl`, `sir_simulate.jl`, `rcateg.jl` |
| **acceptable** | `mgp_macro.jl`, `mgp_mers.jl`, all `seir_*.jl` filters, `mers_tree.jl`, `Examples.jl`, `PhyloPOMP.jl`, `coloring.jl`, `demes.jl`, `indicator.jl`, `newick.jl`, `cblv.jl` | `simulate_checks.jl`, `cblv.jl`, `fsmarkov.jl`, `guide.jl`, `newick.jl`, `parse.jl`, `runtests.jl`, `seir_{naive,soft,guided,hard}.jl`, `mers_{naive,guided,hard}.jl` |

Line-by-line cuts with replacement text: `docs/review_docstrings_2026-10-05_detail.md`.

**Applied 2026-10-05.** 52 files, 3,154 lines deleted and 586 added. `src` went from 9,250 to 7,410 lines and `test` from 6,035 to 5,308. A parse-level comparison found no code change other than the planned strings and the removed `@test true`, `demes_sampled` and load-time asserts. The regression grep below now gives 0 hits.

Approximate totals: about 525 lines to delete. About 2,650 lines sit in blocks that shrink to 1–4 lines each.

### Worst blocks

1. `src/examples/mgp_seir_filter.jl:31-117`: the 87-line "decay discrepancy" essay. It names scratch probe files, says "ALL MATCH: true", and cites tex and naive line numbers. Keep one line in the `compiled_decay` docstring: decay = λ + Σ leftover over DEATH events = γ·ℓ.
2. `src/examples/mgp_filter_ir.jl:1-125`: rename history ("CORRECTION (M06, Part 0)"), milestone scope lists, and handoff paths.
3. `src/examples/mgp_mers_filter.jl:27-178`: FINDING 1/2/3. These are verification claims, a bug narrative, and a workaround with its resolution.
4. `src/examples/mgp_decay.jl:82-160`: "DISCREPANCY FOUND … reported, not silently resolved". It contradicts item 1, which says the question is resolved.
5. `test/kli_kingman_moran_test.jl:1-74`: the M00–M12 narrative and a textbook derivation.
6. `src/examples/mers_soft.jl:49-63`: a proof that a "former on-disk version" equals the current code. `mers_hard.jl:46-52` points back to it.

### Non-comment problems the sweep found

- `test/kli_full_transitions_test.jl:237`: `@test true  # documents the cross-check above`. This test checks nothing.
- `src/examples/mgp_mers_filter.jl:442` says it matches "`mgp_seir_filter.jl`'s `mers_compiled_filter_pomp`". That function is defined only in `mgp_mers_filter.jl`.
- `src/examples/mgp_phi.jl:111-121`: `@assert`s run at package load act as tests. Move them to `test/kli_phi_test.jl`.
- `test/seir_hard.jl:17`: the testset is named "SEIR model with guided proposals" but tests the hard kernel.
- Printed strings carry milestone tags: `mgp_mgpaudit.jl:164-176, 243, 251`, `mgp_explain.jl:248`, and the `@testset`/`@info` titles ending in `(Mxx)` in the `kli_*` tests. Check whether any test matches the text before changing them.
- `test/kli_si2r_test.jl` header contains an absolute home-directory path.
- `test/kli_mers_compiled_test.jl:137`: `demes_sampled` is never used.

### Regression grep after cleanup

```bash
grep -rnE 'cited inline|VERIFIED|STATUS:|RESOLVED|SCOPING|[Mm]ilestone|\bM[01][0-9]\b|handoff|scratch|probe|Gate [0-9]|prior version|was fixed|documented, not hidden|^# ={10,}' src test --include=*.jl
```

At review time this gives 405 hits in 26 files. The worst are `mgp_filter_ir.jl` (41 hits), `mgp_mers_filter.jl` (39), `mgp_seir_filter.jl` (38), and `mgp_mgpaudit.jl` (36).

Line ranges above are as of 2026-10-05. They shift after the first edit, so apply the cuts bottom-up within each file, or re-read each file before editing it.
