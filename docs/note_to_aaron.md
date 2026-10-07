# Two likelihood defects in the guided MERS and SEIR filters

To: Aaron King. From: Deepan Islam. 2026-10-07.
Applies to `origin/mers_guide` at c80a4ce (`src/examples/mers_guided.jl`, `seir_guided.jl`, `src/guide.jl`).

## Summary

Two functions in the guided filters give log-likelihood estimates that are too high or too low.

1. `regular_transmission_cc!` and `regular_transmission_hh!` in `mers_guided.jl` return a log-weight of zero. The target weight is `log(1 - C(ℓ,2)/C(I,2))`, so the estimate is too high. On 2026-09-18 the difference on one 4-tip tree was about 6 log-units.
2. `regular_transmission_hc!` and `_ch!` in `mers_guided.jl` and `regular_transmission!` in `seir_guided.jl` draw the outcome "no lineage moves" with probability proportional to `n − ℓ`, through `choose_branch` (`guide.jl:253`). The weight is zero when every host in the source deme is tracked (`ℓ = n`), although the target mass of that outcome is positive. The estimator therefore omits part of the likelihood. On six small simulated trees its estimates were 0.14 to 0.78 log-units lower (MERS) and 2.3 to 2.6 log-units lower (SEIR) than those of the reference filter described below.

I have not reported item 2 before. Item 1 is unchanged since 5861849. A patch for item 2 is in `aaron_choose_move.patch`, made against f501185 in my checkout of your repository.

## Item 1: missing no-fork factor

In an interval with no observed genealogy event, a transmission within one deme must not have joined two of the `ℓ` tracked lineages. Its target weight is `1 − C(ℓ,2)/C(I,2)`, with `I` counted after the event. This equals `Φ_id + ℓ·Φ_inl`, the sum over saturations of the event's production slots. `mers_naive.jl` applies it. Commit f501185 in my checkout adds it:

```julia
log(1 - ellC*(ellC-1)/Ic/(Ic-1))     # regular_transmission_cc!, Ic after the increment
```

`regular_transmission_hh!` is symmetric.

## Item 2: zero proposal probability for "no move"

For a transmission from deme `i` to a new host in deme `j`, let `ℓ_i` and `ℓ_j` be the tracked lineages, and `n_i` and `n_j'` the host counts after the event. The target factors are:

- No lineage moves: `f0 = 1 − ℓ_j/n_j'`. This sums two cases, an untracked parent and a tracked parent that keeps its own lineage (`Φ_id + ℓ_i·Φ_inl`).
- Tracked lineage `b` moves from `i` to `j`: `f1 = (1 − (ℓ_i − 1)/n_i)/n_j'`.

`f0 > 0` whenever `ℓ_j < n_j'`, including `ℓ_i = n_i`. `choose_branch` draws "no move" with probability `(n − ℓ_i)/s`, where `s` also contains the relative hazards of the tracked lineages (`guide.jl:253-267`). At `ℓ_i = n_i` it never draws that outcome, so the term `log(f0) − log(q)` is never added. The expectation of the importance weights is then below the likelihood by the mass of the omitted outcome.

The fix draws the outcome with `w0 = f0` for "no move" and `wmove = f1` for each moved lineage, multiplied by `relhaz`:

```julia
b, q = choose_move(t, guide, node, cols, i, j, f0, f1)   # q == 0: return -Inf
ll   = b == 0 ? log(f0) - log(q) : log(f1) - log(q)      # after swap! if b != 0
```

`choose_move` is about 20 lines in `guide.jl` and is included in the patch. SEIR progression (`regular_progression!`, which calls `choose_branch(..., E, cols, Expos, Infec)`) needs no change. The progressing host carries its lineage, so "no move" requires an untracked host and its target mass is zero at `ℓ = n`.

## Evidence

Each entry is the `logmeanexp` of 40 independent filter runs with Np = 5000 (`trigger = 0.2, target = 0.8`), with its Monte Carlo standard error. The trees were simulated from the model tables, have one root and 3 to 5 tips. The reference filter is a filter with a naive proposal, built from the event table. Its singular steps agree with `terminal_sample!`, `inline_sample!`, `singular_branch!` and `singular_root!` to 1e-15 on 770 outcome comparisons (320 MERS, 450 SEIR). "Patched" is `origin/mers_guide` with the attached patch applied and nothing else changed.

| Model, tree | Reference filter | `origin/mers_guide` | Patched |
|---|---|---|---|
| MERS 1 (3 tips) | −1.228 ± 0.009 | −1.363 ± 0.007 | −1.177 ± 0.035 |
| MERS 2 (4 tips) | −3.180 ± 0.017 | −3.425 ± 0.017 | −3.151 ± 0.037 |
| MERS 3 (4 tips) | −8.312 ± 0.057 | −9.095 ± 0.122 | not run |
| SEIR 1 (3 tips) | −10.30 ± 0.05 | −12.55 ± 0.03 | −10.25 ± 0.08 |
| SEIR 2 (4 tips) | −15.13 ± 0.23 | −17.51 ± 0.09 | −15.39 ± 0.13 |
| SEIR 3 (5 tips) | −19.01 ± 0.19 | −21.58 ± 0.08 | −19.08 ± 0.12 |

The guide for MERS was `fsmarkov(Camel=>0.5, Human=>0.5, (Camel,Human)=>0.01)`, and for SEIR `fsmarkov(Expos=>0.1, Infec=>1, (Expos,Infec)=>1)`. The MERS columns include the item 1 fix (f501185), so the MERS differences isolate item 2.

In the five trees run with the patch, the patched estimate lies within 0.3 of the reference filter. In all six trees the unpatched estimate lies 0.14 to 2.6 below it. The largest patched difference is SEIR tree 2, 0.26, about 1.0 combined standard errors. The patch changes only the proposal, and the reference filter shares no proposal code with the guided filters. The agreement is therefore consistent with item 2 as the cause of the difference. It does not rule out other, smaller differences.

## Limits

- Six simulated trees, each with at most 5 tips, one parameter set per tree, one guide per model.
- The standard error of `logmeanexp` understates the uncertainty when the weights are heavy-tailed. On 6 to 9 tip MERS trees with the guide above, the runs were too noisy to compare, so they are not reported.
- The patch changes only the proposal of the three functions named above. I did not run your test suite with it.
- Item 1's size (about 6 log-units) comes from one tree measured on 2026-09-18, before the other changes described here.
- A regression test for item 2: set `ℓ_i = n_i`, call `regular_transmission_hc!` 10⁴ times, and require that "no move" occurs at least once.
