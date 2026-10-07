# Handoff — M15: generic `singular_update!` (MERS and SEIR), checked against Aaron's `mers_guided.jl` and `seir_guided.jl`

## Status
MERS and SEIR singular parts done and checked. `apply_move!`, `kli_select`, `kli_decay` are still stubs.
Nothing committed.

## What changed
- `src/examples/mgp_filter.jl`: `singular_update!(cols, geneal, node, x, θ, model) -> (Δll, x′)` replaces the stub
  (new `geneal` argument; nothing referenced the old signature). Root, Sample (no children), Node (fork).
  `_orderings` gives the distinct child orderings of an event's slot demes.
- `src/examples/mgp_mers_filter.jl`, `src/examples/mgp_seir_filter.jl`: `mers_compiled_filter_pomp` and `compiled_filter_pomp` take `generic_singular = false`. Default paths are unchanged.
- `test/kli_singular_update_test.jl`, wired in `test/runtests.jl` after `kli_mers_compiled_test.jl`.

Fork proposal weight is `α_u / (number of distinct orderings of u)`, with outcomes in event order and then lexicographic
ordering, one `rcateg` call. Charge `log α_u + log φ_u − log p`, φ_u from `full_transitions` at the post-fork ℓ and post-event n.
Sample: the deme is `n.deme`, or the deme holding the lineage when `n.deme` is missing (SEIR). Candidates are the SAMPLE events with
`from` that deme (`sum(r) == 1` for a sampled ancestor). With several, one is drawn with weight α. Charge `log α + log φ − log q`;
φ = 1 for a destructive sample, `(n − ℓ_post)/n` for a non-destructive tip, `1/n` for a sampled ancestor. The two non-destructive φ
are the rules in `NaiveSEIR.singular_part!`; the IR does not derive them.

## Reference used
`../PhyloPOMP.jl_aarons`, branch `fix/guided-mers-noop-factor`, commit f501185 (authored dislam5, one ahead of
`origin/mers_guide` 5861849). 5861849 still has the missing no-fork factor in `regular_transmission_cc!/_hh!`; f501185 adds it.

## Evidence
1. Exact Rationals, (ℓ, n) grid, 10,736 checks, 0 failures: TCC/THH no-fork `Φid + ℓΦinl = 1 − C(ℓ,2)/C(I,2)`;
   THC/TCH no-move `= 1 − ℓ_other/I_other`; cross `Φcr = (1 − ℓ′_src/I_src)/I′_dst`. These are Aaron's regular-move terms.
2. Per-outcome targets (`ll + log p`) of Aaron's `singular_branch!` and `terminal_sample!` run on his guide, against
   the generic code: 320/320 outcomes (6 fork outcomes + 2 sample demes, 40 states), worst |Δ| 8.9e-16.
3. Generic vs `NaiveMERS.singular_part!` inside the compiled filter, seed-matched: 120/120 Np=200 log-likelihoods,
   worst |Δ| 3.6e-15 (20 simulated trees, 3 to 10 tips, B = 0 and B > 0).
4. Filter-level, Np=5000 × 40 reps, three 3–4 tip simulated trees (logmeanexp ± SE):

   | tree | generic singular | our `mers_guided.jl` | Aaron f501185 | Aaron f501185 + our hc/ch patch |
   |---|---|---|---|---|
   | 1 | −1.228 ± 0.009 | −1.200 ± 0.016 | **−1.363 ± 0.007** | −1.177 ± 0.035 |
   | 2 | −3.180 ± 0.017 | −3.200 ± 0.020 | **−3.425 ± 0.017** | −3.151 ± 0.037 |
   | 3 | −8.312 ± 0.057 | −8.364 ± 0.128 | **−9.095 ± 0.122** | not run |

## SEIR
- Per-outcome targets of Aaron's `terminal_sample!`, `inline_sample!`, `singular_branch!` (`seir_guided.jl`) against the generic code:
  450/450 match, worst |Δ| 4.4e-16 (both fork orderings, destructive and non-destructive tip, sampled ancestor, ℓ_I − 1 = 0, 1, 2).
- Generic vs `NaiveSEIR.singular_part!` inside `compiled_filter_pomp`, seed-matched Np=100: 212/212 finite match, worst |Δ| 5.3e-15
  (27 trees with 22 sampled ancestors; χ = 0 and χ > 0).
- Filter level, Np=5000 × 40, trees of 3, 4, 5 tips, guide `fsmarkov(Expos=>0.1, Infec=>1, (Expos,Infec)=>1)`:

  | tree | generic singular | our `seir_guided.jl` | Aaron `seir_guided.jl` | Aaron + our `regular_transmission!` |
  |---|---|---|---|---|
  | 1 | −10.30 ± 0.05 | −10.24 ± 0.08 | **−12.55 ± 0.03** | −10.25 ± 0.08 |
  | 2 | −15.13 ± 0.23 | −14.92 ± 0.26 | **−17.51 ± 0.09** | −15.39 ± 0.13 |
  | 3 | −19.01 ± 0.19 | −18.70 ± 0.21 | **−21.58 ± 0.08** | −19.08 ± 0.12 |

  Aaron's SEIR guided filter has the same no-move zero-support problem as MERS (`regular_transmission!` uses `choose_branch`), and the gap is
  larger, 2.3 to 2.6 log-units here.

## Finding: a third upstream difference (MERS)
Aaron's `regular_transmission_hc!/_ch!` use `choose_branch`, whose no-move weight is `n − ℓ`. When every host in the source
deme is tracked (`ℓ = n`) that is zero, so "tracked parent keeps its own lineage" is never proposed. Its target mass
`ℓ·Φinl` is positive (identity 1 above), so the estimator is biased low. Our copy uses `choose_move` with the target factors
`f0`, `f1` (`src/guide.jl`). Patching only those two functions and adding `choose_move` to a scratch copy of Aaron's repo moves
his estimates onto ours (table, last column). The singular functions and the f501185 no-fork fix agree with ours.
To relay to Aaron together with the two bugs in memory `aarons-mers-guided-rewrite-bugs`.

## Not covered
- Real 274-tip tree (tied camel tip times give −Inf in every kernel, see `mers_kernel_diagnostic_findings.md`).
- SEIR `apply_move!`/regular part: still the compiled version; only the singular part is generic.
- Aaron's guide was `m0 = fsmarkov(Camel=>0.5, Human=>0.5, (Camel,Human)=>0.01)`, uncalibrated; at 6 to 9 tips his estimator
  is heavy-tailed and its SE is unreliable (trees 1 and 4 of the 6 to 8 tip set). The 3 to 4 tip set above is the usable comparison.

Scripts are in the session scratchpad only (`pernode_aaron*.jl`, `pernode_generic*.jl`, `aaron_run.jl`, `ours_run.jl`, `seir_*_run.jl`, `gate*.jl`).
