# Handoff: M16, generic regular step (`regular_step!`, naive proposal) for MERS and SEIR

## Status
Done for MERS and SEIR. With the naive proposal, the generic regular step reproduces `mers_compiled_regular_part!` and
`compiled_regular_part!`: same final state, same coloring, same RNG stream, ll within 1e-13.
The default filter paths are unchanged. Nothing committed.

## What changed
- `src/examples/mgp_filter.jl`: the stubs are replaced.
  - `kli_slots(ev)` and `kli_slots(model)`: a regular BIRTH or MIGRATION whose `r` has a slot in a deme `d ≠ from`
    gets two slots, `[0, d]`: no move, then move to `d`. Every other regular event gets one slot. Over the
    model's regular events in table order this gives the compiled slot order (MERS 12 slots, SEIR 6).
    It throws `ArgumentError` for three event kinds: a BIRTH with `r[from] == 0`, a BIRTH with products in two
    other demes, and a regular SAMPLE.
  - `kli_select(ev, cols, x, model) -> Vector{Float64}` returns the naive π for each slot:
    - BIRTH split: the target share `T0/(T0+T1)`, from `full_transitions`. Here `T0 = Φ_id + ℓ_a·Φ_inl` and
      `T1 = ℓ_a·Φ_cr`, with Φ_cr computed after one lineage moves. This is `no_move_share`.
    - MIGRATION: `[1 − ℓ_a/n_a, ℓ_a/n_a]`.
    - DEATH: `[1 − ℓ_a/n_a]`, or 0 when `n_a = ℓ_a`.
    - Everything else: `[1]`.
  - `kli_decay(alpha, pi, slots, cols, x, θ, model)` is `total_decay` plus, for each regular event, its target
    exit rate minus `α·Σπ`. The target exit rate is `α·1{n>ℓ}` for DEATH and `α` otherwise. With naive π this
    reduces to `compiled_decay`.
  - `apply_move!(cols, ev, slot, x1, model; rng)` returns `log Φ − log q`:
    - Slot 0: `log(Φ_id + ℓ_a·Φ_inl)`.
    - Move slot: draw `b = rand(rng, cols[from])`, `swap!` it, then return `log Φ_cr − log(1/ℓ_pre)`.
    - DEATH and NEUTRAL: 0. The DEATH weight `1/π` enters through the step's `−log π`.
  - `kli_rates!` fills the per-slot α and π and returns the decay.
  - Helpers: `_deme_counts`, `_phi`, `_no_move_target`, `_cross_phi`. Nothing is looked up by event name.
  - `regular_step!(cols, ll, t, dt, x, model, θ; rng = default_rng()) -> (ll, x, t + dt)` makes its calls in the
    compiled order: `rcateg`, then always `step = −log(rand())/s` (even when `s = 0`), then
    `ll −= decay·step + log π_k`, then `apply_pop`, then `apply_move!`.
- `src/examples/mgp_mers_filter.jl`: `mers_compiled_filter_pomp(...; generic_regular = false)`, plus `_mers_theta(args)`.
- `src/examples/mgp_seir_filter.jl`: `compiled_filter_pomp(...; generic_regular = false)`.
- `test/kli_regular_step_test.jl`: new file, registered in `test/runtests.jl` after `kli_singular_update_test.jl`.
- `test/kli_filter_ir_test.jl`: line 212 asserted that `kli_decay` throws (the old stub). It now checks that
  `kli_decay` ≈ `compiled_decay` at one SEIR state.

## Evidence (all run on 2026-10-07)
`julia --project=test <driver> kli_regular_step_test.jl` passed 25/25 in 28 s. The driver defines `h1`/`h2` and includes the file.

| check | result |
|---|---|
| MERS per-outcome weights and π against closed forms, grid `I_c, I_h ∈ 1:6`, all ℓ | 19,764 checks, 0 failures; 324 states with ℓ_src = n_src |
| SEIR per-outcome (infection, progression, recovery), grid `E ∈ 0:6, I ∈ 1:6`, all ℓ | 4,854 checks, 0 failures; 168 states with ℓ_I = I |
| `kli_decay` vs `compiled_decay`, 300 random states per model | worst \|Δ\| 3.6e-15 |
| isolation vs `mers_compiled_regular_part!`, 400 draws | 0 mismatches; worst \|Δll\| 9.9e-14; a lineage moved in 250 |
| isolation vs `compiled_regular_part!`, 400 draws (ψ, χ sometimes > 0, E = 0 allowed) | 0 mismatches; worst \|Δll\| 4.3e-14; a lineage moved in 311 |
| MERS Np=100 seed-matched, default vs `generic_regular` and vs generic regular + singular | 6 trees, 23 finite of 24, 23/23 and 23/23 match, worst \|Δ\| 7.1e-15 |
| SEIR Np=100, same comparisons | 7 trees (one of 8 not built), 28 finite, 28/28 and 28/28, worst \|Δ\| 3.6e-15 |

In the isolation checks a match requires all four of: ll within atol 1e-8 / rtol 1e-10, identical state, identical
`cols.cols`, and an identical next `rand(UInt64)` from the global RNG.

The per-outcome checks compare against these closed forms:
- TCC/THH: `1 − C(ℓ,2)/C(n+1,2)`.
- THC/TCH no move: `1 − ℓ_dst/(n_dst+1)`.
- Cross: `ℓ_src·(1 − (ℓ_src−1)/n_src)/(n_dst+1)`. This is Φ_cr divided by q = 1/ℓ.
- π against `no_move_share`, to 1e-14.

The full weights `e^w/π` of the two slots are equal and equal `T0 + T1`. For BIRTH splits at ℓ = n, π₀ > 0 and w₀
is finite. For SEIR progression at ℓ_E = E, the test asserts π₀ = 0: the moving host is tracked, so the lineage must
move. This matches the reference and is not the upstream no-move defect.

Development runs, scratch only and not in the test file:
- Isolation, 3000 draws per model: 0 mismatches; worst \|Δll\| 7.1e-14 (MERS) and 5.7e-14 (SEIR).
- End-to-end, 17 MERS and 18 SEIR trees × 10 seeds, Np=100: 168/168 finite matches per model, worst 7.1e-15.

Existing tests rerun after the change:
- `kli_filter_ir_test.jl` + `kli_singular_update_test.jl`: 4381/4381.
- `kli_mers_compiled_test.jl` + `kli_seir_compiled_test.jl`: 58/58.

The full `test/runtests.jl` did not finish within the 180 s cap, so it has not been run to completion.

## Why ll agrees to about 1e-13 and not bit for bit
Three computations round differently from the compiled code:
- π of a BIRTH split comes from Rational `full_transitions` values. The compiled code uses the Float64 expression in `no_move_share`.
- The DEATH driver is `γI·(1 − ℓ/I)`. The compiled code uses `γ(I − ℓ)`.
- The DEATH decay leftover rounds the same way.

These change `rcateg` or the `t + step < tf` comparison only when a draw lands within an ulp of a boundary.
No run produced a different event sequence. The −Inf particles in the end-to-end runs show no RNG-consumption quirk:
every comparison counted had at least one finite ll, and all of them matched.

## Not covered
- SIR, MTBD, SI2R: their regular events fit the slot rule (`fork(a => a, b)`, `swap`, `chop`), but
  `regular_step!` was not run on them and there is no compiled reference to compare against.
- Proposals other than naive: `kli_decay` takes π as given, but only naive π exists here.
- BIRTH events that do not continue the parent, or that produce into two other demes: not supported (ArgumentError).
- The full test suite (see above).

## Follow-up (2026-10-07)
- Opus review (independent, against Aaron's patched `choose_move` functions): 40,580 per-outcome keys, 0 mismatches, worst difference 3.6e-15. Default paths unchanged.
- Added `_check_demes(cols, model)`: `singular_update!`, `regular_step!`, `kli_select` and `apply_move!` throw `ArgumentError` when the coloring's demeset has a different number of demes than `model.demes`. Only the count can be checked; the generic code still assumes the demeset lists its demes in the order of `model.demes`. Tested in `test/kli_regular_step_test.jl` (6 tests).
- Full suite, `RUN_HEAVY_TESTS=no julia --project=test test/runtests.jl`, run before the deme-count check was added: 45,814 passed, 0 failed, 3m04s. After the change only the five files that call the changed functions were rerun (`kli_regular_step_test`, `kli_singular_update_test`, `kli_filter_ir_test`, `kli_mers_compiled_test`, `kli_seir_compiled_test`): all pass.
