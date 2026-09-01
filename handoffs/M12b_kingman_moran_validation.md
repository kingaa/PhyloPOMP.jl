# Handoff — M12b: Kingman Coalescent / Moran Special-Case Validation

## Status
COMPLETE

## Objective
Inserted between M12 and M13 specifically to answer the project's human
user's question: "convince me this actually works and is not circular
logic." Every prior milestone (M00-M12) validated the compiler against
`seir_naive.jl`/`mers_naive.jl` (hand-coded filters) and
`mers_filter_suite.tex` (a hand derivation) — all three project-authored,
per M11's own audit (`handoffs/M11_model_crossvalidation.md`). This
milestone instead checks the compiler against a source external to the
project entirely: the standard Kingman coalescent / Moran-model
combinatorial result (Kingman 1982; textbook statement, e.g. Wakeley,
*Coalescent Theory*, §3), the first of `mgp_filter.jl`'s "Validation
gates" (§4.3.1, Eqs. 17-20 of the KLI paper). The linear birth-death
(Stadler) special case, §4.3.2, is explicitly left for a follow-up
milestone — not attempted here.

## Which event to check, and why (a judgment call this milestone had to
resolve, not a given)
The task brief's own working hypothesis was correct on inspection: a true
Moran-style "one step, both roles land in the same finite pool" event
needs BOTH of its production slots in the SAME deme. SEIR's `infection`
(r=(1,1), demes (E,I)) and MERS's `transmission_hc`/`transmission_ch`
(r=(1,1), one slot in each of two DIFFERENT demes) do not have this shape
— a "coalescence" there would need one slot drawn from one pool (E or H)
and the other from a different pool (I or C), so there is no single N to
plug into `C(N,2)`. MERS's `transmission_cc`/`transmission_hh` (TCC/THH,
r=(2,0)/(0,2), `move=fork(I_c=>I_c,I_c)`/`fork(I_h=>I_h,I_h)`) are the
only two BIRTH marks in either shipped model whose production vector puts
both slots in one deme — confirmed directly (not assumed) by inspecting
`src/examples/mgp_mers.jl`'s `@event` declarations, and independently
reconfirmed by `full_transitions`' own `ForkTransition.slot_demes`, which
comes out `[1,1]` (TCC) / `[2,2]` (THH) — both entries equal, i.e.
genuinely "same pool" — versus THC/TCH's Fork case, which is never even
reachable at the TWO-slots-tracked level in the same way (THC/TCH each
have only ONE production slot, r=(1,1) means one slot per event, not two
slots in one deme; there is no `s=(2,0)`-style saturation for them at
all). TCC/THH are therefore the correct and only candidates in this
codebase for the Kingman/Moran check.

## The standard result, derived from scratch (not read off this project's
tex)
**Setup.** Moran model, population size `N`. A single birth-death step
picks one individual to reproduce (parent role) and one individual/slot to
be replaced (offspring role). Tracing genealogies backward, this step
causes two specific tracked lineages `i,j` to coalesce exactly when BOTH
roles land on `{i,j}` — i.e., the step's causally-involved pair of
population members is exactly `{i,j}`.

**Derivation.** Treat the `N` individuals as exchangeable (the standard
neutral-model symmetry assumption — no individual is combinatorially
distinguished). Then each of the `C(N,2)` unordered pairs of individuals
is equally likely to be "the pair causally involved" in a given step, so:

    P(a SPECIFIC tracked pair {i,j} is the one realized this step) = 1/C(N,2)

Because at most one pair can be "the" causally-involved pair per step (the
events "pair {i,j} is involved" for different pairs are mutually
exclusive), summing this probability over the `C(ell,2)` pairs drawable
from `ell` tracked lineages is valid (not an overcount), giving:

    P(ANY of the C(ell,2) tracked pairs coalesces this step) = C(ell,2)/C(N,2)

This is an **exact finite-N identity**, not a diffusion/large-N
approximation — the familiar "coalescent time units, rate `C(ell,2)`"
statement of Kingman's theorem is the `N -> infinity` limit of exactly
this formula, not a separate fact requiring separate justification.

This combinatorial fact is verified in `test/kli_kingman_moran_test.jl`'s
"Part 0" testset by brute-force double-loop enumeration (no `binomial()`
call at all) for `N` up to 14 and every `ell` in `0:N` — 390 assertions —
BEFORE `Base.binomial` itself is trusted for the larger randomized sweeps
that follow.

## Mapping to this project's compiled code

`kli_binomial_ratio(r, s, ell, n; Q=1)` (M02, `mgp_phi.jl`) computes

    φ_u(s) = Q · ∏_d C(n_d - ell_d, r_d - s_d) / C(n_d, r_d)

For TCC (`r = (2,0)`, both slots in deme C) at the full saturation
`s = (2,0)` (both slots filled by tracked lineages — a `ForkTransition`,
M03):

    φ_fork = C(I_C - ell_C, 0)/C(I_C, 2) · C(I_H - ell_H, 0)/C(I_H, 0)
           = 1/C(I_C, 2) · 1/1
           = 1 / C(I_C, 2)

— independent of `ell_C`, `ell_H`, and `I_H` entirely (the `C(a,0)=1`
identity kills every factor except the camel-deme denominator). This is
**exactly** the "this SPECIFIC pair is the one realized" probability
derived above, with `N = I_C` (the post-event camel-deme population). The
`C(n_d-ell_d, r_d-s_d)` numerator being fixed at `C(·,0)=1` is precisely
the algebraic reason `φ_u`, as M02/M03 define it, is "which SPECIFIC pair"
rather than "which OF THE ell tracked lineages, i.e. any pair" — the
`C(ell_d, s_d)` "which specific tracked lineages fill the slots"
multiplicity is deliberately NOT baked into `φ_u` at all (verified
directly by evaluating `kli_binomial_ratio` at many `(I_C,ell_C)`
instances and confirming the ratio depends only on `I_C`, never `ell_C`
— not assumed from `mers_filter_suite.tex`'s Step C comment, which makes
the identical claim in prose and is cited only as a secondary,
after-the-fact cross-check, matching the task's explicit instruction).

To get "ANY of the `C(ell_C,2)` tracked pairs coalesces," the same
disjoint-summation argument used above applies: multiply by
`C(ell_C,2)`. Combined with the event's REAL compiled hazard,
`alpha_TCC(x,θ) = tcc.hazard(x,θ) = β_cc·S_c·I_c/N_c` (M01's `Event`,
unmodified), a Poisson process of rate `alpha_TCC`, each occurrence
independently causing "any tracked pair coalesces" with probability
`C(ell_C,2)·φ_fork`, is — by Poisson thinning — exactly a Poisson process
of "tracked-pair-coalescence" events at rate

    alpha_TCC(x,θ) · C(ell_C,2) · φ_fork = alpha_TCC(x,θ) · C(ell_C,2)/C(I_C,2)

which is **exactly** the standard finite-N Moran coalescence rate for a
population of size `N = I_C` undergoing reproduction events at rate
`alpha_TCC`, with `alpha_TCC` playing the role of "the rate at which a
Moran birth-death step occurs." THH mirrors this with `H` in place of `C`.

As a secondary, after-the-fact cross-check only (not the derivation path):
`mers_filter_suite.tex` line 758 independently states the filter's
"no-fork" rate contribution as
`alpha_TCC(t,x',x)[1 - ell_C(ell_C-1)/(I_C(I_C-1))]`, i.e. its own
implicit "fork/any-pair" rate is `alpha_TCC · ell_C(ell_C-1)/(I_C(I_C-1))
= alpha_TCC · C(ell_C,2)/C(I_C,2)` after the standard `C(k,2)=k(k-1)/2`
identity — algebraically identical to the formula derived above from
scratch. This agreement is reported as a bonus, not used as the source of
the target formula (which was derived independently, per the task's
explicit instruction, before this cross-check was even looked up).

## Numerical verification

`test/kli_kingman_moran_test.jl`, registered in `test/runtests.jl`
(one line added, after `kli_properties_test.jl`). 2507 individual
`@test` assertions total, all passing, in exact `Rational{Int}`
arithmetic throughout (no `Float64` anywhere in this file):

| Part | What it checks | Trials | Calls (real project code) |
|---|---|---|---|
| 0 | Brute-force ground truth for `C(N,2)`, `C(ell,2)/C(N,2)`, `N` in `0:14`, every `ell` | 390 asserts | none (pure brute force + `Base.binomial` cross-check) |
| 1 | `kli_binomial_ratio`'s TCC/THH Fork saturation `== 1/C(N,2)`, and independence from the other deme's state | 200+200 trials, 808 asserts | `kli_binomial_ratio` (M02) |
| 2 | `full_transitions`' `ForkTransition` reproduces the same `phi`, `slot_demes` confirms "same deme" shape, Chu-Vandermonde sum-to-1 | 100 trials, 400 asserts | `full_transitions` (M03) |
| 3 (**main claim**) | `alpha_TCC/THH · C(ell,2) · φ_fork == alpha · moran_pair_prob(N,ell)`, target computed independently | 150+150 trials, 902 asserts | `Event.hazard` (M01) + `kli_binomial_ratio` (M02) |
| 4 | Boundary states (`N=0,1`; `ell=N=2`) behave sanely | 7 asserts | `kli_binomial_ratio` |

`I_C`/`I_H` were drawn up to 500-ish (`rand(rng, 0:400)` / `2:500`
depending on the testset); `ell` uniformly from `0:N`; hazard parameters
(`S_c`, `N_c`, `β_cc` etc.) drawn as strictly-positive random
`Rational{Int}` values so `alpha_TCC/THH != 0` and the check is not
vacuous. Every comparison is **exact** `Rational{Int}` equality
(`==`, not `≈`); Part 3 additionally tracks `worst_diff` across all 150
trials per event and asserts it is exactly `0//1` — there is no
floating-point tolerance anywhere in this file, and no discrepancy was
found at any of the (390+808+400+902+7 =) 2507 checked instances.

Full suite: `RUN_HEAVY_TESTS=no julia --project=. test/runtests.jl` was
run before (confirmed baseline `18759/18759`, matching M12's stated
count) and after (`21266/21266` = `18759 + 2507`, i.e. every new
assertion passes and nothing regressed).

## Scope: precisely what this validates, and what it does NOT

This validates that `kli_binomial_ratio`/`full_transitions`'
`ForkTransition` machinery, evaluated at MERS's two within-deme BIRTH
marks (TCC, THH), reduces EXACTLY to the standard finite-N Moran
coalescence-rate formula, **at a fixed, frozen instantaneous state**
`(I_C, ell_C)` (resp. `(I_H, ell_H)`) — an exact algebraic identity
between two Rational numbers computed from the same `(N, ell)`, holding
for every reachable `N >= 0`, `0 <= ell <= N` tested (and, by the
algebraic derivation above, for every such state, not just the sampled
ones). It does NOT hold only in a large-`N` or diffusion limit — it is
exact at every finite state, which is a stronger and more useful fact
than the textbook's usual continuum statement, precisely because the
finite-N Moran model IS the exact, non-limiting object Kingman's theorem
is a limit of.

What this does **not** validate:
- The driver/decay/reduction machinery (`mgp_decay.jl`, `mgp_reduce.jl`,
  `mgp_proposal.jl`) is untouched by this check — TCC/THH's `Identity`
  and `InlineSameDemeTransition` branches (which those pieces DO consume,
  per `mgp_mers_filter.jl`'s Finding 2) are exercised only as a
  Chu-Vandermonde structural sanity check (Part 2), not independently
  validated against an external source here.
- The other four MERS marks (THC/TCH, cross-deme) and both SEIR marks are
  explicitly OUT of scope — they are not Moran-shaped (see "Which event"
  above) and are not addressed by this milestone at all.
- The linear birth-death (Stadler) special case, §4.3.2, is not attempted.
- This says nothing about whether the SINGULAR/branch-point wiring in
  `mers_naive.jl`'s `singular_part!` (the only place `ForkTransition`'s
  `phi` is actually consumed at runtime, per `mgp_mers_filter.jl`'s
  Finding 2: "`ForkTransition` ... is NEVER realized regularly by ANY of
  the four MERS BIRTH marks") correctly implements the KLI likelihood
  contribution for an OBSERVED branch point — this milestone only checks
  the underlying `φ_u` value against the Moran target, not the wiring
  that consumes it. (As a non-primary observation: `mers_naive.jl`'s own
  `singular_part!`, line 77, computes `log(lambda_cc) - log(Ic*(Ic-1)/2)`
  for an observed TCC branch point — literally `alpha_TCC / C(I_C,2)`,
  the "specific observed pair" density this milestone derives — which is
  a third, independent piece of code arriving at the same number, but
  this was noticed in passing while reading for context, not used as a
  derivation source, and is not itself re-verified here.)
- No claim is made about parameter regimes, boundary/degenerate dynamics
  beyond the handful directly tested (Part 4), or about the correctness
  of any event OTHER than TCC/THH.

## Does this escape the "project-authored circularity" concern?

**Partially, and for a specific, well-defined slice of the compiler —
yes, genuinely — but it does not (and cannot, by itself) close the
concern for the whole system.**

Why it genuinely helps: the TARGET formula in this milestone,
`C(ell,2)/C(N,2)`, is derived here from elementary, textbook, pre-project
probability (a symmetric-population birth-death step, exchangeability,
disjoint-event additivity) and cross-checked by brute-force enumeration
that calls none of this project's code. It is then compared, via exact
rational arithmetic, against the REAL `kli_binomial_ratio`/
`full_transitions`/`Event.hazard` functions — not a reimplementation, and
not a re-statement of `mers_filter_suite.tex`. Agreement here is evidence
of a different KIND than M02-M11's agreement: it rules out the
specific failure mode M11 flagged (the whole project inventing a
self-consistent but wrong theory) for the one piece checked, because the
Moran formula did not come from inside the project at all.

Why it does not close the concern in general: (1) it covers exactly one
structural pattern — the "two production slots, one deme" Fork case, at
exactly two of the twelve MERS event marks (TCC/THH) — not the driver,
decay, reduction, cross-deme, or SEIR machinery, all of which remain
validated only by the project's own internal (M02-M11) cross-checks;
(2) the KLI framework's actual novel content is precisely the STRUCTURED
part — multi-deme genealogies, migration, decay/boost corrections for
subsampling, the driver assembly (Eqs. 44-47) — none of which reduces to
a textbook special case, so no external check of comparable rigor is
available for those pieces (this is exactly why the KLI paper itself
only offers TWO special-case reductions, §4.3.1 and §4.3.2, as its own
sanity checks — the framework's generality is inherently
"un-special-case-able" beyond that); (3) a single external check, however
clean, is one data point — it would take an external check at every
structurally distinct piece of the compiler (which does not exist in the
literature, hence this whole project) to fully retire the concern. The
honest summary: this milestone converts "not circular" from an
unaddressed worry into a specific, falsifiable, EXTERNALLY-sourced claim
for the Fork/Moran slice of the compiler, tested at 2507 exact instances
with zero discrepancies — genuine, non-trivial evidence, but explicitly
scoped evidence, not a blanket proof that the whole compiler is correct.

## Things that gave me pause (reported honestly, not smoothed over)
- The `φ_fork` formula's complete independence from `ell_C` (and from the
  other deme entirely) initially looked suspicious — a probability that
  doesn't depend on how many lineages are tracked seems wrong at first
  glance. Working through the hypergeometric identity resolved this: `φ_u`
  is defined (M02) as "probability THIS EXACT identified set of tracked
  lineages fills the slots," not "probability SOME set of size `s` does,"
  so of course it's `ell`-independent once you fix WHICH `ell_C`-sized
  identified subset you're asking about — the `ell`-dependence lives
  entirely in the `C(ell_C,2)` multiplier applied afterward. This matches
  `mgp_mers_filter.jl`'s own Finding 2 language exactly, which is
  reassuring but was independently re-derived here before that section
  was consulted.
- `ForkTransition`'s `phi` is, per `mgp_mers_filter.jl`'s Finding 2,
  "NEVER realized regularly by ANY of the four MERS BIRTH marks" — it is
  the SINGULAR-time code path (`mers_naive.jl`'s `singular_part!`,
  untouched/read-only) that actually consumes this number at runtime.
  This milestone checks the `φ_u` VALUE the compiled M02/M03 machinery
  produces, calling those functions directly, rather than checking the
  currently-wired filter's runtime behavior (which never calls
  `full_transitions` on TCC/THH's Fork case at all, by Finding 2) — this
  is the correct scope per the task brief ("calling the REAL M02/M03
  functions", not "checking the currently-wired filter"), but it is worth
  being explicit that this is a check of the COMPILER'S MATH, not of an
  end-to-end filter execution path, since TCC/THH's Fork case is not
  currently wired into any executable filter's regular-time loop.

## Files added
- `test/kli_kingman_moran_test.jl`
- `handoffs/M12b_kingman_moran_validation.md` (this file)

## Files modified
- `test/runtests.jl` — one line (`include("kli_kingman_moran_test.jl")`,
  after `kli_properties_test.jl`).

No existing source file (`src/**`) was modified. `mers_naive.jl`,
`seir_naive.jl`, and the other hand-coded filter modules were read for
context (see "things that gave me pause") but not modified, per the
task's constraints.
