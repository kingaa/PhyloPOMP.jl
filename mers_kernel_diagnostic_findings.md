# MERS kernel ESS/collapse diagnostic — findings

**Date:** 2026-07-31 (see 2026-08-04 update below)
**Scope:** Follow-up to the mers_private / yang-phylopomp investigation. Goal: check
whether PhyloPOMP.jl's new support-safe kernels (Soft/Guided/Hard) measurably reduce
the particle-filter bottlenecks that `mers_private/ANALYSIS_GUIDE.md` §2.4 documents
for the R bootstrap filter on the real 274-tip MERS tree, since that repo's own notes
flag guided SMC as *"conceptual, not implemented."*

## Update (2026-08-04) — a second, simpler collapse mode under the module's bare defaults

Everything below this point describes one specific parameter regime: **Table 2 of
`Yang-Notes_V2.pdf`, passed explicitly** to `filter_pomp(...)` (in particular
`χ_H > 0`, since the diagnostic explicitly found zero Human-sample failures — that
requires a nonzero human sampling rate). That regime's conclusion (camel-population
extinction, `β_CH=0`/`B_C=0` absorbing state) still stands **for those parameter
values** and is unchanged by this update.

Since then, `mers_naive.jl`'s own commit history moved on ("change default
MERS parameters"), and the module's *bare* defaults — i.e. what you get calling
`NaiveMERS.filter_pomp()` with no overrides — now include **`chi_h = 0.0`**
(zero human sampling rate). Checked directly: the first `-Inf` under these bare
defaults occurs at **node 11 of 548 — the very first Human Sample node in the
tree** — not node 144, and not by extinction:

```julia
first human sample node index: 11 out of 548
first -Inf conditional logLik at step: 11 out of 548
```

The mechanism is a one-line, unconditional identity: at a destructive Human
sample, `ll += log(χ_H · I_h)`; with `χ_H = 0` this is `log(0) = -Inf` for
*any* `I_h > 0`, before any population dynamics — camel extinction, R0, tied
timestamps, none of it gets a chance to matter, because the very first human tip
in temporal order already forces `-Inf`.

**Practical upshot:** if you call `filter_pomp()` bare (no explicit `chi_h`), you
will see `-Inf` at node 11 regardless of which of the four kernels you use, and
this has nothing to do with the extinction analysis below — it's simpler and
fires first. The Table-2 extinction mechanism only becomes visible/relevant once
`chi_h` (and any other zero-by-default rate needed for the tree) is set to a
nonzero value that lets human sampling succeed, at which point camel sampling's
own extinction-driven `-Inf` (node 144, as documented below) becomes the
binding constraint instead.

## Summary

All four Julia kernels (Naive, Soft, Guided, Hard) collapse **identically** on the
full empirical tree at a set of parameters explicitly set to Table 2 of
`mers_private/Yang-Notes_V2.pdf` (i.e. not a disputed fitted point, and *not* the
`filter_pomp()` module's own bare defaults — see the 2026-08-04 update above for a
separate, simpler collapse mode under those): exactly
**100 of 548 filtering steps give `-Inf`**, the first at step 144, for every kernel,
every seed tested, and every parameter variation tested. Mapping all 100 back to the
genealogy: **every single one is a Camel-deme Sample node — zero Human, fork, or root
failures.**

**Definitive mechanism (confirmed by direct state inspection, not just node-type
mapping):** under these parameters, the simulated camel-infected count `I_c` goes
stochastically extinct (drops to exactly 0) early in essentially every trajectory —
confirmed in 15/15 independent free-simulation draws, all showing `I_c = 0` well
before node 144 — and **stays there permanently**, because `β_CH = 0` (no
human→camel transmission) and `B_C = 0` (no camel births) mean nothing can ever
regenerate a camel infection once `I_c` hits zero. Camel transmission
(`β_CC·S_c·I_c/N_c`) is itself proportional to `I_c`, so extinction is a true
absorbing state. Every camel sample from that point onward fails for a direct,
mechanical reason: the destructive-sampling likelihood contribution is
`log(χ_C · I_c)`, which is `log(0) = -Inf` exactly whenever `I_c = 0`. This is not a
proposal-kernel deficiency (every kernel proposes from the same population process)
and not really about tied timestamps either — see below.

## What was ruled out

Tested systematically (script transcripts not preserved, but each swept independently
holding everything else at the Table-2 defaults, `Np=1000`, `seed!(2121916527)` unless
noted):

1. **Kernel choice.** Naive, Soft, Guided, Hard all give exactly 100/548 collapsed
   steps, first at step 144. This alone is a clean, positive confirmation of the
   design: the epsilon-floor mechanism (`mers_funs.jl`) guarantees `π > 0` wherever
   the *population rate* is nonzero, but it cannot and must not manufacture proposal
   mass where the *target* itself is incompatible (`Φ = 0`) — exactly the KLI
   Theorem 5 boundary. A structural `Φ=0` bottleneck is invisible to every kernel
   equally, which is what we see.
2. **Cross-transmission rate** (`β_HC`, camel→human): swept 1, 2, 5, 10 — collapse
   count stayed in the 66–112 range, never near zero, and the *first* collapse at
   step 144 didn't move.
3. **Human-to-human rate** (`β_HH`): swept 3.65 (Table-2 default) through 70
   (near the high end anything in the yang-phylopomp/mers_private profiles
   considered) — no material change.
4. **Camel-to-camel rate** (`β_CC`, i.e. camel R0): swept 1×–5× the default
   (R0 ≈ 1.0 to 5.0) — **exactly** 100/548 collapsed every time. At the time this
   was read as ruling out "camel epidemic stochastically fades out"; **that reading
   was wrong.** Extinction is confirmed as the actual mechanism (see below) — the
   sweep didn't rule it out, it just showed that raising `β_CC` up to 5× doesn't
   prevent it: with `I_C0=50` and 100 destructive camel-sample removals draining the
   same population throughout the tree with no replenishment credit, stochastic
   extinction still occurs before all 100 samples are reached, at every tested R0 in
   this range, and once it occurs it's permanent (`β_CH=0`, `B_C=0`, no way back).
5. **Initial infected camels** (`I_C0`): swept 50, 200, 500, 1500 — **exactly**
   100/548 every time. Consistent with extinction being driven by the *sampling
   drain* (100 fixed removals over the tree) at least as much as by the starting
   population size — a larger `I_C0` delays extinction but doesn't prevent it before
   all 100 camel samples are needed.
6. **RNG seed**: swept 1, 2, 42, 12345, 2121916527, for both Naive and Soft —
   **exactly** 100/548, first@144, every time. This is the key clue: a real
   proposal-luck/Monte-Carlo problem would vary with the seed. Complete
   seed-invariance across a stochastic particle filter means the cause is
   deterministic, not stochastic.

## What was found

**Step 1 — map collapses to node type/deme.** Mapped every one of the 100 collapsed
steps (`cond_logLik(pf) == -Inf`, Naive kernel, `Np=1000`) back to its genealogy node
type and deme, for all 274 sample tips:

```
by node type:             {"Sample": 100}     — zero Node/fork failures, zero Root failures
Sample collapses by deme: {"Camel": 100}      — zero Human sample failures
same 100 node indices fail identically in Naive, Soft, Guided, and Hard (set-equal)
```

Every camel sample in the tree collapses; no human sample ever does.

**Step 2 — a tied-timestamp cluster was found and initially treated as the cause.**
Node 144 (`slate ≈ 3.7938`, the first collapse) sits right before a run of camel
samples with times identical to the 10th–12th significant digit (nodes 145–151, all
`≈ 3.8129357510...`) — a classic artifact of month-resolution collection dates
collapsing onto a shared numeric timescale. Checking all 100 collapses against the
full set of tree timestamps: 85/100 have another node within `1e-6` of the same
time; 15/100 have an ordinary, unique timestamp and still collapse (nodes 144, 154,
168, 170, 175, 240, 249, 330, 351, 385, 405, 477, 482, 485, 497). That 15/100 was the
signal that tied timestamps were not the whole story.

**Step 3 — direct state inspection resolved it.** `simulate(p, nsim=...)` returns a
full per-node state trajectory (`.states`, one entry per genealogy node, including
`I_c`/`ellC`/`cols`), so the actual state right before node 144 can be read directly
rather than inferred:

```julia
sm = simulate(p, nsim=15)
# every one of 15 independent draws:
state_at_node_143 → I_c = 0, ellC = 0 (cols[Camel] empty)
state_at_node_144 → ll = -Inf   (a Camel Sample node)
```

**`I_c` (simulated camel-infected count) is already exactly zero before node 144, in
100% of trajectories checked — not a coloring-tracking artifact, an actual
population-extinction event.** The proximate cause is mechanical: at a destructive
Camel sample, the likelihood contribution is `ll += log(χ_C · I_c)`; with `I_c = 0`
this is `log(0) = -Inf`, unconditionally, regardless of coloring or proposal kernel.

**Why extinction happens so reliably:** `β_CH = 0` (no human→camel transmission) and
`B_C = 0` (no camel births/immigration) mean `I_c = 0` is a true absorbing state —
camel transmission itself (`β_CC·S_c·I_c/N_c`) is proportional to `I_c`, so nothing
can ever restart it. With camel R0 `= β_CC/γ_C = 1.0` exactly (the critical
threshold — textbook near-certain eventual extinction for a critical branching
process) *and* a continuous drain from 100 destructive camel-sample removals against
a starting pool of only `I_C0 = 50`, extinction occurs early and permanently in
essentially every realization, well before all 100 camel samples in the tree are
reached.

**This also corrects the tied/isolated framing from Step 2**: it isn't two
mechanisms. Once extinction occurs (early, in every trajectory), *every subsequent
camel sample fails for the identical reason*, whether its timestamp happens to be
tied to another node's or not. Tied-timestamp clustering is a real, separate
phenomenon in this tree (still visible in the artifact below), but it is not causally
responsible for the collapse — it's a correlate, because camel samples are simply
scattered through the back two-thirds of the tree and post-extinction, all of them
fail regardless of clustering.

**Either way, this is a data/model-regime characteristic, not a bug in the new
kernels, and not something a better proposal kernel can fix.** A camel population
that has gone stochastically (and, under these parameters, permanently) extinct
presents a genuine `Φ=0` target incompatibility to every subsequent camel sample;
the epsilon-floor design correctly does not (and should not) manufacture proposal
mass to paper over that, which is exactly why all four kernels fail identically.

A visualization of this — the full tree, the tied-timestamp cluster exploded, the
per-kernel comparison, and a two-lane timeline mapping every one of the 274 samples
to collapsed/not-collapsed/tied/isolated — is published at
https://claude.ai/code/artifact/769d1032-e69f-4e28-920b-0bb58e345406.

## What this means for the two source projects

- **`mers_private`/`yang-phylopomp`'s own ESS diagnostics** describe a superficially
  similar phenomenon in the R bootstrap filter (~3% of steps catastrophic, invariant
  to `Np`) — this Julia-side investigation independently reproduces and, this time,
  fully localizes the mechanism, which those R notes did not pin down ("we do not
  know what type of genealogical event occurs at those steps" — `RESULTS.md` §4,
  limitation 2). Whether the R-side collapse has the *same* root cause (camel
  extinction under a similarly critical/absorbing parameter regime) or a different
  one hasn't been checked — worth a targeted follow-up given how model-regime
  specific this turned out to be, rather than assuming it transfers.
- **Implication for future fitting work (either R or Julia):** since the cause is
  population extinction under an absorbing-state parameter regime (`β_CH=0`,
  `B_C=0`, camel R0=1.0), not a tied-event-time artifact, the fix is **parameter
  choice, not event-time jitter or simultaneous-event handling** — e.g. allowing
  `β_CH>0` (even a small human→camel rate breaks the absorbing state) or fitting
  `β_CC` high enough, and with enough starting population, that camel extinction
  before all 100 samples is no longer near-certain. mif2-style fitting would need to
  avoid wandering into this absorbing region during optimization, which itself is a
  nontrivial parameter-search constraint worth flagging for any future fitting
  attempt (R or Julia).
- **Soft/Guided/Hard are not disproven here** — they were validated on their intended
  target (proposal-support gaps at otherwise-feasible outcomes) via the small
  hand-built fixture in `test/mers_soft.jl`/`test/mers_guided.jl`/`test/mers_hard.jl`
  and the TeX documentation's worked math. This diagnostic instead answers a
  *different*, previously-open question (does better proposal design help on the real
  tree's specific bottlenecks?) with a precise "no, because this specific class of
  bottleneck isn't a proposal problem" — which is itself a useful, falsifiable result.

## Suggested next steps (not started)

Root cause is now confirmed (camel-population extinction under an absorbing-state
parameter regime), so the open work is about *fixing* it rather than further
diagnosis:

1. **If finite likelihood on the full tree is wanted, change the parameter regime,
   not the tree encoding.** Tied timestamps and simultaneous-event handling are no
   longer believed to be the blocker (see above) — the fix is a `β_CH>0` (breaks the
   absorbing state) and/or a `β_CC`/`I_C0` combination where camel extinction before
   all 100 samples is no longer near-certain. Re-check the tied-timestamp cluster's
   effect only *after* extinction is no longer occurring — it may still cause
   smaller, more ordinary bottlenecks worth handling with jitter, but that's now a
   secondary question, not the primary blocker.
2. **Re-run the R-side ESS diagnostic with this hypothesis in mind.** If
   `mers_private`'s bootstrap-filter collapse is caused by the same camel-extinction
   mechanism (plausible, since their default parameters are the same Table 2 this
   diagnostic used), their "which event type" gap (`RESULTS.md` §4, limitation 2) has
   the same answer, and their next step is also parameter choice, not proposal
   design or event-time handling.
3. **Consider whether a non-destructive-sampling or waning-immunity variant of the
   model is more appropriate for camels** — Yang's original SIRS formulation (with
   `ω_1>0`, per `Yang-Notes_V2.pdf` Table 2) explicitly exists to prevent exactly this
   kind of endemic-reservoir extinction over a multi-year tree; this repo's model
   currently has no waning/immigration term for either deme.
