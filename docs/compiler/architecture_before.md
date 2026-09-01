# Architecture as it exists today (M00 reconnaissance)

Branch: `atpabuser-devel`, HEAD at time of writing: `931af3b` ("Remove
epsilon-floor proposal machinery from MERS kernels"). This document is a
factual map of what is actually in the repository, not a design proposal.
File:line references are to that commit.

**Naming collision warning.** The repo root already has
`check_milestone1.md`, `check_milestone2.md`, `check_milestone3.md`,
documenting a *prior*, unrelated milestone sequence (forward simulator work:
M1 = generic Gillespie simulator, M2 = MERS simulator generalization, M3 =
Newick zero-length-branch edge case). The `M00`/`M01`/... numbering used by
*this* compiler project is a fresh, independent sequence and should not be
confused with those. `handoffs/M00_reconnaissance.md` (this milestone's own
handoff) is the first file in the new sequence's directory.

## 1. Module structure / entry point

`src/PhyloPOMP.jl` is the package root. It declares three package-wide type
aliases (`Name=UInt64`, `Size=UInt64`, `Time=Float64`, `Prob=Float64`,
lines 10-13) and `include`s, in order: `demes.jl`, `genealogy.jl`,
`parse.jl` (Newick parser), `newick.jl` (Newick writer), `cblv.jl`,
`coloring.jl`, `fsmarkov.jl`, `guide.jl`, `rcateg.jl`, `indicator.jl`,
`examples/Examples.jl`, `simulate.jl` (lines 21-56).

`src/examples/Examples.jl` (10 lines) is a flat list of further `include`s:
`seir_{naive,soft,guided,hard}.jl`, `mers_{naive,soft,guided,hard}.jl`,
then `mgp.jl`, then `mgp_filter.jl`. Note `mgp.jl` and `mgp_filter.jl` are
loaded *last*, after all eight hand-coded model modules — the MGP compiler
scaffold is architecturally bolted on, not the thing everything else is
built from.

## A. Model DSL parsing

Two DSLs exist, both macro-based, in different files with different scope:

- **`@demes`** (`src/demes.jl:10-12`) — thin wrapper around `EnumX.@enumx`
  that builds a deme enumeration module. Used everywhere (SEIR/MERS hand
  filters and the MGP scaffold alike) to name the lineage-carrying demes.
- **`@mgp` / `@event`** (`src/examples/mgp_macro.jl`) — the actual
  "model DSL" for the compiler project. `@mgp Name begin ... end`
  (`mgp_macro.jl:144-186`) parses `compartments=`, `demes=`, `params=`
  declarations plus a sequence of `@event name rate=... pop=... move=...
  kind=...` calls, and expands to a `const Name = MGPModel(...)` binding
  plus an auto-generated `@demes` call for the lineage-carrying demes.
  `_decode_move` (`mgp_macro.jl:69-106`) interprets the small `move=`
  vocabulary — `fork(src => targets...)`, `swap(src => target)`,
  `chop(deme)`, `sample(deme)`, `sample_remove(deme)`, `none` — into an
  `EventType` plus `from`/`into`/production-vector `r` triple.
  `_rewrite` (`mgp_macro.jl:9-18`) does simple symbol substitution so bare
  compartment/parameter names in `rate=` expressions become `x.name` /
  `θ.name` field accesses; this is literally the DSL's only "compilation"
  step — no algebra, no IR beyond the resulting `Event`/`MGPModel` structs.
  There is no separate `parse_model` pass distinct from macro-expansion:
  parsing and lowering to `MGPModel` happen in one macroexpansion.
- Two concrete models are declared this way: `@mgp SEIR begin ... end`
  (`mgp.jl:188-198`) and `@mgp MERS begin ... end`
  (`mgp_mers.jl:3-24`, included from the end of `mgp.jl:120` /
  `mgp_macro.jl:199`).
- `mgp.jl:91-107` also hand-writes an *explicit* second copy of the SEIR
  event table, `SEIR_REFERENCE :: MGPModel`, purely as an audit oracle for
  the macro (see item N below) — not a second DSL, just a literal struct
  literal checked equal to the macro's output by
  `test/seir_macro_equivalence.jl`.

There is no parser for the eight hand-coded SEIR/MERS filter modules
(`seir_naive.jl` etc.) — those are ordinary Julia functions with the
population/hazard/coloring math typed out inline; they predate `@mgp` and
are not derived from it in any way (see item M).

## B. Event representation

`struct Event` (`mgp.jl:59-69`): `name::Symbol`, `Δ::Vector{Pair{Symbol,Int}}`
(state jump), `hazard::Function` (closure `(x,θ)->Float64`), `r::Vector{Int}`
(production vector, one entry per lineage-carrying deme), `type::EventType`
(one of `BIRTH, MIGRATION, DEATH, SAMPLE, NEUTRAL` — `mgp.jl:32`, the five
*pure* event types of KLI §3.2; compound events §3.2(f) are explicitly out
of scope, `mgp.jl:30`), `from::Int` (source deme index, 0 = none),
`into::Vector{Int}` (target deme indices), `regular::Bool` (regular vs.
singular), `observed::Bool` (can it produce an observed genealogy feature).

This is exactly the *static per-mark data* (Δ_u, α_u, r_u, W_u/from/into)
the task background describes. The docstring at `mgp.jl:52-57` explicitly
states saturation `s` and lineage count `ℓ` are *not* stored on `Event` —
they are runtime-genealogy-dependent and computed by the filter — and that
the specific coloring operator (χ/σ/κ) is *not* a stored field either,
because one mark is compatible with multiple coloring outcomes via the
compatibility indicator `Q_u`. So the static/dynamic split from the task
background is already explicitly recognized in the `Event` docstring, even
though the dynamic side (`kli_select`/`apply_move!`, see D/H/I/J below) is
unimplemented.

There is no explicit field for τ_u (event time) on `Event` — that is
naturally the simulation/filter clock `t`, not part of the static mark
description, consistent with `Event` describing an *event type*, not an
occurrence.

## C. MGPModel representation

`struct MGPModel` (`mgp.jl:77-82`): `name::Symbol`, `compartments::Vector{Symbol}`
(full state vector), `demes::Vector{Symbol}` (the lineage-carrying subset —
KLI §2.6's "coloring set"), `events::Vector{Event}`. That's the entire
model representation: no separate initial-condition spec, no time-varying
covariates, no compound-event support. `SEIR` (`mgp.jl:91-107` /
`mgp.jl:188-198`) and `MERS` (`mgp_mers.jl:3-24`) are its only two
instances.

## D. Hazards

Stored as an opaque Julia closure `Event.hazard :: Function`, built by
`@mgp`'s `_rewrite` (`mgp_macro.jl:9-18`, `133`) from the `rate=` expression
verbatim (arbitrary Julia syntax is accepted — `_rewrite` recurses through
any `Expr`, substituting only bare compartment/parameter symbols). Evaluated
mechanically by `kli_hazard(ev,x,θ) = Float64(ev.hazard(x,θ))`
(`mgp_filter.jl:36`) — this is the one function in `mgp_filter.jl` explicitly
marked as *not* a stub (mechanical helpers section, `mgp_filter.jl:25-27`).
There is no separate hazard IR or symbolic representation; hazards are
never introspected, only called.

## E. Δ (state change) representation

`Vector{Pair{Symbol,Int}}` on `Event.Δ` (`mgp.jl:61`), e.g.
`[:S=>-1, :E=>+1]`, built by `_pairs` (`mgp_macro.jl:20-38`) from the
`pop=(S=-1, E=+1)` `@event` keyword. Applied mechanically by two separate
(structurally identical) functions: `apply_pop` in `mgp_filter.jl:44-50`
(used by the filter scaffold) and `apply_delta` in `simulate.jl:79-88`
(used by the forward simulator) — same logic duplicated in two files, not
shared.

## F. Production vectors r

`Event.r :: Vector{Int}` (`mgp.jl:63`), indexed by position in
`MGPModel.demes`, built by `_production` (`mgp_macro.jl:61-67`) by counting
how many times each target deme appears in a `fork`/`sample` move's target
list. Per the docstring (`mgp.jl:43-44`) this is the KLI §3.3.1 per-mark
constant `r^u`. It is stored but, per grep, never actually *read* anywhere
in `mgp_filter.jl` or `simulate.jl` — no consumer of `Event.r` exists yet in
the scaffold (it's populated for future use, not wired into any
computation). The hand-coded SEIR/MERS filters do not use an `r` vector at
all; they have the equivalent information hand-inlined into each event's
rate-splitting logic (e.g. `seir_naive.jl:111-119`'s `alpha[1]=alpha[2]=...`
pairs for the two orientations of a birth event).

## G. Source/target deme wiring

`Event.from::Int` / `Event.into::Vector{Int}` (`mgp.jl:65-66`), deme
*indices* into `MGPModel.demes` (not the enum values themselves), set by
`_decode_move` (`mgp_macro.jl:69-106`) per move type: `fork` → `from` =
source deme, `into` = target demes (self excluded if present, e.g.
`fork(I => E, I)` gives `into=[E]` — see MERS's `fork(I_c => I_c, I_c)`
giving `into=[I_c]`, a self-loop kept as one element after de-duplication);
`swap` → single `into` entry; `chop`/`sample`/`sample_remove` → `from` only,
empty `into`; `none` → `from=0`, `into=[]`. This is the closest existing
analogue to the task background's "W_u" (source/target deme wiring) concept
— it exists, but only as plain integer indices, with no richer wiring
structure (e.g. no explicit orientation/labeling of which child continues
which parent lineage — that's left to be decided at runtime by
`apply_move!`, currently a stub).

## H. Genealogy / coloring representation

Two *separate* representations exist for two different purposes, and the
module header of `simulate.jl:15-21` explicitly documents why they are kept
separate rather than unified:

- **`Genealogy{D}`** (`src/genealogy.jl:53-67`) — the full, fixed, given
  tree (from data or from `simulate`): a flat time-ordered
  `Vector{GenealNode}` (`genealogy.jl:18-42`) with `type::NodeType`
  (`Root`/`Sample`/`Node`, `genealogy.jl:9`), `slate::Time`, `deme`,
  `lineage`, `parent`, `children`. `repair!` (`genealogy.jl:108-136`)
  re-sorts/renames/re-traces lineages after mutation.
- **`Coloring{D,N}`** (`src/coloring.jl:10-19`) — the *particle filter's*
  live state: `N` `BitSet`s (one per lineage-carrying deme), each holding
  the currently-tracked lineage IDs. Mutators `swap!`, `chop!`, `fork!`,
  `plant!` (`coloring.jl:36-111`) are the mechanical primitives a coloring
  operator would call; `ell(y)` (`coloring.jl:27-29`) returns per-deme
  lineage counts ℓ_d. This is a *reduced* coloring: it tracks only `d`
  (which deme each tracked lineage is currently colored), never an
  event-indicator/mark component `m`. Per the task background's Y=(d,m)
  distinction, **the full KLI coloring is not represented anywhere in this
  codebase** — everything here (hand-coded filters and the MGP scaffold
  alike) works directly with the reduced `d`-only `Coloring`, i.e. the
  marginalization over `m` is assumed to have already happened, not
  computed.
- **`SimInventory`** (`simulate.jl:54-57`) — a *third*, simpler
  bookkeeping structure, one `Vector{Name}` of live lineages per deme, used
  only by the forward simulator (tracks every extant individual, not just a
  filter-relevant sparse subset — a different growth regime from
  `Coloring`, per the module header's explicit rationale,
  `simulate.jl:15-21`).

## I. Coloring operators χ/σ/κ

`swap!`/`chop!`/`fork!`/`plant!` in `coloring.jl` are the *mechanical*
primitives (mutate a `Coloring` in place, return new ℓ). They are not
labeled χ/σ/κ in code, but the docstrings of the MGP scaffold's stub
(`apply_move!`, `mgp_filter.jl:92-118`) explicitly map: chop=χ,
swap=σ, fork=κ (`mgp.jl:24-31`'s `EventType` docstring; `mgp_filter.jl:127`
comment). What is genuinely missing is the *dynamic selection logic*: given
an event `ev` and the current `cols`, which of possibly-several compatible
coloring outcomes (i.e. the actual Q_u/φ_u computation) to apply, and with
what importance weight. That logic:
  - **Exists, hand-derived per event, in the eight SEIR/MERS filter
    modules** (`seir_naive.jl:102-178`'s `singular_part!` and
    `NaiveSEIR`'s `regular_part!` embedded rcateg calls — e.g.
    `seir_naive.jl:161` `rcateg([n.present[1,1]*n.present[2,2],
    n.present[1,2]*n.present[2,1]], true)` picking one of two fork
    orientations with hand-written weights).
  - **Does not exist generically in the MGP scaffold**: `apply_move!`
    (`mgp_filter.jl:116-118`) is a bare `error(...)` stub. No φ_u/Φ_u
    computation, no `Q_u` compatibility check, no saturation enumeration
    exists in `mgp_filter.jl` at all.

## J. Proposal logic (π_u)

Exists richly, but only hand-derived per model, never as a generic
π_u-from-model-data computation:

- **SEIR**: four hand-written proposal families, each a separate module —
  `NaiveSEIR` (`seir_naive.jl`, no guide — uses raw lineage-count ratios,
  e.g. `pi[1]=@indicator(I>0,1-ellI/I)` at `seir_naive.jl:115`),
  `SoftSEIR` (`seir_soft.jl`), `GuidedSEIR` (`seir_guided.jl`), `HardSEIR`
  (`seir_hard.jl`) — Soft/Guided/Hard build an anticipatory π from a
  reverse-time guide sweep (`src/guide.jl`'s `Guide`/`relhaz`/
  `choose_branch`), the Soft/Guided/Hard split matching KLI's "reverse-time
  finite-state Markov guide" discussion referenced in `mgp_filter.jl:153-156`.
- **MERS**: the same four-way split — `mers_naive.jl`, `mers_soft.jl`,
  `mers_guided.jl`, `mers_hard.jl` — sharing common code
  (`knowledge!`, root/fork `singular_part!`) via `mers_funs.jl`, generalized
  to MERS's two-host cross-deme structure.
  `mers_filter_suite.tex` §"From Math to Code" (lines 1153-1174) explicitly
  documents each MERS kernel's SEIR analogue.
  A prior version of these kernels (see commit history, `handoff.md`
  "Follow-up round 4" item 4) added an ε-floor mixing mechanism
  (`floored_shares`, `floored_branch_law`, `hard_branch_law`,
  `cross_proposal`, a `proposal_floor` parameter) as a *proposal-only*
  robustness hack. **Commit `931af3b` (current HEAD) removed all of it** —
  confirmed by `grep -rln "floored\|proposal_floor\|epsilon" src/ test/`
  returning only LaTeX build artifacts (`.toc`/`.aux`/`.log`) and one stale
  comment in `test/mers_naive.jl:17` ("and no `proposal_floor`") describing
  what does *not* exist. `mers_filter_suite.tex:1144-1151` has a "v31
  update" notebox recording the removal and stating MERS now matches SEIR's
  unfloored pattern. **No target-math epsilon-floor hack exists in the
  current tree** — the task background's specific concern (proposal-epsilon
  baked into the *target* math) does not apply to HEAD; the floor was
  always in π (the importance kernel), never in Φ, and has since been
  deleted outright.
- **MGP scaffold**: `kli_select` (`mgp_filter.jl:56-72`) is a bare
  `error(...)` stub — no generic π computation exists there at all; the
  docstring states it "for SEIR ... contains the lineage-count factors
  represented by `pi` in `seir_naive.jl`" but this is aspirational
  commentary, not implemented code.

## K. Driver/boost/decay logic

KLI's β_u = α_u·π_u (driver) and B_u = Φ_u/π_u (boost) terminology appears
explicitly only in `mgp_filter.jl` docstrings (`kli_select`'s docstring at
line 61, `apply_move!`'s at line 105) as forward references to what the
stubs must eventually compute — not as implemented functions. `kli_decay`
(`mgp_filter.jl:74-89`) is likewise a stub citing "SEIRS decay KLI Eq. 47"
and Appendix B Eq. B2. In the hand-coded filters, the *decay* accumulator
is a real, working, per-model return value: `NaiveSEIR.event_rates!`
(`seir_naive.jl:103-121`) returns a `decay::Prob` folded into
`ll -= decay*step` in `regular_part!` (`seir_naive.jl:145,179`); the MERS
analogue is `NaiveMERS.event_rates!` (`mers_funs.jl` equivalents /
`mers_naive.jl:132-162`). So "driver/boost/decay" as a *concept* is fully
realized per-model by hand, and named/cited but not implemented generically.

## L. Filter execution loops

Two independent loop implementations:

- **Hand-coded, per model** (the thing actually used/tested): each of
  the eight SEIR/MERS `{naive,soft,guided,hard}.jl` modules builds a
  `POMP.pomp(...)` object (`PartiallyObservedMarkovProcesses` package,
  reexported via `PhyloPOMP.jl:19`) with `rinit`, `rprocess = onestep(...)`,
  `logdmeasure` closures. `rprocess`'s closure calls the model's
  `singular_part!` then `regular_part!` per genealogy interval
  (`seir_naive.jl:223-246` is representative). Advancement to the next
  regular event is via `rcateg` (`src/rcateg.jl`, categorical draw over
  unnormalized weights) inside a `while t<tf` loop drawing exponential
  waiting times — same Gillespie-style direct-method pattern used by
  `simulate.jl`. This is the actual, tested, working filter path (SMC via
  `POMP.pfilter`, exercised in `test/*_naive.jl` etc. and
  `mers_filter_suite.tex`'s benchmark table).
- **Generic MGP scaffold**: `regular_step!` (`mgp_filter.jl:163-188`) — a
  single, model-agnostic function implementing the same
  hazard/decay/exponential-step/`rcateg` structural pattern, parameterized
  by `model::MGPModel` instead of hand-written per-model code. It calls
  `kli_hazard`, `kli_select`, `kli_decay`, `apply_pop`, `apply_move!` —
  three of those five (`kli_select`, `kli_decay`, `apply_move!`) are
  stubs, so **`regular_step!` cannot currently run to completion**: it
  would immediately throw on the first regular event once `total>0`.
  There is no singular-event dispatch loop at all in `mgp_filter.jl` —
  `singular_update!` (`mgp_filter.jl:120-137`) is declared (also a stub)
  but nothing calls it; no analogue of `POMP.pomp(...)`/`rinit`/
  `logdmeasure` wiring exists for `MGPModel` yet. There is no test that
  exercises `regular_step!` end-to-end (searched `test/`: no reference to
  `regular_step!`, `kli_select`, `kli_decay`, `apply_move!`, or
  `singular_update!` anywhere in `test/*.jl`).

## M. Model-specific hand-written likelihood terms

The bulk of the working likelihood math lives in:
- `src/examples/seir_funs.jl` (shared transmission!/progression!/recovery!/
  waning!/sampling rate-splitting helpers + `event_rates!`, used by
  `SoftSEIR`/`GuidedSEIR`... actually only referenced by
  `filter_pomp`-style code — see file header) and `seir_naive.jl` (its own
  independent, structurally-identical-but-separately-typed
  `singular_part!`/`event_rates!`/`regular_part!`).
- `src/examples/mers_funs.jl` (shared `knowledge!`, `singular_part!` for
  root/sample/fork events across Soft/Guided/Hard MERS) and
  `mers_naive.jl` (independent `singular_part!`/`event_rates!`/
  `regular_part!`).
- Each of `seir_soft.jl`/`seir_guided.jl`/`seir_hard.jl` and
  `mers_soft.jl`/`mers_guided.jl`/`mers_hard.jl` supplies its own
  `regular_part!` (the proposal-specific piece) while sharing the rest via
  `include`.

These are all *fully worked, closed-form* per-event weight formulas (e.g.
`seir_naive.jl:160` `ll += log(β*S*I/pop)` then `:173` `ll -= log(E*I)` for
a birth/fork event) — i.e. exactly the "hand-derived Φ_u formulas, skipping
the full/reduced distinction" the task background describes checking for.
There is no code path anywhere that computes a general φ_u via the stated
formula `φ_u = Q_u · Π_d C(n_d-ℓ_d,r_ud-s_d)/C(n_d,r_ud)` or a general Φ_u
by summing over the marginalized mark `m` — every one of these closed forms
was derived by hand per model/per event and is simply asserted correct by
citation to `mers_filter_suite.tex` equation numbers (see §"Event-Specific
Derivations", `mers_filter_suite.tex:468-541`) or to the (external,
not-in-repo) KLI paper directly.

## N. Existing tests / reference objects / oracles

`test/runtests.jl` includes, in order: `parse.jl`, `newick.jl`, `cblv.jl`,
`fsmarkov.jl`, `guide.jl`, `rcateg.jl`, `seir_macro_equivalence.jl`,
`seir_simulate.jl`, `seir_naive.jl`, `seir_soft.jl`, `seir_guided.jl`,
`seir_hard.jl`, `mers_naive.jl`, `mers_simulate.jl`, `mers_soft.jl`,
`mers_guided.jl`, `mers_hard.jl`.

Oracles currently in use, in increasing order of independence:
1. **Structural self-consistency**: `test/seir_macro_equivalence.jl`
   compares the `@mgp`-generated `SEIR.events` against the hand-typed
   `SEIR_REFERENCE` table (exact struct-field equality,
   `seir_macro_equivalence.jl:86-96`) and against `NaiveSEIR`'s
   independently-typed `event_rates!` driver/decay values over 1000
   randomized states (`seir_macro_equivalence.jl:113-134`). This validates
   the *DSL → Event table* step and the population-hazard/driver-rate
   arithmetic only — explicitly documented as *not* validating
   `mgp_filter.jl`'s move/singular-weight stubs
   (`seir_macro_equivalence.jl:11-12`: "This file does not claim
   full-filter equivalence").
2. **Cross-consistency across proposal kernels**: `test/mers_soft.jl`
   checks Soft vs. Guided diverge under a deliberately nonuniform guide
   (as expected, since they use different importance kernels); the
   suite's benchmark table (`mers_filter_suite.tex:1209-1244`) reports
   Soft/Guided/Hard converging to the same log-likelihood value (within
   Monte Carlo error) on a shared fixture — an informal three-filter
   agreement check, not a formal automated test asserting numeric
   closeness across all three files.
3. **Simulator/filter round-trip**: `test/seir_simulate.jl` and
   `test/mers_simulate.jl` — simulate a genealogy from the `@mgp`-declared
   model with `simulate(model,θ;...)`, then feed it into the
   *independently hand-coded* `NaiveSEIR.filter_pomp` /
   `SoftMERS.filter_pomp` and assert a finite log-likelihood
   (`seir_simulate.jl:81-90`). This is the strongest *cross-encoding*
   check that currently exists — but it only asserts finiteness, not a
   specific numeric likelihood value, and it validates the *simulator*
   against the *hand-coded filter*, not the MGP filter scaffold (which
   cannot run at all — see item L).
4. **External R oracle**: `scripts/seir_crossvalidate.jl` /
   `scripts/seir_crossvalidate.R` /
   `scripts/gen_seir_for_crosscheck.R` /
   `scripts/simulate_seir_forest.R` cross-check the Julia forward
   simulator's *output distribution* (tree shapes/CBLV encodings) against
   R `phylopomp::runSEIR` (commit `e5d6b30`'s message: "SEIRS/MERS sim
   test + check against R phylopomps runSEIR"). This is a simulator-only
   check (not a filter/likelihood check) and lives outside `test/
   runtests.jl` (run manually via the `scripts/` files, with results
   recorded in `seir_crossvalidate_results.md` and `results_7.23.md`), so
   it is not part of CI (`.github/workflows/CI.yml` runs
   `julia --project test/runtests.jl`-style Pkg tests only — not
   verified in detail here, but no `scripts/*.R` invocation appears in
   the CI workflow file names found).
5. No independent-Julia/R oracle exists yet for *filter log-likelihood
   values* (as opposed to simulator tree shapes) — flagged explicitly as
   unfinished in `handoff.md`'s "Remaining work" section: "Cross-language
   gold-standard validation (never implemented)."

No test anywhere references `MGPModel`'s filter path (`regular_step!`,
`kli_select`, etc.) — confirmed by grep across `test/*.jl`.

## O. Experimental / unverified / stub code

Explicitly self-labeled as unverified stubs (grep for `error("...stub...")`
in `src/examples/mgp_filter.jl`):
- `kli_select` (`mgp_filter.jl:70-72`)
- `kli_decay` (`mgp_filter.jl:87-89`)
- `apply_move!` (`mgp_filter.jl:116-118`)
- `singular_update!` (`mgp_filter.jl:135-137`)

All four are documented (each has a docstring citing specific KLI equation
numbers) but literally `error(...)` if called — this is the entire
"Reduction over event indicator m → Reduced KLI IR → Filter IR" portion of
the target compiler pipeline: **it does not exist yet**, only its intended
call sites and equation citations do.

The file's own header (`mgp_filter.jl:10-16`) states this directly: "STATUS:
compiler-checked scaffold... UNVERIFIED STUBS... deliberately not filled
in... Filling a stub does NOT make it correct until it is (a) verified
term-by-term against the cited equation and (b) cross-validated."
`mgp_filter.jl:190-201` lists an explicit, unimplemented "Validation gates"
checklist (Kingman coalescent / linear birth-death special cases, SIRS/SEIRS
worked-example reproduction, three-filter agreement, Chu-Vandermonde
collapse) that nothing in the repo currently satisfies for the MGP filter
path (items 2-4 partially exist for the *hand-coded* filters, per §N, not
for `mgp_filter.jl`).

`mgp.jl`'s own header (`mgp.jl:10-14`) is more optimistic and somewhat in
tension with `mgp_filter.jl`'s: it claims "STATUS: compiler-checked" for the
*event-table* generation (true, per `seir_macro_equivalence.jl`) but this
should not be read as saying the *filter* is checked — only the DSL→`Event`
table lowering is.

Other loose ends found (grep for `TODO`/`FIXME`):
- `src/genealogy.jl:10`: `## FIXME: inclusion of Root in NodeType introduces
  some inelegant redundancy`.
- `src/guide.jl:47,80`: two `FIXME: inelegant to store redundant
  information` (duplicated `target`/`dtarget` matrices on `GuideNode`).
- `TODO.md` (repo root, 5 lines): "iterated filtering", "parallelization of
  filters (with explicit control of RNG)", "general phylopomp filter
  construction", "simple Euler multinomial filter", "genealogical
  simulations" — the last item is now done (`src/simulate.jl`, milestone 1
  in the *other* numbering, see the warning at the top of this document);
  the rest remain open and are orthogonal to the compiler project.
- `handoff.md` "Remaining work" (end of file) flags the pre-existing
  `mers_naive.jl` double-susceptible-depletion crash (`AssertionError:
  S_c > 0`/`S_h > 0` instead of `-Inf`) as a known, out-of-scope defect —
  still present at HEAD (not touched by commit `931af3b`).
- `mers_kernel_diagnostic_findings.md` documents a genuine
  target-incompatibility finding (Φ=0 from near-simultaneous camel sample
  tips causing deterministic `-Inf` on the real 274-tip empirical tree) —
  a data/tree artifact, not a code defect, but relevant context for anyone
  validating filter correctness against the real MERS tree specifically.

## Source-of-truth material found (not in the git repo)

- The KLI paper itself, `StructuredMGPs.pdf` (cited throughout `mgp.jl`,
  `mgp_filter.jl`, and `mers_filter_suite.tex` as "King, Lin & Ionides,
  'Exact phylodynamic likelihood via structured Markov genealogy
  processes,' arXiv:2405.17032"), is **not checked into this repository**.
  It was located on this machine at
  `/home/dislam/University of Michigan Dropbox/Deepan Islam/Summer-2026/
  phylodynamics_companion/Literature/StructuredMGPs.pdf` — outside the git
  working tree, so `git log`/`git blame` cannot see it and a fresh clone of
  this repo would not have it. A stray `/tmp/StructuredMGPs_layout.txt`
  (outline/layout notes, presumably from a prior session) also exists
  outside the repo and outside the scratchpad convention; not part of the
  codebase.
- `src/examples/mers_filter_suite.tex` / `.pdf` (1253 lines / 23 pages) is
  the repo's own from-scratch derivation of the MERS filter, citing KLI
  equations throughout. It is the most complete worked derivation
  in-repo, but it derives *closed-form, model-specific* Φ_u/proposal
  formulas by hand (§"Event-Specific Derivations",
  `mers_filter_suite.tex:468-541`) rather than deriving them mechanically
  from the general φ_u/Φ_u formulas in the task background — i.e. it is a
  worked example to check a compiler against, not a compiler itself, and
  not obviously structured to make that mechanical derivation easy to
  extract.
- No other KLI-adjacent notes/derivations were found in the repo beyond
  `mers_filter_suite.tex` and the `mgp_filter.jl` docstring citations.

## Summary diagram of what exists today

```
Hand-coded path (WORKING, TESTED):
  seir_funs.jl / mers_funs.jl (hand-derived Φ_u, π_u, decay per event)
    -> seir_{naive,soft,guided,hard}.jl / mers_{naive,soft,guided,hard}.jl
       (POMP.pomp rinit/rprocess/logdmeasure)
    -> POMP.pfilter (SMC loop, from PartiallyObservedMarkovProcesses pkg)
    -> tested in test/{seir,mers}_{naive,soft,guided,hard}.jl

MGP scaffold path (PARTIAL, UNTESTED PAST STEP 1):
  @mgp DSL (mgp_macro.jl)
    -> MGPModel/Event table (mgp.jl) [VERIFIED vs SEIR_REFERENCE + NaiveSEIR
       driver/decay, via seir_macro_equivalence.jl]
    -> simulate(model,θ;...) (simulate.jl) [WORKS — cross-checked against
       NaiveSEIR/SoftMERS filters for finiteness, and against R runSEIR for
       simulator output shape]
    -> mgp_filter.jl's regular_step! [CANNOT RUN — kli_select, kli_decay,
       apply_move! are error() stubs]
    -> no singular-event loop, no pomp()-style wiring, no test coverage
       past the Event-table level
```
