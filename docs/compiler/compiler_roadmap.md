# Compiler roadmap: desired passes vs. what exists today (M00)

This maps the target pipeline

```
Model DSL -> Population IR -> Full Genealogical/KLI IR -> Reduction over m
  -> Reduced KLI IR -> Filter IR -> {Julia code, audit text, LaTeX/table, verification}
```

onto concrete passes and states, against what was found during the M00
reconnaissance (see `docs/compiler/architecture_before.md` for file:line
detail). Pass names below are descriptive placeholders, per the task
prompt ("names are flexible") — nothing in the repo currently uses these
names except where noted.

For each pass: **EXISTS** (file/function given), **PARTIAL** (what's there,
what's missing), or **NOT PRESENT**.

## 1. `parse_model` (DSL text/macro -> raw model description)

**EXISTS**, conflated with lowering. `@mgp`/`@event` in
`src/examples/mgp_macro.jl` parse a `begin...end` block's declarations
(`compartments=`, `demes=`, `params=`) and `@event` calls
(`mgp_macro.jl:108-136`) directly into `Event(...)` constructor call
expressions — there is no intermediate "raw AST" or "parsed-but-not-lowered"
representation distinct from the final `MGPModel`/`Event` structs. For this
compiler project's purposes, `parse_model` and `lower_population` (item 2)
are currently the *same pass*, executed at macroexpansion time. Splitting
them (e.g. into a real intermediate parse tree that a `validate_model` pass
can inspect before lowering) is new work, not a refactor of something
broken — the current single-pass version works for its narrow scope
(verified structurally equal to a hand-written reference table by
`test/seir_macro_equivalence.jl`).

Two independent hand-written models exist as `@mgp` output: `SEIR`
(`mgp.jl:188-198`) and `MERS` (`mgp_mers.jl:3-24`). No third model has been
run through the DSL, so its generality beyond these two is unproven.

## 2. `lower_population` (raw model -> Population IR: Δ, α, r per event)

**EXISTS**, as the `Event`/`MGPModel` structs themselves (`mgp.jl:59-82`).
This *is* the Population IR the task background describes: per-event static
data Δ_u (`Event.Δ`), α_u (`Event.hazard`), r_u (`Event.r`), from/into
wiring W_u (`Event.from`/`Event.into`), event type, regular/singular flag,
observed flag. τ_u (event time) is deliberately absent from the static
struct (it's a runtime quantity, correctly not modeled here — see
`architecture_before.md` §B).

What's missing for a *complete* Population IR: no explicit initial-condition
IR (x0 is passed ad hoc as a `NamedTuple` kwarg to `simulate`/`filter_pomp`,
never validated against `MGPModel.compartments`); no representation of
compound events (KLI §3.2(f), explicitly out of scope per `mgp.jl:30`); no
IR-level distinction between the population-level demes (`compartments`)
and lineage-carrying demes (`demes`) beyond the two separate `Vector{Symbol}`
fields — i.e. no validation pass confirming every event's `r`/`from`/`into`
indices are consistent with `demes` (this is currently done ad hoc inside
`_decode_move`/`_production` at macro-expansion time,
`mgp_macro.jl:40-45,61-67`, with no separate audit step).

## 3. `validate_model` (structural checks on Population IR)

**PARTIAL**, folded into macro-expansion, not a separate inspectable pass.
`@mgp`'s expansion does perform validation — duplicate-declaration checks
(`mgp_macro.jl:161`), missing-declaration checks (`mgp_macro.jl:170`),
deme-must-be-compartment checks (`mgp_macro.jl:176`),
compartment/parameter-name-disjointness (`mgp_macro.jl:177-178`), and
per-event `rate=`/`move=` presence (`mgp_macro.jl:121-122`) — but these are
`error()` calls thrown *during macroexpansion*, not a callable
`validate_model(model::MGPModel) -> Vector{Issue}` function that could be
run against an already-constructed `MGPModel` (e.g. one built by hand, like
`SEIR_REFERENCE`) or produce a structured audit report. There is no check
anywhere that hazards are non-negative for all reachable states, that Δ
never drives a compartment negative, that every `BIRTH`/`MIGRATION` event's
`from` deme is actually in `demes`, etc. — i.e. no *semantic* validation
beyond what the DSL's own syntax-level checks catch.

## 4. `enumerate_saturations` (per event, enumerate compatible s given ℓ, n)

**NOT PRESENT** anywhere in generic form. The concept (saturation s(Y_t),
KLI §3.4.2) is named only in `Event`'s docstring (`mgp.jl:52-54`, explicitly
as something *not* stored, computed only at runtime) and cited in
`mers_filter_suite.tex` §"Step B: enumerate saturations"
(`mers_filter_suite.tex:443-449`), which walks through the enumeration *by
hand* for each MERS event type rather than implementing a generic
enumeration procedure. No Julia function exists that, given an `Event` and
a lineage-count/coloring state, enumerates the set of compatible
saturations. The hand-coded filters implicitly enumerate small, fixed sets
of outcomes per event via literal `rcateg([...])` calls over 2-3 explicit
options (e.g. `seir_naive.jl:161`, `mers_naive.jl:65,95`) — these are the
*results* of an enumeration a human did on paper, not a runtime enumeration
procedure.

## 5. `derive_full_compatibility` (Q_u, full-coloring compatibility indicator)

**NOT PRESENT** generically. `Q_u` (KLI Eq. 38) is named in
`apply_move!`'s docstring (`mgp_filter.jl:107`) and derived by hand per
MERS event in `mers_filter_suite.tex` §"Step D: the coloring operator"
(`mers_filter_suite.tex:459-467`). No Julia code computes a general
compatibility indicator from `Event` + genealogy state. The hand-coded
filters' `-Inf`-on-incompatibility branches (e.g. `seir_naive.jl:47-52`,
`n.lineage ∉ cols[Infec]` triggering an inconsistency correction) are
per-event hard-coded special cases of what a general Q_u check would be,
not an implementation of Q_u itself.

## 6. `compute_phi` (φ_u, and Φ_u after marginalizing m)

**NOT PRESENT** generically; **EXISTS as hand-derived closed forms** per
model/event. This is the single largest gap relative to the task
background's stated concern. Every one of the eight SEIR/MERS filter
modules computes something numerically equivalent to Φ_u — e.g.
`seir_naive.jl:160,173` (`ll += log(β*S*I/pop)`, `ll -= log(E*I)` for a
birth/fork) or `mers_naive.jl:77` (`ll += log(lambda_cc) -
log(Ic*(Ic-1)/2)`) — but every one of these was derived on paper (see
`mers_filter_suite.tex` §"Event-Specific Derivations",
`mers_filter_suite.tex:468-541`, and the Complete Reference Table at
`mers_filter_suite.tex:541-638`) using the Chu-Vandermonde-collapsed,
already-marginalized-over-m formula for that *specific* event type, then
hand-transcribed into Julia. No code computes the general
`φ_u = Q_u · Π_d C(n_d-ℓ_d, r_ud-s_d)/C(n_d,r_ud)` from `Event` data and a
saturation, and no code performs the sum over full-coloring transitions
`Φ_u(d,d') = Σ_m φ_u` generically — both are asserted correct by citation
to a hand-derivation, not computed. This is exactly the situation flagged
as a risk in the task background: "check whether the current code computes
something like Φ_u directly... via hand-derived model-specific formulas" —
confirmed true for all eight filter modules; `mgp_filter.jl`'s `apply_move!`
stub is the *only* place in the codebase where a generic Φ_u/boost
computation is even attempted, and it is unimplemented.

## 7. `reduce_event_indicator` (marginalize m; Full KLI IR -> Reduced KLI IR)

**NOT PRESENT** as a distinct pass; implicitly assumed complete everywhere.
No code in this repository ever represents or manipulates the full
coloring `Y=(d,m)` — `Coloring{D,N}` (`coloring.jl:10-19`) tracks only `d`.
Every consumer (hand-coded filters and the `mgp_filter.jl` scaffold alike)
works directly in the reduced space, i.e. the marginalization over `m` is
assumed to have already been carried out (correctly, by hand, per the
`mers_filter_suite.tex` derivations) rather than performed or checked by
any code. There is no representation of `m` to reduce *from* — so there is
currently no way to even state, in code, what "before reduction" would
look like, let alone verify the reduction step (e.g. via the
Chu-Vandermonde identity, KLI Eq. 8, which `mgp_filter.jl:199` lists as an
unimplemented "Validation gate": "Chu-Vandermonde collapse (Eq. 8) verified
where m/s sums are marginalized — not left as permanently separate
'inline-node inserted' states").

## 8. `build_filter_ir` (regular/singular/decay classification with
   provenance)

**PARTIAL**. `Event.regular::Bool` (`mgp.jl:67`) already gives a
per-event regular/singular classification, and `Event.observed::Bool`
(`mgp.jl:68`) flags which events can produce genealogy features — this is
real, working IR-level metadata, checked by
`test/seir_macro_equivalence.jl:86-96` against the hand table. What's
missing is the finer three-way split the task background asks for
(regular flow / singular flow / decay term) *with provenance back to the
source event* as an explicit, inspectable structure: currently "decay" is
just a running scalar accumulator (`decay::Prob` in
`seir_naive.jl:103-121`'s `event_rates!`, similarly in `mers_naive.jl:132-
162`) that sums contributions from several events without retaining which
event contributed how much — i.e. decay exists as a *number*, not as a
*term with provenance*. `mgp_filter.jl`'s `kli_decay` stub
(`mgp_filter.jl:74-89`) is meant to fill this generically but is
unimplemented. No "Filter IR" data structure (as opposed to a bare `Float64`
log-likelihood accumulator) exists anywhere in the repo.

## 9. `build_proposal` (π_u construction, importance-kernel design)

**EXISTS**, richly, but per-model/per-kernel, not from Population IR
generically. Four independent proposal families exist for each of SEIR and
MERS (naive/soft/guided/hard — see `architecture_before.md` §J), sharing
the reverse-time guide machinery in `src/guide.jl`
(`Guide`/`GuideNode`/`relhaz`/`choose_branch`,
`src/fsmarkov.jl`'s `FSMarkovProc`). This is real, tested, working
importance-sampling design — but it is hand-selected and hand-coded per
event per model, not derived mechanically from `Event`/`MGPModel` data by a
`build_proposal` pass. `mgp_filter.jl`'s `kli_select` stub
(`mgp_filter.jl:56-72`) is the intended generic entry point and is
unimplemented. Note the four-kernel taxonomy (naive/soft/guided/hard)
itself is a solid, reusable *design vocabulary* already validated
end-to-end for two models — a `build_proposal` pass should probably
formalize this existing taxonomy rather than invent a new one.

## 10. `compile_filter` (assemble Filter IR -> executable Julia code)

**PARTIAL**, in two disconnected forms. The hand-coded path compiles
"by hand" into `POMP.pomp(...)` objects with `rinit`/`rprocess`/
`logdmeasure` closures (e.g. `seir_naive.jl:195-253`) — this is real,
working, generated-per-model Julia code, but the "compilation" step is a
human writing Julia directly, not a pass consuming an IR. The MGP scaffold
has the skeleton of a *generic* compiled filter —
`regular_step!` (`mgp_filter.jl:163-188`) is genuinely model-agnostic
(parameterized by `model::MGPModel`), so structurally it is exactly what
`compile_filter` should produce, but it does not run (three of its five
callees are stubs) and there is no singular-event loop or `POMP.pomp(...)`
wiring around it at all — i.e. the "assembled executable" half of this pass
does not exist yet.

## 11. `audit` (audit text output)

**PARTIAL**, as documentation/docstrings rather than a generated artifact.
`mgp.jl`'s and `mgp_filter.jl`'s extensive docstrings (citing specific KLI
equation numbers per function) function as a manually-written audit trail,
and `test/seir_macro_equivalence.jl` is itself a form of executable audit
(macro output vs. hand-written reference). No code generates an audit
report as an artifact (e.g. "for model X, event `infection` maps to KLI
Eq. Y via Q_u=..., φ_u=..."). `mers_filter_suite.tex` is the closest thing
to a generated/curated audit document but it is hand-written LaTeX, not
compiler output tied programmatically to `MGPModel`/`Event` data (i.e. if
`SEIR`'s hazard formula changed, nothing would flag that
`mers_filter_suite.tex` — or its SEIR references — is now stale).

## 12. `LaTeX/table` output

**NOT PRESENT** as generated output; **EXISTS** as hand-written LaTeX.
`mers_filter_suite.tex` (1253 lines, compiles to a 23-page PDF, see its
own "From Math to Code" §13, `mers_filter_suite.tex:1133-1245`) is a
complete, carefully-maintained, hand-written document with a real
changelog discipline (v26 through v31 tracked in-document) — a strong
existing convention to build a generator *against*, but nothing generates
LaTeX/tables from `MGPModel`/`Event` data programmatically.

## 13. `verify` (validation gates / cross-checks)

**PARTIAL**, strong for the hand-coded path, essentially absent for the
MGP scaffold. See `architecture_before.md` §N for the full oracle
inventory: structural macro-equivalence (real, automated), simulate/filter
round-trip finiteness (real, automated), Soft-vs-Guided divergence and
Soft/Guided/Hard convergence (real, but the convergence check is informal
— a benchmark table, not an automated numeric-closeness assertion), R
`runSEIR` simulator cross-check (real, but simulator-only, and outside
`test/runtests.jl`/CI). `mgp_filter.jl:190-201`'s own "Validation gates"
list (Kingman coalescent/Moran, linear birth-death, SIRS/SEIRS worked
examples, three-filter agreement, Chu-Vandermonde collapse) is an
explicit, unimplemented checklist for what a `verify` pass would need to
run once the filter path exists — currently zero of these five gates has
any code behind it for the `MGPModel`/`mgp_filter.jl` path specifically
(items 2-4 are informally satisfied for the *hand-coded* filters only).
No cross-language (R) oracle exists for filter *log-likelihoods*
specifically (only for simulator output shape) — flagged as "never
implemented" in `handoff.md`.

## Summary table

| Pass | Status | Where |
|---|---|---|
| parse_model | EXISTS (fused w/ lowering) | `mgp_macro.jl` |
| lower_population | EXISTS | `mgp.jl` (`Event`/`MGPModel`) |
| validate_model | PARTIAL (inline, not callable/inspectable) | `mgp_macro.jl:108-186` |
| enumerate_saturations | NOT PRESENT (hand-done in tex) | `mers_filter_suite.tex:443-449` (by hand) |
| derive_full_compatibility (Q_u) | NOT PRESENT (hand-done in tex) | `mers_filter_suite.tex:459-467` (by hand) |
| compute_phi (φ_u, Φ_u) | NOT PRESENT generically; hand-coded per model | `seir_naive.jl`, `mers_naive.jl`, etc. |
| reduce_event_indicator (marginalize m) | NOT PRESENT (assumed done, unverified) | n/a — no full-coloring repr. anywhere |
| build_filter_ir (regular/singular/decay+provenance) | PARTIAL | `Event.regular`/`.observed`; decay is a bare scalar |
| build_proposal (π_u) | EXISTS per-model, 4 kernels x 2 models | `seir_{soft,guided,hard}.jl`, `mers_{soft,guided,hard}.jl`, `guide.jl` |
| compile_filter | PARTIAL (hand path works; generic path stubbed) | `seir_naive.jl` (hand); `mgp_filter.jl:163-188` (stub-blocked) |
| audit | PARTIAL (docstrings only) | `mgp.jl`, `mgp_filter.jl` docstrings |
| LaTeX/table generation | NOT PRESENT (hand-written doc exists) | `mers_filter_suite.tex` |
| verify | PARTIAL (strong for hand path, absent for MGP path) | `test/seir_macro_equivalence.jl`, `test/*_simulate.jl`; `mgp_filter.jl:190-201` checklist unimplemented |

## Where the real work is

Reconnaissance suggests the compiler project's hardest, highest-value work
is items 4-7 (`enumerate_saturations`, `derive_full_compatibility`,
`compute_phi`, `reduce_event_indicator`) — nothing generic exists for any
of them, they are the actual mathematical core of the KLI reduction, and
`mers_filter_suite.tex`'s hand derivations (§"Construction of Compatibility
Terms", `mers_filter_suite.tex:434-541`) are the best available worked
example to derive a general procedure against and check it. Items 1-3
(parsing/lowering/validation) are close to done and mostly need
restructuring (splitting fused steps, making validation callable) rather
than new mathematics. Items 8-13 (filter IR, proposal, compile, audit,
LaTeX, verify) have strong *per-model* existing implementations
(especially proposal design) that a generic pass should aim to reproduce
and unify, not replace wholesale — the four-kernel naive/soft/guided/hard
taxonomy in particular looks like a design to keep, not redo.
