# Docstring sweep: line-level detail (2026-10-05)

Companion to `review_simulator_2026-10-05.md`, Part 2. Four read-only reviewers, one rubric.

**Applied 2026-10-05.** Line numbers below refer to the pre-cut files, now stale. Pre-cut state: `git diff refs/review/pre-docstring-cuts -- src test`.

## Removed note worth keeping

From `test/mtbd.jl` (274-tip MERS tree, MTBD filters return -Inf). The test comment now keeps only the first two sentences.

> Exact log likelihood at mers_mle with a 1-year stem: +381.8742. Both filters return -Inf here (0 of 6 runs finite at Np = 2000, for each of several guides). Measured: the first non-finite observation depends on the seed (102, 144 and 147 in three runs); effective sample size falls to about 1 at branch points well before the first tip; particles carry about 70 untracked humans. Hypothesis, not yet isolated: untracked hosts are simulated without sampling and pay for it through the exp(-psi*I*dt) survival weight, which with psi2 = 3.05/yr makes the weights very uneven. If so, the remedy is to propose untracked hosts conditioned on leaving no sampled descendants, using the ODE's p_i(t).

Format: `path:start-end | delete|compress | replacement (≤2 lines, "—" for delete) | rubric item(s)`

Rubric items: 1 history · 2 milestone/Gate/handoff/scratch · 3 verification claims · 4 R/C++ or other-file/tex line refs · 5 ALL-CAPS, banners, self-justification · 6 restates code · 7 derivation essay.

---

## A. `src/examples/` compiler files with the heaviest comments

### mgp_filter_ir.jl — rewrite — ≈277 of 348
- `1-48` | delete | — | 1,2,5 (M05/M06 bucket spec, "CORRECTION (M06, Part 0)", tex cites, handoff path, banner)
- `50-125` | compress | `# Fork outcomes of a regular event are never regular flow: the filter realizes forks only at observed branch points. They are classed as OutflowImbalanceTerm (KLI Eq. 45, 47).` | 1,2,4,5,7
- `130-140` | compress | `Abstract supertype of RegularFlow, SingularFlow and OutflowImbalanceTerm. Each carries event, reduced and Φ.` | 2
- `143-152` | compress | `A :noop or :cross outcome of a regular BIRTH/MIGRATION event (KLI Eq. 45). :fork is never regular flow.` | 2,5
- `159-182` | compress | `An outcome of a singular event (event.regular == false). Realized only at the observed data-event time.` | 2,4
- `189-240` | compress | `Fork mass of a regular event that cannot be regular flow. Φ is copied from the ReducedTransition. mechanism is always :inflow_outflow_imbalance; reason is :fork_unobserved_at_regular_time.` | 1,2,4,5,7
- `249-265` | compress | `One event's classified terms: regular, singular and outflow_imbalance vectors.` | 1,2
- `273-303` | compress | `Classify each reduced transition. If event.regular, :noop/:cross go to RegularFlow and :fork to OutflowImbalanceTerm; otherwise all go to SingularFlow. Throws ArgumentError if a transition's event is not event.` | 2,4
- `339-345` | compress | `Compose reduced_transitions and classify_filter_terms for (event, ℓ, n).` | 6

### mgp_mers_filter.jl — rewrite — ≈269 of 511
- `1-26` | compress | `# MERS filter built from the compiler IR. The regular part is compiled. The singular part reuses NaiveMERS.singular_part!.` | 2,5
- `27-64` | delete | — | 1,2,3,4,5 (FINDING 1; "CONFIRMED, NOT JUST ASSUMED")
- `65-141` | compress | `# TCC/THH: no-move weight is Φ_id + ℓ_d·Φ_inl (C(ℓ,s)-weighted sum). THC/TCH cross: boost(Φ_cross, 1/ℓ).` | 1,2,3,4,7 (FINDING 2; bug narrative at 92-98)
- `142-178` | delete | — | 1,2,3,4,5 (FINDING 3: workaround, "RESOLUTION", "v28", "bit-exact Gate-5")
- `182-198` | compress | `Return the .phi of the T transition in ts, or 0//1 if absent. A slot-bearing transition is absent when ℓ_d == 0.` | 7
- `204-219` | compress | `Fill alpha/pi_ for the 12 MERS events and return the total decay (compiled_decay). Slots: 1 TCC, 2 THH, 3/4 THC no-move/cross, 5/6 TCH no-move/cross, 7/8 removal, 9/10 birth, 11/12 death.` | 1,2,4
- `251-254` | delete | — | 1,2
- `275-308` | compress | `Compiled counterpart of NaiveMERS.regular_part!. Same control flow and RNG call order, so a fixed seed gives the same trajectory. Returns (ll, Sc, Ic, Sh, Ih).` | 2,4,6
- `347-348` | compress | `# Φ_id + ℓ_c·Φ_inl` | 2,6
- `367-368` | delete | — | 6
- `430-445` | compress | `Build the MERS filter POMP for genealogy gen. Singular part is NaiveMERS.singular_part!; regular part is mers_compiled_regular_part! with model = MERS.` | 2,4 (line 442 also names a function in the wrong file)

### mgp_seir_filter.jl — rewrite — ≈199 of 400
- `1-30` | compress | `# SEIR filter built from the compiler IR. The regular part is compiled. The singular part reuses NaiveSEIR.singular_part!.` | 1,2,5
- `31-117` | delete | — | 1,2,3,4,5,7 (decay-discrepancy essay, scratch probes, "ALL MATCH: true")
- `122-132` | compress | `Return total_decay plus, for each DEATH event, the leftover α_full·1{n>ℓ} − α_reduced, where α_reduced = α_full·(n−ℓ)/n. Equals γ·ℓ, the naive decay. ℓ and n are per-deme post-event counts in model.demes order.` | 1,2,4
- `150-162` | compress | `Fill alpha/pi_ for the 6 SEIR events and return the decay. Slots: 1/2 infection no-move/cross, 3/4 progression identity/cross, 5 recovery, 6 waning. χ is the culling rate; culling is singular, so it enters only the decay.` | 4,5 ("cited inline" at 158)
- `187-189` | compress | `# NaiveProposal: infection splits no-move : move by no_move_share; progression keeps the source-deme split.` | 4
- `201-227` | compress | `Compiled counterpart of NaiveSEIR.regular_part!. Same control flow and RNG call order, so a fixed seed gives the same trajectory. Returns (ll, S, E, I, R).` | 2,3,4,6
- `259-266` | compress | `# No move: Φ_id + ℓ_I·Φ_inl = 1 − ℓ_E/E'. Divided by pi_[1] in the ll update above.` | 1,4
- `279-283` | compress | `# CrossDemeTransition is a singleton per event and state.` | 2,4
- `292-295` | compress | `# progression has no slot in its own deme, so only IdentityTransition occurs.` | 6
- `330-340` | compress | `Build the SEIR filter POMP for genealogy gen. Singular part is NaiveSEIR.singular_part!; regular part is compiled_regular_part! with model = SEIR.` | 4,5

### mgp_mgpaudit.jl — rewrite — ≈217 of 304
- `1-53` | compress | `# mgpaudit / @mgpaudit: per-mark derivation report (audit_model -> full_transitions -> reduce_event_indicator -> classify_filter_terms -> explain). Included after mgp_explain.jl because the report structs name its types.` | 2,5,7
- `58-94` | compress | `Default (ℓ, n) for mgpaudit: SEIR ([2,2],[6,5]), MERS ([2,2],[5,5]), else fill(2,D)/fill(5,D). Chosen so every BIRTH mark reaches a Fork saturation.` | 2,4,5
- `106-113` | compress | `One event's row in an MgpAuditReport: MarkAuditInScope or MarkAuditOutOfScope.` | 2,5
- `116-125` | compress | `Derivation of a BIRTH/MIGRATION event at (ℓ, n): event_audit, reduced explanations, the FilterSpec and one TermExplanation per term.` | 2
- `135-143` | compress | `Record for a DEATH/SAMPLE/NEUTRAL event. note says why the event is outside the saturation chain.` | 2,4
- `149-156` | compress | `Result of mgpaudit: model_name, the (ℓ, n) state, and one MarkAudit per event in model order.` | 5
- `164-176` | compress | runtime string → `out of scope for the saturation chain: $(ev_type) events are decay/rate-driven, not saturation-driven.` | 2,4 (check tests for text matches first)
- `178-202` | compress | `Derive every event of model at (ℓ, n) (default_audit_state if omitted). Returns an MgpAuditReport. Throws ArgumentError unless length(ℓ) == length(n) == length(model.demes) and ℓ .<= n.` | 2,5,7
- `243, 251` | compress | runtime strings → `reduced transitions: ` / `classification: ` | 2
- `277-292` | compress | `@mgpaudit Model [ℓ=[...] n=[...]] expands to mgpaudit(Model; ℓ=..., n=...) and returns its MgpAuditReport.` | 5

### mgp_explain.jl — trim — ≈130 of 259
- `1-36` | compress | `# explain(x) returns a printable provenance record for a FilterTerm, ReducedTransition or KLITransition (FilterTerm -> ReducedTransition -> KLITransition -> Event). Read-only.` | 2,5,7
- `40-45` | compress | `# Kind tag and detail string per KLITransition subtype.` | 5
- `62-64`, `113-115`, `169-171` | delete | — (banners) | 5
- `69` | compress | drop "(M03)" | 2
- `95-96` | compress | `t.event is the source Event.` | 2
- `120-123` | compress | drop "(M04)" and "(M03's phi_u values ...)" | 2
- `141-143` | compress | `... collapsed into it by reduce_event_indicator (Φ_u = Σ φ_u).` (drop tex lines) | 2,4
- `176-186` | compress | `Structured explanation of a FilterTerm: bucket, source mark, Φ and the ReducedExplanation. mechanism/reason are nothing except for OutflowImbalanceTerm.` | 1,2
- `211-217` | compress | `For OutflowImbalanceTerm, also reports mechanism (:inflow_outflow_imbalance) and reason (:fork_unobserved_at_regular_time).` | 2
- `248` | compress | runtime string → `" (not KLI's λ)"` | 4

---

## B. Other `src/examples/mgp*.jl` files

### mgp_decay.jl — rewrite — ≈224 of 281
- `1-31` | delete | — | 1,2,4,7
- `32-80` | compress | `# Decay λ(t,x,y), KLI Eq. 47 / B2. DEATH: α·1{n[d] ≤ ℓ[d]}, d = event.from. SAMPLE: the full α. BIRTH, MIGRATION, NEUTRAL: 0.` | 4,5,7
- `82-120` | delete | — | 2,3,4,7
- `122-160` | delete | — | 1,2,3,4,5,7 ("DISCREPANCY FOUND"; contradicts mgp_seir_filter.jl, which records it as resolved)
- `164-182` | compress | `Contribution of one DEATH/SAMPLE event to λ. Fields: event, alpha (hazard), gate (1//1 or 0//1), value = alpha*gate, reason (:sampling_hazard_full or :sub_threshold_removal).` | 5,6
- `191-221` | compress | `SAMPLE: gate = 1. DEATH: gate = 1 if n[d] ≤ ℓ[d], d = event.from, else 0. Exact if x, θ are Rational. Throws ArgumentError for other event types, event.from == 0, d out of range, ℓ[d] > n[d].` (keep signature lines 192-193) | 2,4,5,7
- `258-271` | compress | `λ(t,x,y): sum of decay_contribution(...).value over DEATH and SAMPLE events; other types are skipped.` | 4,5,7

### mgp_transitions.jl — rewrite — ≈170 of 253
- `1-49` | compress | `# Full KLI transitions (Eq. 9) for BIRTH and MIGRATION events. DEATH, SAMPLE, NEUTRAL events have no saturation structure (see mgp_decay.jl).` | 2,4,5,7
- `54-78` | compress | `Classified outcome of one saturation s of an event. Fields: event; s (saturation); phi (exact Rational{Int}). IdentityTransition and InlineSameDemeTransition both leave the reduced coloring unchanged but stay distinct (event indicator differs).` | 2,4
- `85-86` | compress | `Y' = Y.` | 4
- `101-102`, `117-118`, `135-136` | delete | — | 4
- `152-173` | compress | `Placeholder subtype. Never constructed: full_transitions throws for DEATH, SAMPLE, NEUTRAL.` (or remove the type) | 1,2,4,5,7
- `189-197` | compress | `Classify each saturation s by sum(s): 0 -> Identity; 1 -> InlineSameDeme if slot deme == event.from, else CrossDeme; ≥2 -> Fork.` | 2,4,6
- `205-208` | delete | — | 3,4,5
- `213-219` | compress | error text → `"full_transitions: event `$(event.name)` has type $(event.type); only BIRTH/MIGRATION are supported"` | 2,4

### mgp_reduce.jl — rewrite — ≈137 of 181
- `1-63` | compress | `# Group full_transitions output by reduced (color-only) outcome and sum phi into Φ_u (KLI Eq. 9; Chu-Vandermonde, Eq. 8).` | 1,2,4,5,7
- `70-74` | compress | `Result of grouping KLITransitions that give the same outcome on the color-only state.` | 2,4
- `86-95` | compress | `kind: :noop, :cross or :fork (same as key[1]). Φ: exact Rational{Int} sum of phi over the group. transitions: the grouped KLITransitions, in input order.` | 2,5
- `107-113` | compress | `Grouping key for t; see ReducedTransition for the shapes.` | 5,7
- `129-135` | compress | `Group transitions (full_transitions output for one event/state) by reduced_key. Φ is the exact sum of phi per group.` | 2,4
- `144-149` | compress | `Example: MERS THC gives 4 transitions; Identity and InlineSameDeme merge into (:noop, 1), giving 3 reduced transitions.` | 6,7
- `174-177` | compress | `full_transitions(event, ℓ, n; Q) followed by reduce_event_indicator.` | 2

### mgp_phi.jl — rewrite — ≈132 of 190
- `1-30` | compress | `# Saturations S_u(ℓ) (KLI Sec. 3.4.2) and the production-slot binomial ratio φ_u (Eq. 9).` | 1,2,4,5,7
- `40-42`, `62-63`, `108-110` | delete | — | 1,2,3,5
- `67-71` | compress | `ℓ_d = 0 or r_d = 0 collapses deme d's range to {0}.` | 2,4
- `92-101` | compress | `Base.binomial extends to negative a (binomial(-1,1) == -1); this returns 0 there, as for b < 0 and b > a.` | 3,7
- `111-121` | delete (move to `test/kli_phi_test.jl`) | — | 3,6 (load-time `@assert`s)
- `135-138` | compress | `Q is the compatibility indicator (KLI Eq. 38), supplied by the caller.` | 2,5
- `146-165` | compress | `Returns 0 if n_d < r_d for some d or if any s_d < 0. Throws ArgumentError if the vectors differ in length.` | 1,2,3,5

### mgp_proposal.jl — rewrite — ≈103 of 132
- `1-38` | compress | `# Proposal-strategy tags, and the driver (KLI Eq. 45) and boost (Theorem 5) compositions.` | 1,2,4,5,7
- `46-58` | compress | `Tag for a hand-designed proposal family: NaiveProposal, SoftProposal, GuidedProposal, HardProposal. π_u is not computed generically.` | 2,5,7
- `62-64` | compress | `Naive family: lineage-count-proportional selection.` | 4
- `70-74` | compress | `Guided family: π from a reverse-time sweep over the deme process.` | 4
- `83-94` | compress | `Driver rate β_u = α_u · π_u (KLI Eq. 45).` | 4,5,7
- `101-121` | compress | `Boost B_u = Φ_u / π_u (KLI Theorem 5). π is the total probability of proposing the reduced outcome, marginalized over within-outcome choices.` | 2,3,4,5
- `128-130` | compress | `Same as boost(rt.Φ, π).` | 2

### mgp_filter.jl — rewrite — ≈145 of 201
- `1-23` | compress | `# Generic KLI filter over an MGPModel. The weight-bearing slots (kli_select, kli_decay, apply_move!, singular_update!) are not implemented and throw.` | 3,5
- `25-27`, `52-55`, `139-143` | delete | — (banners) | 3,5
- `32-34` | compress | `Population hazard αᵤ(t,x) of ev as Float64 (KLI §2.5).` | 5
- `41-42` | compress | `Return x with ev.Δ applied.` (duplicates `apply_delta` in `src/simulate.jl`) | 5
- `59-68` | compress | `Selection factor πᵤ for regular event ev (KLI Eqs. 44-46); driver is αᵤ·πᵤ. Not implemented; throws.` | 3,4,5
- `77-85` | compress | `Decay rate λ(t,x) between genealogy events (KLI Eq. 47, B2). Not implemented; throws.` | 3,4
- `94-114` | compress | `Apply the coloring move for ev; returns log(phi_u) - log(q). Not implemented; throws. See KLI Eqs. 7-9, 37, 38, 46.` | 1,3,7
- `123-133` | compress | `Singular update at an observed genealogy event (KLI Eq. 14, B3): fork κ, chop χ or swap σ by node type. Roots enter through the initial condition. Not implemented; throws.` | 3,4
- `152-161` | compress | `Forward SMC filter over [t, t+dt) (KLI Theorem 5, Algorithm B1). Proposals differ only in π and the coloring proposal.` | 5,7
- `190-201` | delete | — | 2,3,5,7 ("Validation gates")
- error strings at `71, 88, 117, 136` say "unverified stub" | 3

### mgp_audit.jl — trim — ≈101 of 219
- `1-34` | compress | `# Static audit and validation of an MGPModel. Reads the event table only; never calls a hazard closure.` | 1,2,3,7
- `44-48` | compress | `Static facts about one Event; deme names are resolved from from/into indices.` | 5
- `69-71` | compress | `Plain printable value; show renders a table.` | 1
- `89-93` | compress | `Never calls a hazard or throws; see validate_model for checks.` | 2,4
- `142-150` | delete | — | 1,2,5
- `167-170` | delete | — | 1,2,4

### mgp.jl — trim — ≈81 of 125
- `1-19` | compress | `# Model layer: an MGPModel is a table of Events.` / `# Source: King, Lin & Ionides (StructuredMGPs.pdf); equation numbers refer to it.` | 3,5,7 ("STATUS: compiler-checked")
- `23-31` | compress | `# Pure event types (KLI §3.2). DEATH has no coloring op; SAMPLE is chop χ; NEUTRAL leaves lineages unchanged.` | 5
- `43-44` | compress | `r: production rᵘ, a per-mark constant (KLI §3.3.1), indexed by deme.` | 5
- `52-57` | compress | `s and ℓ are not stored; the filter computes them from the pruned genealogy (KLI §3.4.2). The coloring operator is derived at runtime from type and Qᵤ (Eq. 38).` | 5,7
- `84-90` | compress | `# SEIR written out by hand. Demes are (E, I).` | 4,5
- `92, 95, 98, 101, 104, 107` | delete | — | 6
- `112-120` | compress | `# @mgp (mgp_macro.jl) builds an MGPModel from declarative syntax.` | 3,5 ("VERIFIED:")

### mgp_si2r.jl — trim — ≈31 of 56
- `1-6` | compress | `# SI2R superspreading model. Demes: I_L (low-rate spreader), I_H (super-spreader).` | 4
- `7-18` | delete | — | 6 (event table restates the `@event` lines)
- `26-27`, `30-31`, `34`, `37`, `40`, `43`, `46`, `49`, `52` | delete | — | 6
- `56` | delete | — | 7 (orphan sentence; file has no trailing newline)

### mgp_sir.jl — trim — ≈13 of 23
- `1-13` | compress | `# SIR: one deme (I). Sampling is non-destructive (pop=()), so the host stays infectious.` / `# No hand-coded filter; src/simulate.jl runs it from the event table alone.` | 6

### mgp_macro.jl — acceptable
- `176-177` | optional compress | `Build an MGPModel from a declarative event specification.` | 5

### mgp_mers.jl — acceptable
- `1-2` | optional compress | `# Demes are (I_c, I_h): camel, human.` | 4

---

## C. Filters and core `src/`

### examples/mers_funs.jl — rewrite — ≈105 of 447
- `1-30` | compress | `## Shared by SoftMERS and HardMERS: rates, singular part, filter_pomp. Each module supplies its own regular_part!.` / `## Samples carry the host deme. Coalescences are not deme-fixed, so knowledge! fixes the guide only at Samples.` | 4,5,7
- `35-41` | compress | `Fix the guide at Sample tips, where the host deme is known; leave Nodes free. Returns true if v was filled.` | 4,6
- `71-101` | compress | `## pi[k] is the share of the population rate proposed at slot k; alpha[k] = rate*pi[k].` / `## Unclaimed rate is returned as decay; a firing slot is charged -log(pi[k]). no_move_share sets the no-move share.` | 1,4,7
- `71, 115` | delete | — (banner dash lines; keep table 102-114) | 5
- `143-144` | compress | `## Leftover rate is returned as decay (zero when pi[3]+pi[4] = pi[5]+pi[6] = 1).` | 4
- `202-204` | delete | — | 4,5
- `376-377` | delete | — | 6

### examples/mers_soft.jl — rewrite — ≈52 of 150
- `1-27` | compress | `SoftMERS: filter for the two-host MERS model with soft proposals. Population-event rates equal the model's; the guide only picks among tracked branches.` | 4,5,7
- `42-48` | compress | `## Soft: onC = I_c-offC, so pi[3]+pi[4] = 1 and transmission adds no decay. The guide enters only through choose_branch.` | 4,7
- `49-63` | delete | — | 1,7 ("former on-disk version" proof)
- `79` | delete | — | 6

### examples/mers_hard.jl — rewrite — ≈40 of 137
- `1-24` | compress | `HardMERS: filter for the two-host MERS model with hard proposals. Branch intensities use the raw relative hazards; the shortfall or excess is returned as decay.` | 4,5,7
- `39-45` | compress | `## Hard: onC = sum_relhaz(...), so pi[3]+pi[4] need not be 1; the leftover goes to decay. relhaz! must run before event_rates!.` | 7
- `46-52` | delete | — | 1,4

### examples/mtbd_funs.jl — trim — ≈46 of 79
- `13-16` | compress | `## There are no susceptibles: every rate is per infected host.` | 5,7
- `70-75` | compress | `## Default parameters: MLE for the 274-tip MERS tree (stem 1.0, sampling proportions 0.05).` | 2,3,4

### examples/mtbd_naive.jl — trim — ≈47 of 312
- `10-15` | compress | `Naive proposals: colour changes ignore downstream tip types, so the filter degenerates on large trees. See GuidedMTBD.` | 3,4
- `124-137` | compress | `## Rates of unobserved events and the naive split pi. A cross-type birth proposes "no move" with share 1/(ell+1), so it stays possible when ell = I.` / `## A migration carries its lineage, so "no move" needs an untracked migrant.` | 1,4,7

### examples/mtbd_guided.jl — trim
- `20-22` | delete | — | 3,5

### examples/mers_naive.jl — trim — ≈23 of 334
- `152-162` | compress | `## Cross-deme births: no move (k = 3, 5) leaves the new host untracked; move (k = 4, 6) passes a tracked lineage. no_move_share keeps no move positive when every source host is tracked.` | 1,4,7
- `188-191` | compress | `## Loop on t < tf so the final partial interval is charged decay.` | 1,5
- `264-265` | delete "Parameter names and event order follow R phylopomp's MERS model." | 4

### examples/mers_guided.jl — trim
- `28-30` | compress | `Fix the guide at Sample tips, where the host deme is known.` | 4,6
- `340, 350` | delete trailing `# preboost` | 5
- `450-451` | delete | — (example names a binding not in this module) | 4,6

### guide.jl — trim
- `47`, `80` | delete | — (`## FIXME: inelegant to store redundant information`) | 5
- `235-244` | compress | `Target probability that no tracked lineage passes to the new host when a deme-i host infects a deme-j host. Positive even when every deme-i host is tracked.` | 7
- `262-270` | compress | `Pass the target factors as w0 and wmove. Unlike choose_branch, w0 > 0 even when every host in deme i is tracked. For a migration, set w0 = 0 in that case.` | 1,7

### simulate.jl — trim (small)
- `1-2` | compress | `# Forward simulation of any MGPModel by the Gillespie algorithm.` | 4
- `11-13` | compress | `# Deme is stored only when the caller passes demeset and samplemap; then only Sample nodes carry a deme.` | 4
- `33-35` | compress | `## Lineages are addressed by position; order within a deme is irrelevant.` | 5
- `94-98` | compress | `No extra node is created at the same instant, so Newick round trips are lossless.` | 4,7
- `168-171` | compress | `Pass 1 marks Samples and their ancestors. Pass 2 splices out degree-1 Nodes, parents first.` | 7
- `182`, `189` | delete | — | 6
- `332-334` | compress | `## Δ as tuples in keys(x0) order, and each deme's position in x0.` | 5,6

### genealogy.jl — trim
- `10` | delete | — (FIXME) | 5
- `109`, `142` | delete | — | 6
- `114-117` | compress | `## Sort by time; ties go ancestor before descendant (depth from root), then by name.` | 1,7

### rcateg.jl — trim
- `38-40` | delete | — | 5,7

### parse.jl — trim
- `104-105` | delete | — | 6
- `136-137`, `165-166` | compress | `Leaves G needing repair!; see [`repair!`](@ref).` | 1

### fsmarkov.jl — trim
- `122` | delete | — | 6

Acceptable: `seir_funs.jl`, `seir_naive.jl`, `seir_soft.jl`, `seir_hard.jl`, `seir_guided.jl`, `seir_trees.jl`, `mers_tree.jl`, `Examples.jl`, `PhyloPOMP.jl`, `coloring.jl`, `demes.jl`, `indicator.jl`, `newick.jl`, `cblv.jl`. Terminology note: `seir_guided.jl:1-9` says "soft"; GuidedMERS and GuidedMTBD say "semisoft".

---

## D. `test/`

### kli_kingman_moran_test.jl — rewrite — ≈148 of 364
- `1-74` | compress | `"""Kingman/Moran special case (KLI 4.3.1, Eqs. 17-20). TCC and THH put both fork slots in one deme, so the fork weight is 1/C(N,2) and the any-pair rate is alpha*C(ell,2)/C(N,2). Expected values come from brute-force pair counts, not the compiler."""` | 1,2,5,7
- `79` | compress | `@info h1("Kingman coalescent / Moran special case (KLI 4.3.1)")` | 2
- `89-92` | compress | `## Expected values, computed without the compiler.` | 5
- `94-99` | compress | `"""Count unordered pairs of 1:N by double loop."""` | 5,6
- `108-113` | compress | `"""Count pairs of 1:N lying inside the first ell elements."""` | 5
- `122-132` | compress | `"""C(ell,2)/C(N,2) as an exact Rational; 0 when N < 2."""` | 5,6
- `139` | compress | testset name `"Kingman/Moran special case"` (drop "M12b") | 2
- `141-144`, `153`, `160-163`, `301` | delete | — | 3,5,6
- `170-172` | compress | `## Anchor values: 1/C(5,2), 1/C(2,2), 1/C(100,2).` | 2,3,4
- `189-192` | compress | `## The H-deme factor is 1, so phi_fork must not depend on the other deme.` | 3
- `226-234` | compress | `## Same weight via full_transitions; Identity + ell*Inline + C(ell,2)*Fork sums to 1.` | 2,5,7
- `259-266` | compress | `## Any-pair rate from the compiled hazard and fork weight equals alpha*C(ell,2)/C(N,2).` | 2,5
- `267` | compress | `@info h2("Part 3: compiled any-pair rate == alpha * C(ell,2)/C(N,2)")` | 5
- `279-280` | compress | `## Positive random parameters keep alpha nonzero.` | 5
- `290, 291, 294, 297` | compress | keep `## alpha_TCC matches the declared MERS rate.`; drop "REAL ... (M01/M02)" | 2,5
- `338-346` | compress | `## Boundary states: fewer than 2 individuals gives 0; ell = N = 2 gives 1.` | 5,7
- `356-357` | compress | `## ell_C == I_C == 2: the only pair must coalesce.` | 5

### kli_full_transitions_test.jl — rewrite — ≈140 of 369
- `1-18` | compress | `"""Tests for full_transitions: each saturation is classified Identity, InlineSameDeme, CrossDeme or Fork with an exact Rational phi. MERS cases follow the worked table in mers_filter_suite.tex; SEIR cases use the same binomial ratio."""` | 2,4
- `23`, `41` | compress | drop "(M03)" | 2
- `45-47` | compress | `## TCC, I_C=5, ell_C=2. The H deme is irrelevant (r_H=0).` | 2,4
- `55-59` | compress | drop "(tex:479 ...)" | 4
- `62-65` | compress | `## Inline slot is in deme C; both fork slots are in C.` | 4
- `75-80` | compress | `## THH needs ell_H >= 2 so the fork saturation is enumerated.` | 2,4
- `85-87` | compress | `## I_H=6, ell_H=2: phi = 2/5, 4/15, 1/15; the ell-weighted sum is 1.` | 7
- `108-110` | compress | `## C slot is the continuing parent; H slot is the new child, whose ancestral deme is C.` | 4
- `128-131` | compress | `## (0,1) must be CrossDeme (C to H) and (1,0) InlineSameDeme, never the reverse.` | 2,5
- `147-148` | compress | `## Mirror of THC with C and H swapped.` | 4
- `177-179` | compress | `@info h2("SEIR infection (r=(1,1) at parent deme I): same shape as THC")` | 3
- `181-198` | compress | `## infection: one slot is the new E child (ancestral deme I), the other the continuing parent in I.` / `## n_E=6, ell_E=2, n_I=5, ell_I=2; phi = product of per-deme ratios.` | 4,7
- `215-237` | delete | — | 2,3,4,6,7 (ends in vacuous `@test true`)
- `240-241` | compress | `@info h2("SEIR progression (r=(0,1), MIGRATION)")` | 2,3
- `243-247` | compress | `## progression: single slot at the target deme I, so s_I=1 is CrossDeme and s_I=0 is Identity.` | 2,3
- `252-256` | compress | `## n_I=5, ell_I=2: phi(s_I=0)=3/5, phi(s_I=1)=1/5, with no ell dependence in the second.` | 4
- `269-290` | compress | `## phi(s_I=1) = 1/n_I, with no ell dependence (r=(0,1) gives C(n-ell,0)/C(n,1)).` | 1,2,3,4,7
- `340-347` | compress | `## Events are rejected by type (DEATH/SAMPLE/NEUTRAL), not by r. SEIR sampling has r=(0,1) and is still rejected.` | 4

### kli_properties_test.jl — rewrite — ≈96 of 309
- `1-23` | compress | `"""Randomized property tests over full_transitions, reduce_event_indicator and kli_binomial_ratio for every BIRTH/MIGRATION event of SEIR and MERS. Hundreds of random (ell, n) states per property."""` | 1,2,7
- `28`, `64` | compress | drop "(M12)" | 2
- `41-47` | compress | `## All BIRTH/MIGRATION events, collected by type so new events are covered.` | 5
- `53-57` | compress | `## Random state with ell_d <= n_d by construction.` | 6
- `66, 113, 144, 182, 252` | delete | — (banners) | 5
- `77-78` | compress | `## Generator check.` | 5
- `82-85` | compress | `## enumerate_saturations bounds s_d by min(r_d, ell_d).` | 2,3
- `127-130` | compress | `## Both sides sum the same phi_u values, partitioned differently.` | 2
- `156-159` | compress | `## phi_u is a ratio of non-negative binomials.` | 7
- `163-168` | compress | `## Phi_u is non-negative, so proposal weights are well defined.` | 2,4
- `171-173` | compress | `## Check kli_binomial_ratio directly, bypassing classification.` | 6
- `217-219` | compress | `## enumerate_saturations ignores n, so it does not error here.` | 3
- `231-234` | compress | `## s_d < 0 is never enumerated; kli_binomial_ratio must still return 0.` | 1,2,4

### kli_filter_ir_test.jl — rewrite — ≈71 of 270
- `1-14` | compress | `"""Tests for classify_filter_terms and filter_spec: reduced transitions go into RegularFlow, SingularFlow and OutflowImbalanceTerm buckets. (ell, n) instances and Phi values match kli_reduce_test.jl."""` | 2,4,6
- `19`, `31` | compress | drop "(M05)" | 2
- `83-86` | compress | `## No fork term and no decay-bucket term exist for progression.` | 4
- `115-123` | compress | `## Phi(noop)+Phi(fork) = 7/10, not 1: the sum of phi over s is not a probability.` | 4,7
- `126-128` | compress | `@info h2("MERS transmission_hc (THC): noop+cross to RegularFlow, fork to the decay bucket")` | 6
- `145-147` | delete | — | 4
- `165` | compress | `# traces to full_transitions' output` | 2
- `183-186` | compress | `@info h2("Singular BIRTH/MIGRATION branch, via a synthetic Event (no shipped model has one)")` | 2
- `215-221` | compress | `## A non-regular event has no imbalance bucket; all its outcomes are singular.` | 4
- `224-228` | compress | `@info h2("OutflowImbalanceTerm.mechanism is never :lambda")`; drop "(M06 correction)" from testset | 1,2,3,4
- `235-241` | compress | `## The mechanism tag never claims to be KLI's lambda; lambda is not implemented.` | 2,7
- `251-254` | compress | `## FilterSpec has no field named decay.` | 2,7
- `258-263` | compress | `## kli_decay is a separate code path and is still an error stub.` | 2,3,4

### kli_mers_compiled_test.jl — rewrite — ≈68 of 380
- `1-38` | compress | `"""Compiled MERS filter vs NaiveMERS: seed-matched Np=1 log-likelihoods must agree on simulated genealogies. naive_oracle_filter_pomp is NaiveMERS.filter_pomp taking gen as an argument. Genealogies use NaiveMERS.Demes so the oracle's deme comparisons see the same enum type."""` | 2,5,7
- `43`, `182` | compress | drop "(M09 Gate 5)" / "(M09)" | 2
- `56-61` | compress | `## Saturations are absent when ell is too small; phi_of returns 0 for them.` | 7
- `67-71` | compress | `## NaiveMERS.filter_pomp with gen as an argument.` | 5
- `184-186` | compress | `@info h2("TCC/THH: naive's drawless no-fork term equals Phi_id + ell*Phi_inline, not the plain noop sum")` | 3
- `191` | compress | `# camel side: ell_C=2, I_C=5 (post-event)` | 2
- `206-212` | compress | `## Other (ell, I) cases; ell_C=1 has no Fork saturation, so phi_of returns 0.` | 1,3
- `254-256` | compress | `@info h2("Decay: this MERS instance reduces to 18/5")`; testset `"decay instance"` | 2
- `267-271` | compress | `## compiled_decay computes in Float64, so compare with isapprox.` | 4
- `331-334` | compress | `@info h2("End-to-end log-likelihood over many (params, genealogy, seed) combinations; beta_hc/beta_ch nonzero")` | 2
- `371` | compress | drop "Gate 5" from the log string | 2
- `137` | `demes_sampled` defined, never used

### kli_proposal_test.jl — rewrite — ≈51 of 126
- `1-17` | compress | `"""Tests for ProposalStrategy, boost and driver: driver = alpha*pi, boost = Phi/pi. boost(Phi(cross), pi) is compared with seir_naive.jl's tracked-parent infection branch in exact Rational arithmetic."""` | 2,4
- `22`, `31` | compress | drop "(M07 Part B)" | 2
- `44`, `50-51` | compress | drop "trivial, but exact" | 5
- `64-65` | compress | drop "seir_naive.jl:150-156" | 4
- `67-82` | compress | `## n_E=6, ell_E=2, n_I=5, post-event ell_I=2; Phi(cross)=1/10. Pre-event tracked I count is 3, since the swap removes one.` | 2,4
- `88` | compress | drop "# seir_naive.jl:116, pi[2] = ellI/I" | 4
- `91-101` | compress | `## Naive correction factor = (1/pi2) * a * (1 - ellI_post/I) * (1/E), without the decay term.` | 4,7
- `102-110` | compress | `## pi_u(cross) = pi2 * (1/a) = 1/n_I; the ell factors cancel.` | 4,7
- `121` | delete | — | 2,5

### kli_si2r_test.jl — trim — ≈204 of 745
- `1-17` | compress | `"""SI2R verification: runs phi_u, Phi_u grouping, Chu-Vandermonde and decay on the SI2R superspreading model and compares with values hand-derived in si2r_model.qmd. SI2R was not used while developing the compiler."""` | 2,4,5 (absolute home path)
- `30-32`, `57-59` | delete | — | 5,6
- `39-42` | compress | `## Binomial ratio from the QMD, in exact Rational{Int}.` | 5
- `65, 72, 82, 89, 96, 103, 110, 116, 123` | delete | — | 6
- `131-138` | compress | `## Hazards for all nine marks match the QMD alpha_u column. A wrong hazard, e.g. TH using I_L, would pass every other test here.` | 5,7
- `161-167` | compress | `## kli_binomial_ratio equals the QMD formula for every event, swept over random states.` | 5
- `189, 201, 213, 225, 237` | delete | — | 4
- `246`, `257`, `290`, `308`, `334`, `347` | compress | drop "QMD lines N" | 4
- `271-279` | compress | `## Classification: TL gives Identity/Inline/Fork; TH adds CrossDeme; L and H give Identity/CrossDeme.` | 5,7
- `314-316`, `320-321`, `326-327`, `339-340`, `352-353` | compress | drop "QMD line N: swap(..)" | 4
- `358-373` | compress | `## Identity and InlineSameDeme form one :noop group. Phi_noop is the unweighted sum; the QMD closed form (Eq 132/135) weights Inline by ell, so test phi_id + ell*phi_inl.` | 2,4,5,7
- `413-417` | compress | `## InlineSameDeme is absent when ell=0; treat it as phi=0.` | 7
- `493-501` | compress | `## Per event: Phi equals the member phi sum, members match the QMD, and no mass is lost.` | 5 (lists steps (a)-(d); code labels the last (c))
- `549-556` | compress | `## lambda = psi*(I_L+I_H) + gamma*I_L*1{I_L<=ell_L} + gamma*I_H*1{I_H<=ell_H}.` | 5
- `614-620` | compress | `## TL fork weight is 1/C(I_L,2), independent of ell_L (same Kingman/Moran identity as TCC).` | 2,5
- `644-645` | delete | — | 3
- `651-660` | compress | `## TH fork is an ordered pair, 1/(I_L*I_H), unlike TL's unordered 1/C(I_L,2). A compiler that pattern-matched MERS TCC/THH could get this wrong.` | 4,5,7
- `686-693` | compress | `## Simplified closed forms from the QMD, independent of qmd_binomial_ratio's product form.` | 5,7
- `707-708`, `721-724` | compress | drop "boost line N" | 4

### kli_phi_test.jl — trim — ≈70 of 252
- `1-21` | compress | `"""Reference-value tests for production_slots, enumerate_saturations and kli_binomial_ratio (phi_u). Values are hand-computed in exact rationals; MERS cases follow the worked table in mers_filter_suite.tex."""` | 2,4
- `34` | compress | drop "(M02)" | 2
- `90-94`, `103-105`, `124-128`, `138-142` | compress | drop tex line and table-row cites; keep formulas | 4
- `221-225` | compress | `## s=-1 gives r-s=2, which the b<0 rule alone would not zero; kli_binomial_ratio guards s < 0 explicitly.` | 1,3
- `231-234` | compress | `## r=(0,0): s=(0,0) is forced and phi=1 for removal and sampling alike.` | 2,4

### mgpaudit_test.jl — trim — ≈43 of 274
- `1-19` | compress | `"""Tests for explain and mgpaudit/@mgpaudit. explain walks FilterTerm to ReducedTransition to KLITransition to Event; mgpaudit gives a derivation per BIRTH/MIGRATION event and an out-of-scope note for the rest."""` | 2,4,6
- `24`, `43` | compress | drop "(M06)" | 2
- `104-106` | compress | `@info h2("explain(::OutflowImbalanceTerm): event name, kind and mechanism tag, MERS transmission_cc")` | 2
- `124-126` | compress | `## explain never labels this mechanism :lambda.` | 6
- `166-167` | compress | `## progression has no OutflowImbalanceTerm.` | 2
- `171-172` | compress | `## The printed report names every event.` | 5
- `265-266` | compress | `@info h2("default_audit_state per-model values")` | 2

### kli_reduce_test.jl — trim — ≈54 of 234
- `1-14` | compress | `"""Tests for reduce_event_indicator and reduced_transitions: full transitions are grouped by color-only outcome and phi is summed into Phi_u. Instances and phi values are those of kli_full_transitions_test.jl."""` | 2,4
- `19`, `31` | compress | drop "(M04)" | 2
- `89-91` | compress | `@info h2("THC: 4 full to 3 reduced; (1,0) Inline collapses into noop, (0,1) Cross does not")` | 5
- `109-111` | compress | `## Collapsing follows deme match, not saturation index.` | 5

### kli_decay_test.jl — trim — ≈39 of 198
- `1-15` | compress | `"""Tests for decay_contribution and total_decay (KLI Eq. 47 / App. B Eq. B2) for DEATH and SAMPLE events, in exact Rational arithmetic. Expected values are hand-computed."""` | 2,4,5
- `20`, `29` | compress | drop "(M07 Part A)" | 2
- `51-56` | compress | `## lambda = chi_C I_C + chi_H I_H + gamma_C I_C 1{I_C<=ell_C} + gamma_H I_H 1{I_H<=ell_H} = 13/5.` | 4
- `117` | compress | drop "# M03's finding" | 2
- `170-173` | compress | `## total_decay skips BIRTH/MIGRATION/NEUTRAL: compare with a manual sum over DEATH/SAMPLE events.` | 5

### kli_seir_compiled_test.jl — trim — ≈28 of 162
- `1-28` | compress | `"""Compiled SEIR filter vs NaiveSEIR: seed-matched Np=1 log-likelihoods agree to atol 1e-6, rtol 1e-8 over many parameter/genealogy/seed combinations. compiled_regular_part! makes the same random draws in the same order as NaiveSEIR.regular_part!."""` | 1,2,5
- `33`, `78` | compress | drop "(M08 Gate 5)" / "(M08)" | 2
- `122-123` | compress | `@info h2("End-to-end log-likelihood, seed-matched Np=1 filters")` | 2
- `153` | compress | drop "Gate 5" from the log string | 2

### population_ir_test.jl — trim — ≈25 of 139
- `1-17` | compress | `"""Structural check of the Population IR (Event/MGPModel) for SEIR and MERS: each event carries delta, hazard, production vector and from/into wiring. audit_model(SEIR) matches SEIR_REFERENCE; both models pass validate_model."""` | 2,4
- `51-52` | compress | `## Event stores no runtime event time.` | 5
- `118` | compress | `@info h2("validate_model catches a broken model")` | 5

### mers_simulate.jl — trim — ≈45 of 227
- `31-35` | compress | `## Small N and R0 near 1 make cross-species, few-sample trees common within the retry budget.` | 2,4
- `72-74` | compress | `## Filter parameters match the simulation; mers_rinit rescales each species by N/(S0+I0).` | 4
- `106-111` | compress | `## Require both species sampled (exercises transmission_hc/ch) and at most 12 samples (keeps pfilter cheap).` | 7
- `122-124` | compress | `## MERS sampling is destructive: Sample nodes have no children.` | 4
- `172-176` | compress | `## demeset/samplemap only label nodes and use no randomness, so the same seed replays g; the newick equality checks that.` | 4,5
- `193-197` | compress | `## Independent realization. GuidedMERS.check returns nothing or throws.` | 5,6
- `198` repeats the `@info` at `192`

### mers_soft.jl (test) — trim — ≈28 of 121
- `54-64` | compress | `## Shared defaults (mers_funs.jl). N_c = N_h = 10000 makes a 3-tip tree hard to filter, so assert only that the estimate over 5 replicates is finite.` | 1,3
- `72-83` | compress | `## Soft, hard and guided are proposals for one likelihood, so their pfilter estimates agree within Monte Carlo error (z < 3). 15 replicates at Np=2000 keep the heavy-tailed estimator stable.` | 1,3,7

### mtbd.jl (test) — trim — ≈38 of 121
- `54-61` | compress | `## Any guide gives an unbiased filter; it only changes variance. Conductance 3 keeps the guide informative and the replicate spread small.` | 1,7
- `98-111` | compress | `## Both filters return -Inf on the 274-tip MERS tree (known failure). Exact log likelihood at mers_mle with a 1-year stem: +381.8742.` | 1,3,7

### seir_macro_equivalence.jl — trim — ≈16 of 146
- `1-13` | compress | `"""Equivalence of the @mgp SEIR event table with NaiveSEIR, model layer only: compartments, demes, hazards, increments, event types, regular driver rates, decay (N = pop). Full-filter equivalence is not tested."""` | 5,7
- `135-136` | compress | `## With the culling event in the table, equivalence holds for chi > 0.` | 1

### seir_simulate.jl — trim
- `68-69` | compress | `## NaiveSEIR draws psi against chi at each sample.` | 4
- `87` | compress | `## With chi = 0 the culling hazard is 0.` | 1

### sir_simulate.jl — trim
- `272`, `276` | compress | drop "(M3 bug)" and "(M4 fork bug)" | 1,2

### rcateg.jl (test) — trim
- `58-62` | compress | `## A pinned RNG hits the measure-zero branches: u -> 1 overshoots the last weight by rounding; u == 0 lands on a zero-weight first category.` | 1

### seir_hard.jl (test)
- `17` testset named "SEIR model with guided proposals" but tests the hard kernel.

Acceptable: `simulate_checks.jl`, `cblv.jl`, `mers_guided.jl`, `mers_hard.jl`, `mers_naive.jl`, `seir_naive.jl`, `seir_soft.jl`, `seir_guided.jl`, `fsmarkov.jl`, `guide.jl`, `newick.jl`, `parse.jl`, `runtests.jl`.
