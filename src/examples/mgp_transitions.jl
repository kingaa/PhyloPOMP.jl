# mgp_transitions.jl
# =============================================================================
# M03: Full KLI Compatibility Lowering.
#
#   Generic derivation of the FULL (uncollapsed) set of KLI-compatible
#   coloring transitions Y -> Y' for a single production mark u, given its
#   current per-deme lineage counts ell and total counts n. This consumes
#   M02's `production_slots` / `enumerate_saturations` / `kli_binomial_ratio`
#   (src/examples/mgp_phi.jl) -- it does NOT reimplement the saturation
#   enumeration or the binomial-ratio math, only Step A/D of the tex's
#   "Construction of Compatibility Terms" (src/examples/mers_filter_suite.tex,
#   lines 434-467): classifying each feasible saturation `s` by how many
#   production slots it fills and whether the filled slot(s)' ancestral deme
#   matches their post-event deme.
#
#   Explicitly OUT of scope for this milestone (M03):
#     - Q_u itself as a genealogy-state-dependent 0/1 gate (e.g. whether a
#       Fork saturation is actually reachable at a REGULAR event, as opposed
#       to only at a SINGULAR/observed branch-point time) -- `full_transitions`
#       returns the full, ungated compatibility structure; gating by Q_u/the
#       decay machinery is later-milestone (driver/proposal) territory.
#     - DEATH/SAMPLE/NEUTRAL events: their compatibility is decay/rate-driven
#       (mers_filter_suite.tex, Sec. "RC and RH"/"SC and SH", lines 511-539),
#       not saturation-driven, regardless of whether r_u happens to be all
#       zero. This is true for MERS's removal_c/h, sampling_c/h (built from
#       `sample_remove`, r=(0,0)) and SEIR's recovery/waning, but NOT
#       uniformly true of every SAMPLE event: SEIR's singular `sampling`
#       (built from `move=sample(I)`, mgp.jl:105-106) has r=(0,1) -- a real
#       production slot for the inserted leaf -- yet is still out of scope
#       here, because the scope boundary is the event's TYPE
#       (BIRTH/MIGRATION only), not whether its r happens to be zero.
#       `full_transitions` throws on all DEATH/SAMPLE/NEUTRAL events;
#       `ChopTransition` is defined as a forward-compatible placeholder type
#       only (see its docstring).
#     - Reduction/marginalization over the event-indicator m (Phi_u =
#       sum_m phi_u for full transitions collapsing to the same reduced
#       d->d'), i.e. KLI's explicit-reduction step -- that is M04.
#     - Any wiring into mgp_filter.jl's stubs (kli_select, kli_decay,
#       apply_move!, singular_update!) -- this file produces IR/derivation
#       only, not an executable filter.
#
# Primary sources: King, Lin & Ionides, "Exact phylodynamic likelihood via
#   structured Markov genealogy processes" (StructuredMGPs.pdf), Sec. 3.4-3.5
#   / Eq. 9 (as already restated/vetted by M02), and this repo's own worked
#   derivation `mers_filter_suite.tex`, Secs. "Construction of Compatibility
#   Terms" (Steps A-D, lines 434-467) and "Event-Specific Derivations"
#   (TCC/THH/THC/TCH, lines 468-509), whose Complete Reference Table
#   (lines 541-577) is the term-by-term oracle for `test/kli_full_transitions_test.jl`.
# =============================================================================

export KLITransition, IdentityTransition, InlineSameDemeTransition,
       CrossDemeTransition, ForkTransition, ChopTransition, full_transitions

"""
    KLITransition

Abstract supertype for the classified outcomes of a single production mark
`u`'s saturation `s`: a full (uncollapsed) `Y -> Y'` coloring transition,
per `mers_filter_suite.tex`'s Step D (lines 459-467). Every concrete
subtype except `ChopTransition` carries, at minimum:

- `event` : the source `Event` (`u`).
- `s`     : the saturation `Vector{Int}` that produced this transition
            (one element of `enumerate_saturations(event, ell)`).
- `phi`   : `kli_binomial_ratio(event, s, ell, n; Q)`, the production-slot
            compatibility ratio for this saturation (exact `Rational{Int}`).

Concrete subtypes distinguish the four KLI-meaningful classes the project
mandate (Sec. 2.2) requires kept separate -- in particular
`IdentityTransition` and `InlineSameDemeTransition` are DISTINCT even though
both leave the reduced coloring unchanged (`y' = y`): only the former is
"no slot occupied", the latter is "exactly one slot occupied, but its
ancestral deme equals its post-event deme", which is still a distinct full
KLI transition because the event-indicator component `m` changes
(`mers_filter_suite.tex:462-463`, "on KLI's full coloring it is still a
distinct inline-event transition because the event-indicator component
changes").
"""
abstract type KLITransition end

"""
    IdentityTransition(event, s, phi)

No production slot is occupied by a tracked lineage (`s` is all zeros):
`Y' = Y`, including at the full-coloring level (`mers_filter_suite.tex:465-466`,
"No slot occupied: `y'=y`").
"""
struct IdentityTransition <: KLITransition
    event :: Event
    s     :: Vector{Int}
    phi   :: Rational{Int}
end

"""
    InlineSameDemeTransition(event, s, phi, deme)

Exactly one production slot is occupied by a tracked lineage `b`, and that
slot's ancestral deme equals its post-event deme (`deme`). On the reduced
color-only state this is a no-op (`y' = y`), but it is a DISTINCT full-KLI
transition from `IdentityTransition` -- see `KLITransition`'s docstring and
`mers_filter_suite.tex:460-463` (`sigma^b_{d,d} y = y` when
`d_pre = d_post`). `deme` is the (single) deme index common to both the
ancestral deme and the occupied slot's post-event deme.
"""
struct InlineSameDemeTransition <: KLITransition
    event :: Event
    s     :: Vector{Int}
    phi   :: Rational{Int}
    deme  :: Int
end

"""
    CrossDemeTransition(event, s, phi, ancestral_deme, target_deme)

Exactly one production slot is occupied by a tracked lineage `b`, and that
slot's ancestral deme differs from its post-event deme: `y' = swap!(y,
ancestral_deme, target_deme, b)` (`mers_filter_suite.tex:460-461`,
`sigma^b_{d_pre,d_post} y`, `d_pre != d_post`). `ancestral_deme` is the
lineage's deme immediately before the event (`event.from`); `target_deme`
is the deme of the one occupied slot.
"""
struct CrossDemeTransition <: KLITransition
    event          :: Event
    s              :: Vector{Int}
    phi            :: Rational{Int}
    ancestral_deme :: Int
    target_deme    :: Int
end

"""
    ForkTransition(event, s, phi, ancestral_deme, slot_demes)

Two (or more) production slots are occupied by tracked lineages `b, b'`
sharing a common ancestor deme (a branch point): `y' = fork!(y,
ancestral_deme, b, slot_demes, [b, b'])` (`mers_filter_suite.tex:464-465`,
`kappa^{bb'}_{d_anc,d_b,d_b'} y`). `slot_demes` lists the post-event deme
of every occupied slot, WITH multiplicity (e.g. TCC's `s=(2,0)` gives
`slot_demes = [1, 1]`, both slots landing in the camel deme; THC's
`s=(1,1)` gives `slot_demes = [1, 2]`).
"""
struct ForkTransition <: KLITransition
    event          :: Event
    s              :: Vector{Int}
    phi            :: Rational{Int}
    ancestral_deme :: Int
    slot_demes     :: Vector{Int}
end

"""
    ChopTransition(event)

Placeholder `KLITransition` subtype ONLY -- its derivation is a documented
STUB, not implemented in M03. `DEATH` events (and MERS's `sample_remove`-
based `SAMPLE` events) have `r_u = (0,...,0)`, so `enumerate_saturations`
always forces `s=(0,...,0)` for them and the saturation/phi_u machinery is
degenerate (`phi_u = 1`, the empty product, `mers_filter_suite.tex:530`,
"KLI Eq.(9) gives the empty-product value ... = 1"). SEIR's `sampling`
event is the one exception with a genuinely nonzero `r_u=(0,1)` (built from
`move=sample(I)`, mgp.jl:105-106, inserting a leaf slot) -- but it is
STILL out of scope here, along with every other `SAMPLE`/`DEATH`/`NEUTRAL`
event, because compatibility for this whole event-type class is governed by
DECAY/RATE semantics, not the saturation/phi_u machinery -- e.g. removal's
outflow condition `I_d - 1 >= ell_d`
and the sub-threshold decay term `gamma_d I_d * 1{I_d <= ell_d}`
(`mers_filter_suite.tex:511-517`), and sampling's split between an observed-
leaf singular update and a `lambda`-contributing background rate
(`mers_filter_suite.tex:519-539`) -- driver/decay derivations that are
explicitly OUT of scope for M03 per the milestone spec (deferred to the
later driver/boost/decay milestone). This type exists purely so the
transition-IR vocabulary is complete and forward-compatible; it carries no
`s`/`phi` fields and is never constructed by `full_transitions` (which
throws `ArgumentError` on DEATH/SAMPLE/NEUTRAL events instead of returning
`ChopTransition`s).
"""
struct ChopTransition <: KLITransition
    event :: Event
end

"""
    full_transitions(event::Event, ell, n; Q = 1) -> Vector{KLITransition}

Derive the full (uncollapsed) set of KLI-compatible `Y -> Y'` coloring
transitions for production mark `event`, given current per-deme tracked
lineage counts `ell` and total per-deme population counts `n`. Only
`event.type in (BIRTH, MIGRATION)` is supported (see `ChopTransition`'s
docstring for why `DEATH`/`SAMPLE`/`NEUTRAL` are out of scope); calling this
on any other event type throws `ArgumentError`.

Algorithm (`mers_filter_suite.tex` Steps A-D, lines 434-467):
1. `S = enumerate_saturations(event, ell)` (M02) -- every feasible
   saturation `s`.
2. For each `s`, `phi = kli_binomial_ratio(event, s, ell, n; Q=Q)` (M02).
3. Classify `s` by `sum(s)` (how many production slots it fills) and, for
   the one/two-slot cases, whether the filled slot(s)' post-event deme(s)
   equal the ancestral deme -- which, per Step A (tex lines 438-441), is
   ALWAYS `event.from` (the single parent/source deme), for every occupied
   slot, regardless of which deme that slot's post-event deme is:
     - `sum(s) == 0`             -> `IdentityTransition`
     - `sum(s) == 1`, slot deme `d == event.from` -> `InlineSameDemeTransition`
     - `sum(s) == 1`, slot deme `d != event.from` -> `CrossDemeTransition`
     - `sum(s) >= 2`                              -> `ForkTransition`
4. Returns one `KLITransition` per feasible saturation, in the same order
   as `enumerate_saturations` produced them.

This classification is generic over `event.r`/`event.from` -- it is not a
per-model, per-event lookup table. It is verified, term-by-term, against
`mers_filter_suite.tex`'s TCC/THH/THC/TCH derivations and Complete
Reference Table in `test/kli_full_transitions_test.jl`.
"""
function full_transitions(event::Event, ℓ::AbstractVector{<:Integer},
                           n::AbstractVector{<:Integer}; Q::Real = 1)
    event.type in (BIRTH, MIGRATION) ||
        throw(ArgumentError("full_transitions: event `$(event.name)` has type " *
                             "$(event.type), but only BIRTH/MIGRATION events are " *
                             "derived by this milestone -- DEATH/SAMPLE/NEUTRAL " *
                             "events have r_u=(0,...,0) always and their " *
                             "compatibility is decay/rate-driven, not saturation-" *
                             "driven (see ChopTransition's docstring and " *
                             "mers_filter_suite.tex lines 511-539)."))
    event.from == 0 &&
        throw(ArgumentError("full_transitions: event `$(event.name)` has from=0 " *
                             "(no ancestral/source deme) -- cannot derive an " *
                             "ancestral deme for its production slots."))

    r = production_slots(event)
    ancestral = event.from
    S = enumerate_saturations(r, ℓ)

    transitions = Vector{KLITransition}(undef, length(S))
    for (i, s) in enumerate(S)
        φ = kli_binomial_ratio(r, s, ℓ, n; Q = Q)
        total = sum(s)
        transitions[i] = if total == 0
            IdentityTransition(event, s, φ)
        elseif total == 1
            d = findfirst(!=(0), s)
            if d == ancestral
                InlineSameDemeTransition(event, s, φ, d)
            else
                CrossDemeTransition(event, s, φ, ancestral, d)
            end
        else
            slot_demes = Int[]
            for d in eachindex(s)
                for _ in 1:s[d]
                    push!(slot_demes, d)
                end
            end
            ForkTransition(event, s, φ, ancestral, slot_demes)
        end
    end
    transitions
end
