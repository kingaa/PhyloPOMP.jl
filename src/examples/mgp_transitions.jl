# Full KLI transitions for BIRTH and MIGRATION events.
# DEATH, SAMPLE, NEUTRAL events have no saturation structure (see mgp_decay.jl).

export KLITransition, IdentityTransition, InlineSameDemeTransition,
       CrossDemeTransition, ForkTransition, ChopTransition, full_transitions

"""
    KLITransition

Abstract supertype for the classified outcome of one saturation `s` of an event.
Every concrete subtype except `ChopTransition` has these fields:

- `event` : the source `Event`.
- `s`     : the saturation (`Vector{Int}`).
- `phi`   : the compatibility ratio (exact `Rational{Int}`).

`IdentityTransition` and `InlineSameDemeTransition` both leave the reduced
coloring unchanged (`Y' = Y`). They stay distinct because the event indicator
differs.
"""
abstract type KLITransition end

"""
    IdentityTransition(event, s, phi)

No production slot is occupied by a tracked lineage (`s` is all zeros):
`Y' = Y`.
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
transition from `IdentityTransition` -- see `KLITransition`'s docstring.
`deme` is the (single) deme index common to both the
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
ancestral_deme, target_deme, b)`. `ancestral_deme` is the
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
ancestral_deme, b, slot_demes, [b, b'])`. `slot_demes` lists the post-event deme
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

Placeholder subtype. Never constructed: `full_transitions` throws for DEATH,
SAMPLE and NEUTRAL events.
"""
struct ChopTransition <: KLITransition
    event :: Event
end

"""
    full_transitions(event::Event, ell, n; Q = 1) -> Vector{KLITransition}

Derive the full (uncollapsed) KLI transitions `Y -> Y'` of `event` at per-deme
tracked counts `ell` and population counts `n`. Returns one `KLITransition` per
feasible saturation, in the order of `enumerate_saturations`. Only BIRTH and
MIGRATION events are supported. Throws `ArgumentError` for other event types
and for `event.from == 0`.

Classify each saturation `s` by `sum(s)`:
- 0: `IdentityTransition`
- 1: `InlineSameDemeTransition` if the slot deme is `event.from`, else
  `CrossDemeTransition`
- 2 or more: `ForkTransition`
"""
function full_transitions(event::Event, ℓ::AbstractVector{<:Integer},
                           n::AbstractVector{<:Integer}; Q::Real = 1)
    event.type in (BIRTH, MIGRATION) ||
        throw(ArgumentError("full_transitions: event `$(event.name)` has type " *
                             "$(event.type); only BIRTH/MIGRATION are supported"))
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
