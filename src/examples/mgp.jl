# Model layer: an MGPModel is a table of Events.
export @mgp, @event, EventType, Event, MGPModel, SEIR, MERS, SI2R, SIR, MTBD, LBDP, BDEI, BDSS

# Pure event types.
# DEATH has no coloring op; SAMPLE is chop χ; NEUTRAL leaves lineages unchanged.
@enum EventType BIRTH MIGRATION DEATH SAMPLE NEUTRAL

"""
    Event

One jump mark `u` of the population process, as it bears on the genealogy.

Fields
- `name`     : identifier.
- `Δ`        : population stoichiometry, e.g. `[:S=>-1, :E=>+1]`.
- `hazard`   : `(x, θ) -> Float64`, the population hazard αᵤ(t,x).
- `r`        : production rᵘ, a per-mark constant, indexed by deme.
- `type`     : one of the five pure `EventType`s.
- `from`     : source deme index (birth parent / migration source); 0 = none.
- `into`     : child / target deme indices.
- `regular`  : true = regular (between genealogy events); false = singular
               (at ev(Z)).
- `observed` : true = can produce a genealogy feature; false = background.

`s` and `ℓ` are not stored; the filter computes them from the pruned genealogy.
The coloring operator is derived at runtime from `type` and Qᵤ.
"""
struct Event
    name     :: Symbol
    Δ        :: Vector{Pair{Symbol,Int}}
    hazard   :: Function
    r        :: Vector{Int}
    type     :: EventType
    from     :: Int
    into     :: Vector{Int}
    regular  :: Bool
    observed :: Bool
    rate     :: Any        # the rate expression as written in `@mgp` (`nothing` for a hand-built event)
end
Event(name, Δ, hazard, r, type, from, into, regular, observed) =
    Event(name, Δ, hazard, r, type, from, into, regular, observed, nothing)

"""
    MGPModel

A model: its full compartment set, the deme subset that carries lineages
(the `Coloring` set, e.g. {E,I} ⊂ {S,E,I,R}), and its events.
"""
struct MGPModel
    name         :: Symbol
    compartments :: Vector{Symbol}
    demes        :: Vector{Symbol}
    events       :: Vector{Event}
end

# SEIR written out by hand. Demes are (E, I).
const SEIR_REFERENCE = MGPModel(:SEIR, [:S, :E, :I, :R], [:E, :I], [
    Event(:infection,   [:S=>-1, :E=>+1], (x,θ)-> θ.β * x.S * x.I / θ.N,
          [1, 1], BIRTH,     2, [1],   true,  false),
    Event(:progression, [:E=>-1, :I=>+1], (x,θ)-> θ.σ * x.E,
          [0, 1], MIGRATION, 1, [2],   true,  false),
    Event(:recovery,    [:I=>-1, :R=>+1], (x,θ)-> θ.γ * x.I,
          [0, 0], DEATH,     2, Int[], true,  false),
    Event(:waning,      [:R=>-1, :S=>+1], (x,θ)-> θ.ω * x.R,
          [0, 0], NEUTRAL,   0, Int[], true,  false),
    Event(:sampling,    Pair{Symbol,Int}[], (x,θ)-> θ.ψ * x.I,
          [0, 1], SAMPLE,    2, Int[], false, true),
    Event(:culling,     [:I=>-1], (x,θ)-> θ.χ * x.I,
          [0, 0], SAMPLE,    2, Int[], false, true),
])

# @mgp (mgp_macro.jl) builds an MGPModel from declarative syntax.

include("mgp_macro.jl")
include("mgp_mers.jl")
include("mgp_si2r.jl")
include("mgp_sir.jl")
include("mgp_mtbd.jl")
include("mgp_lbdp.jl")
include("mgp_bdei.jl")
include("mgp_bdss.jl")
