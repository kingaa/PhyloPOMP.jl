# =============================================================================
# Generic forward stochastic simulator for MGPModel-specified genealogical
# processes (see src/examples/mgp.jl for MGPModel/Event/@mgp).
#
#   Implements the Doob-Gillespie direct-method stochastic simulation
#   algorithm, building a Genealogy incrementally as birth/migration/death/
#   sample events fire. This mirrors R phylopomp's simulation kernel
#   (popul_proc_t/master_t/genealogy_t/ball_t/inventory_t, src/*.h in that
#   package) in *semantics*, not code: R needs a compile-time C++ template
#   instantiation per model because C++ rate/jump code must be compiled;
#   here the rate/jump closures already stored on each `Event` are ordinary
#   Julia functions that the compiler specializes at JIT time, so ONE
#   generic engine below drives *any* `MGPModel` with no per-model codegen.
#
#   Distinction from src/coloring.jl's `Coloring`: that type tracks only the
#   lineages a particle FILTER is actively coloring against an already-known,
#   fixed tree (size ~ number of tracked branches ℓ). Forward simulation
#   instead needs one live entry per currently-extant INDIVIDUAL in the
#   simulated population (size ~ population count) -- a different growth
#   regime -- so a separate, simpler `SimInventory` type is used here rather
#   than overloading `Coloring`.
#
#   Compartment/deme membership is, by default, treated as a LATENT variable
#   (as for SEIR): it drives the simulation internally but is not recorded
#   on the returned genealogy's nodes (`deme` stays `missing`, genealogy
#   type is `Unstructured`), matching `NaiveSEIR.filter_pomp`'s convention.
#   Some models instead OBSERVE deme membership at the moment of sampling
#   (e.g. MERS: which host species a sample came from is directly
#   observable, even though which species a lineage passed through earlier
#   is not) -- callers needing that pass `demeset`/`samplemap` to `simulate`
#   (see its docstring) so `Sample`-typed nodes carry a real deme label,
#   matching `NaiveMERS`/`SoftMERS`/etc.'s convention.
# =============================================================================

import PartiallyObservedMarkovProcesses: simulate
using Random: AbstractRNG, default_rng

export simulate

# -----------------------------------------------------------------------------
# SimInventory: per-deme bookkeeping of currently-live lineages.
# -----------------------------------------------------------------------------

"""
    SimInventory(ndeme)

Per-deme bookkeeping of currently-live lineages during forward simulation.
`inv[d]` is the `Vector{Name}` of genealogy-node names currently "open"
(i.e., representing the most recent event so far along that lineage) for
every lineage presently alive in deme `d`. Not to be confused with
[`Coloring`](@ref), which serves a different purpose in the particle-filter
code (see module header).
"""
struct SimInventory
    live::Vector{Vector{Name}}
    SimInventory(ndeme::Integer) = new([Name[] for _ ∈ Base.OneTo(ndeme)])
end

Base.getindex(inv::SimInventory, d::Integer) = inv.live[d]

remove!(inv::SimInventory, d::Integer, name::Name) = begin
    v = inv[d]
    i = findfirst(==(name), v)
    @assert !isnothing(i) "lineage $name not found in deme $d"
    deleteat!(v, i)
    nothing
end

add!(inv::SimInventory, d::Integer, name::Name) = begin
    push!(inv[d], name)
    nothing
end

# -----------------------------------------------------------------------------
# Mechanical helpers.
# -----------------------------------------------------------------------------

"""
    apply_delta(x, ev)

Apply the population stoichiometry `ev.Δ` to state `x`, returning a new
`NamedTuple`. Purely mechanical bookkeeping; no genealogy effect.
"""
apply_delta(x::NamedTuple, ev::Event) = begin
    isempty(ev.Δ) && return x
    delta = Dict(ev.Δ)
    NamedTuple{keys(x)}(map(k -> getfield(x,k) + get(delta,k,0), keys(x)))
end

"""
    push_node!(G, slate, type, parent) -> Name

Appends a new node to `G` with the given time, type, and parent, and
returns its name. Deme is left `missing` (see module header).
"""
push_node!(
    G::Genealogy{D},
    slate::Time,
    type::NodeType,
    parent::Union{Nothing,Name},
) where D = begin
    nm = Name(length(G.nodes)+1)
    push!(G.nodes, GenealNode{D.DemeSet}(nm, slate, missing, type, parent))
    nm
end

# -----------------------------------------------------------------------------
# apply_event!: the one place genealogy topology is built from an Event.
# -----------------------------------------------------------------------------

"""
    apply_event!(G, inv, ev, model, t, rng, samplemap)

Mutates `G` (the genealogy under construction) and `inv` (the live-lineage
inventory) to reflect one firing of event `ev` at time `t`. The lineage
acted upon is chosen uniformly at random (via `rng`) from the currently-live
lineages in deme `ev.from` -- the forward-simulation analogue of R's
`random_ball(i)` (`src/inventory.h`).

- `BIRTH`: the chosen lineage's current node gets one child continuing the
  same lineage (same deme) plus one new child lineage per deme in
  `ev.into` -- e.g. SEIR's `infection` event bifurcates an infectious
  lineage into "stays infectious" + "newly exposed". EXCEPTION: if the
  current node is itself `Sample`-typed (this lineage's most recent event
  was a non-destructive sample, and its very next event is this birth), an
  extra plain `Node` is interposed first to hold the bifurcation --
  `NaiveSEIR.singular_part!` (and R phylopomp's own output convention;
  verified empirically, no R-emitted `Sample` node ever has 2 children)
  requires a `Sample` node to have at most 1 child, so the 2 new children
  cannot be attached to it directly.
- `MIGRATION`: relabels the chosen lineage's deme membership; no new node
  (mirrors `swap!` in `src/coloring.jl` -- a migration is not itself an
  observable genealogy feature).
- `DEATH`: the chosen lineage is removed from the inventory; its current
  node is left as an unsampled dead-end tip, to be pruned (see `prune!`).
- `SAMPLE`: a NEW node is created (as `cur`'s child) at time `t` and marked
  `Sample` -- the sample event fires now, at `t`, not at `cur`'s own
  (earlier) creation time, so it must get its own correctly-timestamped
  node rather than relabeling `cur` in place (mirrors R's `sample()`,
  `src/master.h`, which likewise calls `make_node()` to create a fresh node
  at the current time). If `samplemap` is given, the new node's `deme`
  field is set to `samplemap[ev.from]` -- see `simulate`'s docstring.
  Whether the lineage continues afterward is read off `ev.Δ`: if it
  decrements the sampled deme's own compartment count, the sample was
  destructive (lineage removed); otherwise it is non-destructive and the
  new Sample node ITSELF becomes the new open tip in `inv` (no separate
  continuation node -- deliberately: a same-instant placeholder would be a
  zero-length edge below a `Sample` node, and `parse_newick`'s `clip_zlb!`
  refuses to collapse exactly that case, by design, since collapsing it
  would silently overwrite the `Sample` type; a zero-length edge below a
  plain `Node`, as created by the BIRTH exception above, collapses cleanly
  instead -- verified directly, see `check_milestone3.md`).
- `NEUTRAL`: no lineage is involved at all (`ev.from == 0`); only `ev.Δ`
  (already applied by the caller) has any effect.
"""
apply_event!(
    G::Genealogy,
    inv::SimInventory,
    ev::Event,
    model::MGPModel,
    t::Time,
    rng::AbstractRNG,
    samplemap::Union{Nothing,AbstractVector},
) = begin
    ev.type == NEUTRAL && return nothing
    d = ev.from
    b = rand(rng, inv[d])
    cur = G[findfirst(n -> n.name==b, G.nodes)]
    if ev.type == BIRTH
        remove!(inv,d,b)
        parent = cur
        if cur.type == Sample
            mid = push_node!(G,t,Node,b)
            push!(cur.children,mid)
            parent = G.nodes[end]
        end
        c1 = push_node!(G,t,Node,parent.name)
        push!(parent.children,c1)
        add!(inv,d,c1)
        for j ∈ ev.into
            cj = push_node!(G,t,Node,parent.name)
            push!(parent.children,cj)
            add!(inv,j,cj)
        end
    elseif ev.type == MIGRATION
        j = only(ev.into)
        remove!(inv,d,b)
        add!(inv,j,b)
    elseif ev.type == DEATH
        remove!(inv,d,b)
    elseif ev.type == SAMPLE
        remove!(inv,d,b)
        sn = push_node!(G,t,Sample,b)
        push!(cur.children,sn)
        samplenode = G.nodes[end]
        isnothing(samplemap) || (samplenode.deme = samplemap[d])
        sym = model.demes[d]
        destructive = any(p -> p.first==sym && p.second<0, ev.Δ)
        destructive || add!(inv,d,sn)
    else
        error("apply_event!: unhandled event type $(ev.type)") # COV_EXCL_LINE
    end
    nothing
end

# -----------------------------------------------------------------------------
# prune!: remove unsampled dead-end lineages (extinct or extant-uncensored),
# then collapse the resulting degree-1 internal nodes.
# -----------------------------------------------------------------------------

"""
    prune!(G)

Removes lineages that are neither sampled nor ancestral to a sample --
extinct individuals, and lineages still alive but unsampled at the
censoring time -- from a freshly-simulated genealogy `G`, then collapses
any resulting degree-1 internal (`Node`-type) nodes. `Root` nodes left
childless by pruning are left in place for [`repair!`](@ref) (see
`src/genealogy.jl`) to drop, matching its existing "weed out dead roots"
step -- `prune!` never touches `Root`- or `Sample`-typed nodes itself.
"""
prune!(G::Genealogy) = begin
    byname = Dict(n.name => n for n ∈ G.nodes)
    changed = true
    while changed
        changed = false
        for (nm,n) ∈ collect(byname)
            if n.type==Node && isempty(n.children)
                if !isnothing(n.parent) && haskey(byname,n.parent)
                    filter!(!=(nm), byname[n.parent].children)
                end
                delete!(byname,nm)
                changed = true
            end
        end
    end
    changed = true
    while changed
        changed = false
        for (nm,n) ∈ collect(byname)
            if n.type==Node && length(n.children)==1 &&
                !isnothing(n.parent) && haskey(byname,n.parent)
                c = byname[n.children[1]]
                p = byname[n.parent]
                c.parent = p.name
                p.children[findfirst(==(nm),p.children)] = c.name
                delete!(byname,nm)
                changed = true
            end
        end
    end
    empty!(G.nodes)
    append!(G.nodes, values(byname))
    nothing
end

# -----------------------------------------------------------------------------
# simulate: the user-facing entry point.
# -----------------------------------------------------------------------------

"""
    simulate(model::MGPModel, θ; x0, graft, t0 = 0.0, tmax, rng = Random.default_rng(),
             demeset = Unstructured, samplemap = nothing)

Forward-simulate a genealogy from `model` (an [`MGPModel`](@ref), e.g.
[`SEIR`](@ref)/[`MERS`](@ref) from `src/examples/mgp.jl`) using the
Doob-Gillespie direct-method stochastic simulation algorithm.

Arguments:
- `θ`: the model's parameters (any object with the right properties, e.g.
  a `NamedTuple` -- matches what `model.events[i].hazard(x,θ)` expects).
- `x0::NamedTuple`: initial count of every compartment in
  `model.compartments`.
- `graft::AbstractVector{<:Integer}`: number of founding lineages to plant
  at `t0`, one entry per deme in `model.demes` (same order); each founding
  lineage becomes an independent root of the (possibly multi-root, i.e.
  forest-valued) returned genealogy. `x0`'s counts for lineage-tracked
  compartments must equal `graft` exactly: forward simulation tracks every
  extant individual's lineage (unlike a particle filter, which tracks only
  a sparse, partially-observed subset consistent with a *given* tree).
- `tmax`: the censoring time. Any lineage still alive at `tmax` without
  having been sampled is pruned from the returned genealogy (`prune!`).
- `rng`: an `AbstractRNG` for reproducible, independently-seeded draws.
- `demeset::Module`: the demeset (see [`@demes`](@ref)) for the *returned*
  `Genealogy`'s type parameter. Defaults to `Unstructured`, matching
  SEIR's convention that deme membership is latent and unrecorded.
- `samplemap::Union{Nothing,AbstractVector}`: if given, must have one entry
  per deme in `model.demes`; `samplemap[d]` is the `demeset` enum instance
  recorded as `deme` on any node marked `Sample` for MGPModel-deme index
  `d` (e.g. for MERS, `samplemap = [SoftMERS.Camel, SoftMERS.Human]`,
  matching `model.demes == [:I_c, :I_h]`). Root/internal nodes always keep
  `deme = missing` regardless: only the sampling process observes species
  identity directly, matching `NaiveMERS`/`SoftMERS`/etc.'s convention. The
  default `nothing` leaves every node's deme `missing`, matching SEIR.

Piped continuation of a censored simulation (as R phylopomp's
`simulate(x, time=...)` supports) is out of scope for this version.
"""
simulate(
    model::MGPModel,
    θ;
    x0::NamedTuple,
    graft::AbstractVector{<:Integer},
    t0::Real = 0.0,
    tmax::Real,
    rng::AbstractRNG = default_rng(),
    demeset::Module = Unstructured,
    samplemap::Union{Nothing,AbstractVector} = nothing,
) = begin
    ndeme = length(model.demes)
    length(graft) == ndeme ||
        throw(ArgumentError("`graft` must have one entry per deme ($ndeme), got $(length(graft))"))
    isnothing(samplemap) || length(samplemap) == ndeme ||
        throw(ArgumentError("`samplemap` must have one entry per deme ($ndeme), got $(length(samplemap))"))
    for (d,sym) ∈ enumerate(model.demes)
        getproperty(x0,sym) == graft[d] ||
            throw(ArgumentError(
                "x0.$sym = $(getproperty(x0,sym)) must equal graft[$d] = $(graft[d]): "*
                "forward simulation tracks every extant individual's lineage"
            ))
    end

    G = Genealogy{demeset}(Time(t0))
    inv = SimInventory(ndeme)
    for (d,n) ∈ enumerate(graft), _ ∈ Base.OneTo(n)
        r = push_node!(G,Time(t0),Root,nothing)
        c = push_node!(G,Time(t0),Node,r)
        push!(G[r].children,c)
        add!(inv,d,c)
    end

    x = x0
    t::Time = Time(t0)
    haz = Vector{Float64}(undef, length(model.events))
    while true
        for (k,ev) ∈ enumerate(model.events)
            haz[k] = ev.hazard(x,θ)
            @assert haz[k] ≥ 0 "event $(ev.name): negative hazard $(haz[k])"
        end
        total = sum(haz)
        total ≤ 0 && break
        dt = -log(rand(rng))/total
        if t+Time(dt) ≥ Time(tmax)
            break
        end
        t += Time(dt)
        k, _ = rcateg(haz; rng)
        ev = model.events[k]
        x = apply_delta(x,ev)
        apply_event!(G,inv,ev,model,t,rng,samplemap)
        for (i,sym) ∈ enumerate(model.demes)
            @assert length(inv[i])==getproperty(x,sym) "inventory/population "*
                "mismatch in deme $sym: $(length(inv[i])) tracked lineages vs "*
                "$(getproperty(x,sym)) in state -- indicates a bug in the "*
                "engine, not stochastic noise"
        end
    end

    G.time = Time(tmax)
    prune!(G)
    repair!(G)
    G
end
