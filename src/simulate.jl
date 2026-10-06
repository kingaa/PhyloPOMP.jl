# Forward simulation of any `MGPModel` by the Gillespie algorithm.
# Rates come from each `Event`'s `hazard`;
# the jump is read off its `Δ`, `r` and `type`.
#
# `SimInventory` holds the lineages alive in each deme, one entry per
# living individual.
# A `Coloring` is a different thing: it tracks the few
# branches a filter is coloring on a tree it was given.
#
# Deme is stored only when the caller passes `demeset` and `samplemap`;
# then only Sample nodes carry a deme.

import PartiallyObservedMarkovProcesses: simulate
using Random: AbstractRNG, default_rng, randexp

export simulate, simulate_trajectory

"""
    SimInventory(ndeme)

Lineages alive during a simulation, one vector of node names per deme.
`inv[d]` holds the open node of each lineage currently in deme `d`.
"""
struct SimInventory
    live::Vector{Vector{Name}}
    SimInventory(ndeme::Integer) = new([Name[] for _ ∈ Base.OneTo(ndeme)])
end

Base.getindex(inv::SimInventory, d::Integer) = inv.live[d]

## Lineages are addressed by position; order within a deme is irrelevant.
add!(inv::SimInventory, d::Integer, name::Name) = (push!(inv[d], name); nothing)

swapremove!(inv::SimInventory, d::Integer, i::Integer) = begin
    v = inv[d]
    v[i] = v[end]
    pop!(v)
    nothing
end

"""
    apply_delta(x, ev)

Return `x` with `ev.Δ` applied.
Does not touch the genealogy.
"""
apply_delta(x::NamedTuple, ev::Event) = begin
    isempty(ev.Δ) && return x
    delta = Dict(ev.Δ)
    NamedTuple{keys(x)}(map(k -> getfield(x,k) + get(delta,k,0), keys(x)))
end

"""
    push_node!(G, slate, type, parent) -> Name

Append a node and return its name.
The deme is left `missing`.
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

"""
    apply_event!(G, inv, ev, t, rng, samplemap)

Apply one firing of `ev` at time `t`.
The lineage is drawn uniformly from deme `ev.from`.
Nodes are named `1, 2, 3, ...` as they are created, and nothing is
renamed until `prune!` and `repair!`, so the open node of the chosen
lineage is `G.nodes[name]`.

- `BIRTH`: a new `Node` at time `t`, child of the chosen lineage.
  That lineage is replaced by `r[j]` lineages in deme `j`, all open on
  the new node. `validate_model` requires `sum(r) == 2`, which covers
  `fork(I => E, I)`, `fork(I => I, I)`, and `fork(A => B, B)`.
  Each product gets its own node at its next event.
- `MIGRATION`: move the lineage from `ev.from` to `only(ev.into)`.
  No new node.
- `DEATH`: drop the lineage. The open node is left for `prune!`.
- `SAMPLE`: a new `Sample` node at time `t`. With a `samplemap`, its
  deme is `samplemap[ev.from]`. `r[from] == 1` (`sample`) leaves the
  lineage open on that node. `r[from] == 0` (`sample_remove`) ends it.
  No extra node is created at the same instant, so the tree survives
  a Newick round trip.
- `NEUTRAL`: no lineage (`ev.from == 0`). The caller has already
  applied `ev.Δ`.
"""
apply_event!(
    G::Genealogy,
    inv::SimInventory,
    ev::Event,
    t::Time,
    rng::AbstractRNG,
    samplemap::Union{Nothing,AbstractVector},
) = begin
    ev.type == NEUTRAL && return nothing
    d = ev.from
    v = inv[d]
    isempty(v) && throw(AssertionError(
        "event $(ev.name) fired while its source deme $d has no live lineage: "*
        "its hazard must vanish when that compartment is empty"))
    i = rand(rng, eachindex(v))
    b = v[i]
    cur = G.nodes[Int(b)]
    @assert cur.name == b "apply_event!: node names are not sequential"
    if ev.type == BIRTH
        p = push_node!(G,t,Node,b)
        push!(cur.children,p)
        ## r[d] products stay in this deme; the first reuses the drawn slot.
        if ev.r[d] ≥ 1
            v[i] = p
            for _ ∈ 2:ev.r[d]
                add!(inv,d,p)
            end
        else
            swapremove!(inv,d,i)
        end
        for j ∈ eachindex(ev.r)
            j == d && continue
            for _ ∈ Base.OneTo(ev.r[j])
                add!(inv,j,p)
            end
        end
    elseif ev.type == MIGRATION
        swapremove!(inv,d,i)
        add!(inv,only(ev.into),b)
    elseif ev.type == DEATH
        swapremove!(inv,d,i)
    elseif ev.type == SAMPLE
        sn = push_node!(G,t,Sample,b)
        push!(cur.children,sn)
        isnothing(samplemap) || (G.nodes[end].deme = samplemap[d])
        if ev.r[d] == 0
            swapremove!(inv,d,i)
        else
            v[i] = sn
        end
    else
        error("apply_event!: unhandled event type $(ev.type)") # COV_EXCL_LINE
    end
    nothing
end

"""
    prune!(G)

Drop every `Node` with no `Sample` below it, then collapse internal
nodes left with one child.
`Root` and `Sample` nodes are never dropped here; [`repair!`](@ref)
removes childless roots.

`G` must come from the simulator: `G.nodes[k].name == k`, and every
parent precedes its children.
Pass 1 marks Samples and their ancestors.
Pass 2 splices out degree-1 Nodes, parents first.
"""
prune!(G::Genealogy) = begin
    nodes = G.nodes
    n = length(nodes)
    for k ∈ Base.OneTo(n)
        nd = nodes[k]
        nd.name == k && (isnothing(nd.parent) || nd.parent < k) ||
            throw(ArgumentError("prune!: nodes must be named 1:n with parents before children"))
    end
    ## keep[k]: Sample, Root, or ancestor of a Sample.
    keep = falses(n)
    for k ∈ n:-1:1
        nd = nodes[k]
        (nd.type == Sample || nd.type == Root) && (keep[k] = true)
        keep[k] && !isnothing(nd.parent) && (keep[Int(nd.parent)] = true)
    end
    alive = copy(keep)
    for k ∈ Base.OneTo(n)
        keep[k] || continue
        nd = nodes[k]
        filter!(c -> keep[Int(c)], nd.children)
        if nd.type == Node && length(nd.children) == 1 && !isnothing(nd.parent)
            c = nodes[Int(only(nd.children))]
            p = nodes[Int(nd.parent)]
            c.parent = p.name
            p.children[findfirst(==(nd.name), p.children)] = c.name
            alive[k] = false
        end
    end
    keepat!(nodes, alive)
    nothing
end

"""
    simulate(model::MGPModel, θ; x0, graft, t0 = 0.0, tmax, rng = Random.default_rng(),
             demeset = Unstructured, samplemap = nothing)

Simulate a genealogy from `model` by the Doob-Gillespie method.
Returns the pruned `Genealogy` on `[t0, tmax]`.
It is empty when nothing was sampled.
There is no way to resume a run past `tmax`.

Throws `ArgumentError` for a bad `graft`, `samplemap`, or `x0`, for
`tmax < t0`, or when [`validate_model`](@ref) reports an issue with `model`.
A hand-built `MGPModel` is checked the same way as one from `@mgp`.
Throws `AssertionError` when a hazard is negative, infinite, or `NaN`,
or when a hazard is positive while the compartment it needs is empty.

- `θ`: parameters, in the form `event.hazard(x, θ)` expects.
  A `NamedTuple` is fine.
- `x0`: initial count of every compartment.
- `graft`: founding lineages at `t0`, one integer per deme, in
  `model.demes` order. Each founder is a root, so the result can be a
  forest. `x0` for those demes must equal `graft`: every living
  individual is tracked.
- `t0`: start time. Defaults to 0.
- `tmax`: censoring time. A lineage still alive then, and never
  sampled, is removed by `prune!`.
- `rng`: pass your own if the draw should not use the global RNG.
- `demeset`: deme enumeration of the returned `Genealogy`.
  The default, `Unstructured`, stores no deme.
- `samplemap`: one entry per deme, or `nothing`.
  `samplemap[d]` is an instance of the enum in `demeset`, written on
  `Sample` nodes from deme `d`. Other nodes stay `missing`.
  `demeset` has to be that enum's module.
  The default `Unstructured` has no deme enum, so a `samplemap` cannot
  be used with it.
  For MERS, `demeset = SoftMERS.Demes` and
  `samplemap = [SoftMERS.Camel, SoftMERS.Human]`.
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
) = _simulate(model, θ, x0, graft, t0, tmax, rng, demeset, samplemap, nothing)

"""
    simulate_trajectory(model::MGPModel, θ; x0, graft, t0 = 0.0, tmax, rng, demeset, samplemap)

[`simulate`](@ref), plus the population path.
For the same `rng` state, `genealogy` is identical to what
[`simulate`](@ref) returns.
Returns a `NamedTuple`:

- `genealogy`: the pruned genealogy.
- `times`: `[t0, t1, t2, ...]`. `times[end]` is the time of the last
  event, not `tmax`.
- `states`: `x0`, then the state after each event.
- `events`: the name of the event that took `states[i]` to
  `states[i+1]`. One shorter than `times`.

Censoring adds no entry, so `states[end]` is the state at `tmax`.
`states` counts the whole population. Pruning only removes lineages
from the tree.
"""
simulate_trajectory(
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
    times = Float64[]; states = typeof(x0)[]; events = Symbol[]
    record(G, inv, ev, t, x) = begin
        push!(times, Float64(t)); push!(states, x)
        isnothing(ev) || push!(events, ev.name)
    end
    G = _simulate(model, θ, x0, graft, t0, tmax, rng, demeset, samplemap, record)
    (genealogy = G, times = times, states = states, events = events)
end

## `onevent` is `nothing` or a function `(G, inv, ev, t, x)`.
## It is called once before the first event, with `ev = nothing`, `t = t0`
## and `x = x0`, then after every event. It must not draw from `rng`.
## Either way the random draws are the same.
_simulate(
    model::MGPModel,
    θ,
    x0::NamedTuple,
    graft::AbstractVector{<:Integer},
    t0::Real,
    tmax::Real,
    rng::AbstractRNG,
    demeset::Module,
    samplemap::Union{Nothing,AbstractVector},
    onevent,
) = begin
    issues = validate_model(model)
    isempty(issues) ||
        throw(ArgumentError("model $(model.name): " * join(issues, "; ")))
    ndeme = length(model.demes)
    tmax ≥ t0 ||
        throw(ArgumentError("`tmax` = $tmax must not be earlier than `t0` = $t0"))
    length(graft) == ndeme ||
        throw(ArgumentError("`graft` must have one entry per deme ($ndeme), got $(length(graft))"))
    isnothing(samplemap) || length(samplemap) == ndeme ||
        throw(ArgumentError("`samplemap` must have one entry per deme ($ndeme), got $(length(samplemap))"))
    for c ∈ model.compartments
        haskey(x0, c) ||
            throw(ArgumentError("`x0` has no entry for compartment `$c` of model $(model.name)"))
    end
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

    ## Δ as tuples in keys(x0) order, and each deme's position in x0.
    ks = keys(x0)
    δ = [ntuple(j -> sum((p.second for p ∈ ev.Δ if p.first == ks[j]); init = 0), length(ks))
         for ev ∈ model.events]
    demepos = [findfirst(==(sym), ks)::Int for sym ∈ model.demes]

    _simulate_loop!(G, inv, model, θ, x0, δ, demepos, Time(t0), Time(tmax), rng, samplemap, onevent)

    G.time = Time(tmax)
    prune!(G)
    repair!(G)
    check_simulated!(G)
    G
end

## Own function so the compiler sees the length of `δ`'s tuples.
_simulate_loop!(G, inv, model, θ, x0::NamedTuple, δ::Vector{<:Tuple}, demepos, t0::Time,
                tmax::Time, rng, samplemap, onevent) = begin
    x = x0
    t = t0
    isnothing(onevent) || onevent(G, inv, nothing, t, x)
    haz = Vector{Float64}(undef, length(model.events))
    while true
        for (k,ev) ∈ enumerate(model.events)
            haz[k] = ev.hazard(x,θ)
            ## `≥ 0` misses `Inf`. An infinite hazard gives a zero waiting time
            ## and the loop never reaches `tmax`. `NaN` already fails `≥ 0`.
            (isfinite(haz[k]) && haz[k] ≥ 0) ||
                throw(AssertionError("event $(ev.name): hazard $(haz[k]) is not a finite non-negative number"))
        end
        total = sum(haz)
        total ≤ 0 && break
        dt = randexp(rng)/total
        if t+Time(dt) ≥ tmax
            break
        end
        t += Time(dt)
        k, _ = rcateg(haz; rng)
        ev = model.events[k]
        x = typeof(x)(map(+, Tuple(x), δ[k]))
        all(≥(0), Tuple(x)) || throw(AssertionError(
            "event $(ev.name) drove a compartment negative: $x -- its hazard must "*
            "vanish when a compartment it depletes is empty"))
        apply_event!(G,inv,ev,t,rng,samplemap)
        for i ∈ eachindex(demepos)
            length(inv[i]) == x[demepos[i]] || throw(AssertionError(
                "inventory/population mismatch in deme $(model.demes[i]) after event "*
                "$(ev.name): $(length(inv[i])) tracked lineages vs $(x[demepos[i]]) in state"))
        end
        isnothing(onevent) || onevent(G, inv, ev, t, x)
    end
    nothing
end

"""
    check_simulated!(G::Genealogy)

Check the genealogy [`simulate`](@ref) is about to return.
Throws `AssertionError` at the first problem.

After `prune!` and `repair!`:

- names are `1:n` in time order, and parent/child pointers agree;
- every node time lies in `[t0, tmax]`;
- every edge has positive length;
- a `Root` has one child, a `Node` has two, a `Sample` has at most one;
- `nsample` equals the number of `Sample` nodes.
"""
check_simulated!(G::Genealogy) = begin
    n = length(G.nodes)
    ns = 0
    for (k, nd) ∈ enumerate(G.nodes)
        nd.name == k || throw(AssertionError("simulate: node $k is named $(nd.name) after repair!"))
        G.t0 ≤ nd.slate ≤ G.time ||
            throw(AssertionError("simulate: node $k at time $(nd.slate) outside [$(G.t0), $(G.time)]"))
        deg = length(nd.children)
        if nd.type == Root
            isnothing(nd.parent) || throw(AssertionError("simulate: Root $k has a parent"))
            deg == 1 || throw(AssertionError("simulate: Root $k has $deg children (expected 1)"))
        else
            p = nd.parent
            (p isa Integer && 1 ≤ p ≤ n) ||
                throw(AssertionError("simulate: node $k has invalid parent $p"))
            G.nodes[p].slate < nd.slate ||
                throw(AssertionError("simulate: non-positive edge $p -> $k "*
                    "($(G.nodes[p].slate) -> $(nd.slate))"))
            k ∈ G.nodes[p].children ||
                throw(AssertionError("simulate: node $k not among parent $p's children"))
            if nd.type == Node
                deg == 2 || throw(AssertionError("simulate: Node $k has $deg children (expected 2)"))
            elseif nd.type == Sample
                deg ≤ 1 || throw(AssertionError("simulate: Sample $k has $deg children (expected ≤ 1)"))
                ns += 1
            else
                throw(AssertionError("simulate: node $k has unexpected type $(nd.type)"))
            end
        end
        for c ∈ nd.children
            (1 ≤ c ≤ n && G.nodes[c].parent == k) ||
                throw(AssertionError("simulate: child $c of node $k does not point back"))
        end
    end
    G.nsample == ns ||
        throw(AssertionError("simulate: nsample = $(G.nsample) but $ns Sample nodes"))
    nothing
end
