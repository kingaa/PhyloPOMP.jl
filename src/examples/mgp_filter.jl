# Generic KLI filter over an MGPModel (see mgp.jl).
# Regular part: `regular_step!` with the naive proposal (`kli_slots`, `kli_select`, `kli_decay`, `apply_move!`).
# Singular part: `singular_update!`.
# mgp_seir_filter.jl and mgp_mers_filter.jl use them with `generic_regular = true` / `generic_singular = true`.

"""
    kli_hazard(ev, x, θ) -> Float64

Population hazard αᵤ(t,x) of `ev` as Float64.
"""
kli_hazard(ev::Event, x, θ) = Float64(ev.hazard(x, θ))

"""
    apply_pop(x, ev) -> NamedTuple

Return `x` with `ev.Δ` applied.
"""
function apply_pop(x::NamedTuple, ev::Event)
    isempty(ev.Δ) && return x
    names = keys(x)
    delta = Dict(ev.Δ)
    values = map(name -> getfield(x, name) + get(delta, name, 0), names)
    NamedTuple{names}(values)
end

"""
    _deme_counts(x, model) -> Vector{Int}

Host counts of `model.demes` in `x`.
"""
_deme_counts(x::NamedTuple, model::MGPModel) = Int[getfield(x, s) for s in model.demes]

"""
    _phi(ts, pred) -> Rational{Int}

Sum of `phi` over the transitions in `ts` that satisfy `pred` (0 if none).
"""
_phi(ts, pred) = sum((t.phi for t in ts if pred(t)); init = zero(Rational{Int}))

"""
    _no_move_target(ev, ℓ, n) -> Float64

Target factor of a regular `ev` that moves no tracked lineage: `Φ_id + ℓ_a·Φ_inl`, `a = ev.from`,
from `full_transitions` at `ℓ` and the post-event `n`.
`Φ_inl` counts once per tracked lineage in `a` that could be the parent (M09); it is 0 when `r_a = 0`.
Fork and cross transitions are not regular outcomes of this slot.
"""
function _no_move_target(ev::Event, ℓ::AbstractVector{<:Integer}, n::AbstractVector{<:Integer})
    a = ev.from
    ts = full_transitions(ev, ℓ, n)
    Φid = _phi(ts, t -> t isa IdentityTransition)
    Φinl = _phi(ts, t -> t isa InlineSameDemeTransition && t.deme == a)
    Float64(Φid) + ℓ[a] * Float64(Φinl)
end

"""
    _cross_phi(ev, d, ℓ, n) -> Float64

`Φ_cr` of the CrossDemeTransition of `ev` into deme `d`, at `ℓ` after the swap and the post-event `n`
(0 if absent).
"""
function _cross_phi(ev::Event, d::Integer, ℓ::AbstractVector{<:Integer}, n::AbstractVector{<:Integer})
    ts = full_transitions(ev, ℓ, n)
    Float64(_phi(ts, t -> t isa CrossDemeTransition && t.target_deme == d))
end

"""
    _check_demes(cols, model)

Throw `ArgumentError` unless `cols` has one deme per entry of `model.demes`.
The generic code indexes `cols` by position in `model.demes`, so the demeset must list its demes
in the same order as `model.demes`. Only the count can be checked, because the names of the two differ.
"""
function _check_demes(::Coloring{D,N}, model::MGPModel) where {D,N}
    N == length(model.demes) ||
        throw(ArgumentError("coloring has $N demes ($(join(instances(D), ", "))), but model `$(model.name)` " *
                            "has $(length(model.demes)) ($(join(model.demes, ", ")))"))
    nothing
end

"""
    kli_slots(ev) -> Vector{Int}
    kli_slots(model) -> Vector{Tuple{Int,Int}}

Slots of a regular event: `0` is the slot in which no tracked lineage changes deme; `d > 0` is the slot in which
one tracked lineage moves from deme `ev.from` to deme `d`.
A BIRTH or MIGRATION event whose production vector `r` has a slot in a deme `d ≠ ev.from` has slots `[0, d]`;
every other event has `[0]`.
For a model, the `(event index, slot)` pairs of its regular events, in event order.
This is the slot order of `mers_compiled_regular_part!` and `compiled_regular_part!`.

Throws `ArgumentError` for a BIRTH whose parent does not continue (`r[ev.from] == 0`) or whose
products lie in two other demes, and for a regular SAMPLE event.
"""
function kli_slots(ev::Event)
    ev.type == SAMPLE && throw(ArgumentError("kli_slots: regular SAMPLE event `$(ev.name)` is not supported"))
    ev.type in (BIRTH, MIGRATION) || return [0]
    others = [d for d in eachindex(ev.r) if ev.r[d] > 0 && d != ev.from]
    isempty(others) && return [0]
    length(others) == 1 || throw(ArgumentError(
        "kli_slots: event `$(ev.name)` produces into demes $others besides its source; not supported"))
    ev.type == BIRTH && ev.r[ev.from] == 0 && throw(ArgumentError(
        "kli_slots: BIRTH `$(ev.name)` does not continue its parent (r[from] == 0); not supported"))
    [0, only(others)]
end

kli_slots(model::MGPModel) =
    [(i, s) for (i, ev) in enumerate(model.events) if ev.regular for s in kli_slots(ev)]

"""
    kli_select(ev, cols, x, model) -> Vector{Float64}

Naive selection factors π of the slots of regular event `ev` (`kli_slots(ev)` order), at the pre-event state.
The driver of a slot is `α_u·π`.

- One slot (BIRTH, MIGRATION or NEUTRAL): `[1]`.
- DEATH in deme `a`: `[1 - ℓ_a/n_a]` if `n_a > ℓ_a`, else `[0]` (the dying host is untracked).
- MIGRATION from `a`: `[1 - ℓ_a/n_a, ℓ_a/n_a]`, `[0, 0]` if `n_a = 0`. The moving host is the parent, drawn
  uniformly; its lineage moves iff it is tracked, so no-move has probability 0 when `ℓ_a = n_a`.
- BIRTH from `a` into `d`: the target share `[T0/(T0+T1), 1 - T0/(T0+T1)]`, with `T0 = _no_move_target` and
  `T1 = ℓ_a·Φ_cr` (one lineage moved from `a` to `d`), both at the post-event `n`; `[1, 0]` if both are 0.
  This is `no_move_share` (src/guide.jl). `T0 > 0` when `ℓ_a = n_a`: a tracked parent may keep its lineage.
"""
function kli_select(ev::Event, cols::Coloring, x::NamedTuple, model::MGPModel)
    _check_demes(cols, model)
    slots = kli_slots(ev)
    length(slots) == 1 && ev.type != DEATH && return [1.0]
    a = ev.from
    ℓ = collect(Int, ell(cols))
    na = getfield(x, model.demes[a])
    if ev.type == DEATH
        return [na > ℓ[a] ? 1 - ℓ[a] / na : 0.0]
    elseif ev.type == MIGRATION
        return na > 0 ? [1 - ℓ[a] / na, ℓ[a] / na] : [0.0, 0.0]
    end
    d = slots[2]
    n = _deme_counts(apply_pop(x, ev), model)
    t0 = _no_move_target(ev, ℓ, n)
    t1 = 0.0
    if ℓ[a] > 0
        ℓ1 = copy(ℓ)
        ℓ1[a] -= 1
        ℓ1[d] += 1
        t1 = ℓ[a] * _cross_phi(ev, d, ℓ1, n)
    end
    p0 = t0 + t1 > 0 ? t0 / (t0 + t1) : 1.0
    [p0, 1 - p0]
end

"""
    kli_rates!(alpha, pi, slots, cols, x, θ, model) -> decay

Fill `alpha` (event hazard `α_u`) and `pi` (`kli_select`) for each slot in `slots = kli_slots(model)`,
and return `kli_decay`.
"""
function kli_rates!(alpha, pi, slots, cols::Coloring, x::NamedTuple, θ, model::MGPModel)
    k = 1
    while k <= length(slots)
        i = slots[k][1]
        ev = model.events[i]
        α = kli_hazard(ev, x, θ)
        for p in kli_select(ev, cols, x, model)
            slots[k][1] == i || error("kli_rates!: slots do not match kli_slots(model)")
            alpha[k] = α
            pi[k] = p
            k += 1
        end
    end
    kli_decay(alpha, pi, slots, cols, x, θ, model)
end

"""
    kli_decay(alpha, pi, slots, cols, x, θ, model) -> Float64

Decay rate between genealogy events: `total_decay` (λ of M07: SAMPLE hazards and DEATH hazards at `n_a ≤ ℓ_a`)
plus, for each regular event `u`, its target exit rate minus its proposal exit rate `α_u·Σπ` over its slots.
The target exit rate is `α_u·1{n_a > ℓ_a}` for a DEATH in deme `a` and `α_u` otherwise (M08).
`alpha`, `pi` are per slot of `slots = kli_slots(model)`; `ℓ` and `n` are the current counts.

With the naive selection, Σπ = 1 for BIRTH, NEUTRAL and MIGRATION events (MIGRATION when `n_a > 0`), so only DEATH
events add a term, and the result is `compiled_decay` up to rounding; a DEATH adds `γℓ` when its hazard is `γn`.
"""
function kli_decay(alpha, pi, slots, cols::Coloring, x::NamedTuple, θ, model::MGPModel)
    ℓ = collect(Int, ell(cols))
    n = _deme_counts(x, model)
    λ = Float64(total_decay(model, x, θ, ℓ, n))
    k = 1
    while k <= length(slots)
        i = slots[k][1]
        ev = model.events[i]
        α = alpha[k]
        p = 0.0
        while k <= length(slots) && slots[k][1] == i
            p += pi[k]
            k += 1
        end
        a = ev.from
        exit = ev.type == DEATH ? (n[a] > ℓ[a] ? α : 0.0) : α
        λ += exit - α * p
    end
    λ
end

"""
    apply_move!(cols, ev, slot, x1, model; rng = Random.default_rng()) -> Float64

Apply the naive coloring move of `slot` (see `kli_slots`) of regular event `ev`; `x1` is the post-event state.
Returns `log Φ - log q`, `q` the probability of the lineage choice; the regular step charges `-log π` separately.

- DEATH, NEUTRAL: no coloring change, returns 0 (the DEATH weight `1/π` comes from `-log π`).
- Slot 0: no change; `log(Φ_id + ℓ_a·Φ_inl)` (`_no_move_target`) at the current ℓ.
- Slot `d`: draw `b` uniformly from deme `a = ev.from` (`q = 1/ℓ_a`), `swap!` it to `d`, and return
  `log Φ_cr - log(1/ℓ_a)` with `Φ_cr` at ℓ after the swap.

Same RNG call as the cross slots of `mers_compiled_regular_part!` and `compiled_regular_part!`.
"""
function apply_move!(cols::Coloring{D}, ev::Event, slot::Integer, x1::NamedTuple, model::MGPModel;
                     rng::AbstractRNG = default_rng()) where {D}
    _check_demes(cols, model)
    ev.type in (DEATH, NEUTRAL) && return 0.0
    a = ev.from
    n = _deme_counts(x1, model)
    if slot == 0
        return log(_no_move_target(ev, collect(Int, ell(cols)), n))
    end
    ℓa = ell(cols, D(a))
    b = rand(rng, cols[D(a)])
    swap!(cols, D(a), D(slot), b)
    log(_cross_phi(ev, slot, collect(Int, ell(cols)), n)) - log(1 / ℓa)
end

"""
    singular_update!(cols, geneal, node, x, θ, model) -> (Δll, x′)

Singular update at observed genealogy node `node`; `cols` is mutated.
`x` and `θ` are NamedTuples of compartments and parameters, in the names used by `model`.
Returns `-Inf` for `Δll` when the node is incompatible with the current coloring or state.
The state is then still advanced so that the coloring invariants hold.

- Root: plant the lineage in deme `d` with weight `n_d - ℓ_d` (proposal probability `p`, charge `-log p`).
- Sample: the deme `a` is `n.deme`, or the deme holding the lineage if `n.deme` is missing. Candidates are the
  SAMPLE events with `from == a` (with `sum(r) == 1` if the node has a child). With several, draw one with
  weight `α` (probability `q`). Charge `log α + log φ - log q` at the pre-event state, then apply `Δ`.
  `φ` is 1 for a destructive sample (`Δ` decrements deme `a`); `(n_a - ℓ_a)/n_a` for a non-destructive tip, with ℓ
  after the chop; `1/n_a` for a sampled ancestor, whose child lineage goes to the slot deme of `r`.
- Node (two children): draw an event `u` and an ordering of its slot demes with weight
  `α_u/(number of distinct orderings of u)`, over BIRTH events with `sum(r) == 2` and `from` the
  parent's deme. Charge `log α_u + log φ_u - log p`, with `φ_u` the ForkTransition of
  `full_transitions` at the post-fork ℓ and the post-event n.

Same RNG calls as `NaiveMERS.singular_part!`.
The φ for non-destructive samples is the rule used by `NaiveSEIR.singular_part!`; the IR does not derive it.
"""
function singular_update!(cols::Coloring{D}, geneal, node, x, θ, model::MGPModel) where {D}
    _check_demes(cols, model)
    n = geneal[node]
    demes = model.demes
    bump(y, s, k) = merge(y, NamedTuple{(s,)}((getfield(y, s) + k,)))
    deme_of(b) = something(findfirst(d -> Int(b) ∈ cols[d], instances(D)), 0)
    if n.type == Root
        length(n.children) == 1 || error("root $(n.name) has $(length(n.children)) children, expected 1")
        w = Float64[getfield(x, demes[d]) - ell(cols, D(d)) for d in eachindex(demes)]
        i, _, p = rcateg(w, D, true)
        if ismissing(i)
            plant!(cols, D(1), n.lineage)
            return -Inf, bump(x, demes[1], 1)
        end
        plant!(cols, i, n.lineage)
        return -log(p), x
    elseif n.type == Sample
        nkids = length(n.children)
        nkids <= 1 || error("sample $(n.name) has $nkids children, expected at most 1")
        if ismissing(n.deme)
            a = deme_of(n.lineage)
            a == 0 && return -Inf, x
        else
            a = Int(n.deme)
            if n.lineage ∉ cols[n.deme]
                j = deme_of(n.lineage)
                j == 0 && return -Inf, x
                swap!(cols, D(j), n.deme, n.lineage)
                chop!(cols, n.deme, n.lineage)
                evs = [e for e in model.events if e.type == SAMPLE && e.from == a]
                return -Inf, isempty(evs) ? x : apply_pop(bump(x, demes[a], 1), first(evs))
            end
        end
        cands = [e for e in model.events if e.type == SAMPLE && e.from == a && (nkids == 0 || sum(e.r) == 1)]
        isempty(cands) && return -Inf, x
        if length(cands) == 1
            ev, q = only(cands), 1.0
        else
            k, _, q = rcateg(Float64[kli_hazard(e, x, θ) for e in cands], true)
            k == 0 && return -Inf, x
            ev = cands[k]
        end
        destructive = any(pr -> pr.first == demes[a] && pr.second < 0, ev.Δ)
        na = getfield(x, demes[a])
        if nkids == 0
            chop!(cols, D(a), n.lineage)
            phi = destructive ? 1.0 : (na - ell(cols, D(a))) / na
        else
            j = findfirst(==(1), ev.r)
            chop!(cols, D(a), n.lineage, D(j), geneal[only(n.children)].lineage)
            phi = 1.0 / na
        end
        return log(kli_hazard(ev, x, θ)) + log(phi) - log(q), apply_pop(x, ev)
    elseif n.type == Node
        length(n.children) == 2 || error("node $(n.name) has $(length(n.children)) children, expected 2")
        chillins = map(i -> geneal[i].lineage, n.children)
        a = deme_of(n.lineage)
        a == 0 && error("lineage $(n.lineage) of node $(n.name) is in no deme")
        outcomes = Tuple{Event,Vector{Int}}[]
        w = Float64[]
        for e in model.events
            (e.type == BIRTH && sum(e.r) == 2 && e.from == a) || continue
            orders = _orderings(vcat([fill(d, e.r[d]) for d in eachindex(e.r)]...))
            alpha = kli_hazard(e, x, θ)
            for o in orders
                push!(outcomes, (e, o))
                push!(w, alpha / length(orders))
            end
        end
        k, _, p = rcateg(w, true)
        if k == 0
            fork!(cols, D(a), n.lineage, (D(a), D(a)), chillins)
            return -Inf, bump(x, demes[a], 1)
        end
        e, o = outcomes[k]
        fork!(cols, D(a), n.lineage, Tuple(D.(o)), chillins)
        x1 = apply_pop(x, e)
        ts = full_transitions(e, collect(ell(cols)), Int[getfield(x1, s) for s in demes])
        i = findfirst(t -> t isa ForkTransition && sort(t.slot_demes) == sort(o), ts)
        isnothing(i) && return -Inf, x1
        return log(kli_hazard(e, x, θ)) + log(Float64(ts[i].phi)) - log(p), x1
    end
    error("unknown node type $(n.type)")
end

"""
    _orderings(slots) -> Vector{Vector{Int}}

Distinct orderings of the multiset `slots`, in lexicographic order.
"""
function _orderings(slots::AbstractVector{<:Integer})
    length(slots) <= 1 && return [collect(Int, slots)]
    out = Vector{Vector{Int}}()
    for d in sort(unique(slots))
        rest = collect(Int, slots)
        deleteat!(rest, findfirst(==(d), rest))
        for tail in _orderings(rest)
            push!(out, vcat(d, tail))
        end
    end
    out
end

"""
    regular_step!(cols, ll, t, dt, x, model, θ; rng = Random.default_rng()) -> (ll, x, t)

Advance the filter across [t, t+dt) with no observed genealogy events; `cols` is mutated and `t` is returned as `t + dt`.
Slots are `kli_slots(model)`, with rates `α_u·π` (`kli_rates!`) and decay `kli_decay`.
Each iteration draws the slot `k` (`rcateg`), then the waiting time `-log(rand())/Σα·π`, always, even when the
total rate is 0. If the jump falls before `t + dt`, charge `decay·step + log π_k`, apply `Δ`, and add `apply_move!`;
otherwise charge `decay` over the rest of the interval.
Same RNG calls in the same order as `mers_compiled_regular_part!` and `compiled_regular_part!`.
"""
function regular_step!(cols::Coloring, ll::Real, t::Real, dt::Real, x::NamedTuple,
                       model::MGPModel, θ; rng::AbstractRNG = default_rng())
    _check_demes(cols, model)
    tf = t + dt
    slots = kli_slots(model)
    alpha = Vector{Float64}(undef, length(slots))
    pi = Vector{Float64}(undef, length(slots))
    while t < tf
        all(_deme_counts(x, model) .>= ell(cols)) || error("regular_step!: fewer hosts than tracked lineages")
        decay = kli_rates!(alpha, pi, slots, cols, x, θ, model)
        k, s = rcateg(alpha .* pi; rng = rng)
        step = -log(rand(rng)) / s
        if k > 0 && t + step < tf
            ll -= decay * step + log(pi[k])
            i, slot = slots[k]
            ev = model.events[i]
            x = apply_pop(x, ev)
            ll += apply_move!(cols, ev, slot, x, model; rng = rng)
            t += step
        else
            ll -= decay * (tf - t)
            break
        end
    end
    ll, x, tf
end
