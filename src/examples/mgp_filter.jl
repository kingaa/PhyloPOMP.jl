# Generic KLI filter over an MGPModel (see mgp.jl).
# Regular part: `regular_step!` (`kli_slots`, `kli_select`, `kli_decay`, `apply_move!`): the naive proposal without a guide,
# and with one the `:soft`, `:guided` or `:hard` proposal.
# Singular part: `singular_update!`, whose root and fork weights are multiplied by the guide's `present` when there is one.
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
function apply_pop(x::NamedTuple{names}, ev::Event) where {names}
    isempty(ev.Δ) && return x
    values = map(names) do name
        v = getfield(x, name)
        for (s, k) in ev.Δ
            s === name && (v += k)
        end
        v
    end
    NamedTuple{names}(values)
end

"""
    _deme_counts(x, model) -> Vector{Int}

Host counts of `model.demes` in `x`.
"""
_deme_counts(x::NamedTuple, model::MGPModel) = Int[getfield(x, s) for s in model.demes]

## Positions of `model.demes` in the keys of `x`, and the counts at those positions, as `NTuple{N,Int}`.
## The filter state's values are all `Int`, so `Tuple(x)` is an `NTuple` and indexing it is type-stable,
## unlike `getfield(x, s)` with a `Symbol` known only at run time.
_demepos(model::MGPModel, x::NamedTuple, ::Val{N}) where {N} =
    ntuple(k -> findfirst(==(model.demes[k]), keys(x))::Int, Val(N))
_counts(x::NamedTuple, pos::NTuple{N,Int}) where {N} = (v = Tuple(x); ntuple(k -> Int(v[pos[k]]), Val(N)))

"""
    _phi(ts, pred) -> Rational{Int}

Sum of `phi` over the transitions in `ts` that satisfy `pred` (0 if none).
"""
_phi(ts, pred) = sum((t.phi for t in ts if pred(t)); init = zero(Rational{Int}))

"""
    _no_move_target_ir(ev, ℓ, n) -> Float64

Target factor of a regular `ev` that moves no tracked lineage: `Φ_id + ℓ_a·Φ_inl`, `a = ev.from`,
from `full_transitions` at `ℓ` and the post-event `n`.
`Φ_inl` counts once per tracked lineage in `a` that could be the parent (M09); it is 0 when `r_a = 0`.
Fork and cross transitions are not regular outcomes of this slot.
This is the reference; the filter calls `_no_move_target`, which gives the same value without building
the transitions.
"""
function _no_move_target_ir(ev::Event, ℓ::AbstractVector{<:Integer}, n::AbstractVector{<:Integer})
    a = ev.from
    ts = full_transitions(ev, ℓ, n)
    Φid = _phi(ts, t -> t isa IdentityTransition)
    Φinl = _phi(ts, t -> t isa InlineSameDemeTransition && t.deme == a)
    Float64(Φid) + ℓ[a] * Float64(Φinl)
end

"""
    _cross_phi_ir(ev, d, ℓ, n) -> Float64

`Φ_cr` of the CrossDemeTransition of `ev` into deme `d`, at `ℓ` after the swap and the post-event `n`
(0 if absent). Reference for `_cross_phi`.
"""
function _cross_phi_ir(ev::Event, d::Integer, ℓ::AbstractVector{<:Integer}, n::AbstractVector{<:Integer})
    ts = full_transitions(ev, ℓ, n)
    Float64(_phi(ts, t -> t isa CrossDemeTransition && t.target_deme == d))
end

"""
    _phi_one(r, ℓ, n, j) -> Rational{Int}

`kli_binomial_ratio(r, s, ℓ, n)` for `s = 0` (`j = 0`) or `s` with a single 1 in deme `j`, without allocating.
The product of the binomials is reduced once, whereas `kli_binomial_ratio` reduces at every deme. The two are the same
rational while the products fit in `Int`. Past `2^63 - 1`, `kli_binomial_ratio` stays exact or throws `OverflowError`,
while `_phi_one` silently wraps, possibly to a negative value. For the production vectors of regular events
(`e_d`, `2e_a`, `e_a + e_d`) the largest product is `C(n_a,2)` or `n_a·n_d`, so this needs `n` of about 3·10⁹; for
`r = (2,2)` it is about 7.8·10⁴ and for `r = (2,2,2)` about 2·10³.
"""
function _phi_one(r::AbstractVector{<:Integer}, ℓ, n, j::Integer)
    ## one reduction at the end: the same rational (and Float64) as reducing at every step, while the products fit in Int
    num = 1; den = 1
    for d in eachindex(r)
        denom = safe_binomial(n[d], r[d])
        denom == 0 && return zero(Rational{Int})
        num *= safe_binomial(n[d] - ℓ[d], r[d] - (d == j ? 1 : 0))
        den *= denom
    end
    num // den
end

"""
    _phi_one_f(r, ℓ, n, j) -> Float64

`Float64(_phi_one(r, ℓ, n, j))` without forming the fraction. Julia converts a `Rational` to `Float64` by dividing
its numerator by its denominator in `Float64`; the reduced and the unreduced fraction are the same real number, and
when both integers are below 2^53 they are exact in `Float64`, so the division gives the same correctly rounded
result. Above 2^53 it falls back to `_phi_one`.
"""
function _phi_one_f(r::AbstractVector{<:Integer}, ℓ, n, j::Integer)
    num = 1; den = 1
    for d in eachindex(r)
        denom = safe_binomial(n[d], r[d])
        denom == 0 && return 0.0
        num *= safe_binomial(n[d] - ℓ[d], r[d] - (d == j ? 1 : 0))
        den *= denom
    end
    (0 <= num < 2^53 && 0 < den < 2^53) ? Float64(num) / Float64(den) : Float64(num // den)
end

"""
    _no_move_target(ev, ℓ, n) -> Float64

`_no_move_target_ir` without building the transitions: the identity saturation `s = 0` always exists, the
inline one `s = e_a` when `r_a ≥ 1` and `ℓ_a ≥ 1`.
"""
function _no_move_target(ev::Event, ℓ, n)
    a = ev.from
    Φid = _phi_one_f(ev.r, ℓ, n, 0)
    Φinl = (ev.r[a] >= 1 && ℓ[a] >= 1) ? _phi_one_f(ev.r, ℓ, n, a) : 0.0
    Φid + ℓ[a] * Φinl
end

"""
    _cross_phi(ev, d, ℓ, n) -> Float64

`_cross_phi_ir` without building the transitions: the cross saturation `s = e_d` exists when `d ≠ ev.from`,
`r_d ≥ 1` and `ℓ_d ≥ 1`.
"""
_cross_phi(ev::Event, d::Integer, ℓ, n) =
    (d != ev.from && ev.r[d] >= 1 && ℓ[d] >= 1) ? _phi_one_f(ev.r, ℓ, n, d) : 0.0

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
    d = _cross_dest(ev)
    d == 0 ? [0] : [0, d]
end

"""
    _cross_dest(ev) -> Int

The deme `d ≠ ev.from` that a regular BIRTH or MIGRATION produces into, or 0. Throws as `kli_slots` does.
"""
function _cross_dest(ev::Event)
    ev.type == SAMPLE && throw(ArgumentError("kli_slots: regular SAMPLE event `$(ev.name)` is not supported"))
    ev.type in (BIRTH, MIGRATION) || return 0
    d = 0
    for k in eachindex(ev.r)
        (ev.r[k] > 0 && k != ev.from) || continue
        d == 0 || throw(ArgumentError("kli_slots: event `$(ev.name)` produces into demes " *
            "$([j for j in eachindex(ev.r) if ev.r[j] > 0 && j != ev.from]) besides its source; not supported"))
        d = k
    end
    d == 0 && return 0
    ev.type == BIRTH && ev.r[ev.from] == 0 && throw(ArgumentError(
        "kli_slots: BIRTH `$(ev.name)` does not continue its parent (r[from] == 0); not supported"))
    d
end

kli_slots(model::MGPModel) =
    [(i, s) for (i, ev) in enumerate(model.events) if ev.regular for s in kli_slots(ev)]

"""
    kli_select(ev, cols, x, model; lweights = nothing, proposal = :guided) -> Vector{Float64}

Selection factors π of the slots of regular event `ev` (`kli_slots(ev)` order), at the pre-event state.
The driver of a slot is `α_u·π`. Without `lweights` the selection is the naive one:

- One slot (BIRTH, MIGRATION or NEUTRAL): `[1]`.
- DEATH in deme `a`: `[1 - ℓ_a/n_a]` if `n_a > ℓ_a`, else `[0]` (the dying host is untracked).
- MIGRATION from `a`: `[1 - ℓ_a/n_a, ℓ_a/n_a]`, `[0, 0]` if `n_a = 0`. The moving host is the parent, drawn
  uniformly; its lineage moves iff it is tracked, so no-move has probability 0 when `ℓ_a = n_a`.
- BIRTH from `a` into `d`: the target share `[T0/(T0+T1), 1 - T0/(T0+T1)]`, with `T0 = _no_move_target` and
  `T1 = ℓ_a·Φ_cr` (one lineage moved from `a` to `d`), both at the post-event `n`; `[1, 0]` if both are 0.
  This is `no_move_share` (src/guide.jl). `T0 > 0` when `ℓ_a = n_a`: a tracked parent may keep its lineage.

With `lweights(a, d)` (the `relhaz` weights `h_b` of the tracked lineages of deme `a` moving to `d`), `proposal`
changes the BIRTH and MIGRATION shares:

- `:soft`: the naive shares. The lineage is then drawn by `apply_move!` with `q = h_b / Σh`, not `1/ℓ_a`.
- `:guided`: BIRTH `[T0, Σh·Φ_cr]` normalized (the naive `ℓ_a` in `T1` is replaced by `Σh`); MIGRATION
  `[n_a - ℓ_a, Σh]` normalized (the untracked hosts have weight 1 each).
- `:hard`: MIGRATION as `:guided`; BIRTH `[n_a·p0, Σh]` normalized, `p0` the naive no-move share.

All three are unbiased for the same target provided every `h_b > 0`. For `:hard` the shares always sum to 1, unlike
HardSEIR's `off/n`, `on/n`.
"""
function kli_select(ev::Event, cols::Coloring{D,N}, x::NamedTuple, model::MGPModel; lweights = nothing,
                    proposal::Symbol = :guided) where {D,N}
    _check_demes(cols, model)
    p0, p1, ns = _select(ev, ell(cols), x, _demepos(model, x, Val(N)), model, lweights, proposal)
    ns == 1 ? [p0] : [p0, p1]
end

## Shares of the slots of `ev` as `(π_0, π_1, number of slots)`; the arithmetic of `kli_select`, without allocating.
function _select(ev::Event, ℓ::NTuple{N,Int}, x::NamedTuple, pos::NTuple{N,Int}, model::MGPModel, lweights,
                 proposal::Symbol) where {N}
    d = _cross_dest(ev)
    d == 0 && ev.type != DEATH && return (1.0, 0.0, 1)
    a = ev.from
    na = _counts(x, pos)[a]
    if ev.type == DEATH
        return (na > ℓ[a] ? 1 - ℓ[a] / na : 0.0, 0.0, 1)
    elseif ev.type == MIGRATION
        if lweights === nothing || proposal == :soft
            return na > 0 ? (1 - ℓ[a] / na, ℓ[a] / na, 2) : (0.0, 0.0, 2)
        end
        ## guided: an untracked host (weight 1 each) or tracked lineage b (weight h_b) moves
        h = ℓ[a] > 0 ? sum(lweights(a, d)) : 0.0
        s = (na - ℓ[a]) + h
        return s > 0 ? ((na - ℓ[a]) / s, h / s, 2) : (0.0, 0.0, 2)
    end
    n = _counts(apply_pop(x, ev), pos)
    t0 = _no_move_target(ev, ℓ, n)
    t1 = 0.0
    if ℓ[a] > 0
        ℓ1 = ntuple(k -> ℓ[k] - (k == a) + (k == d), Val(N))
        ## naive and soft: each of the ℓ_a lineages has weight 1; guided: weight h_b
        t1 = (lweights === nothing || proposal != :guided ? ℓ[a] : sum(lweights(a, d))) * _cross_phi(ev, d, ℓ1, n)
    end
    p0 = t0 + t1 > 0 ? t0 / (t0 + t1) : 1.0
    if lweights !== nothing && proposal == :hard
        ## hard: no move with weight n_a·p0 (the naive share times the hosts), a move with weight Σ h_b
        off = na * p0
        on = ℓ[a] > 0 ? sum(lweights(a, d)) : 0.0
        p0 = off + on > 0 ? off / (off + on) : 1.0
    end
    (p0, 1 - p0, 2)
end

"""
    kli_rates!(alpha, pi, slots, cols, x, θ, model) -> decay

Fill `alpha` (event hazard `α_u`) and `pi` (`kli_select`) for each slot in `slots = kli_slots(model)`,
and return `kli_decay`.
"""
kli_rates!(alpha, pi, slots, cols::Coloring, x::NamedTuple, θ, model::MGPModel) =
    _kli_rates!(alpha, pi, slots, cols, x, _alphas(_hazards(model), x, θ), model)

"""
    _hazards(model) -> Tuple

The hazard functions of `model.events` as a tuple. `Event.hazard` is typed `Function`, so a call through it
is dispatched at run time; a function that takes this tuple as an argument is compiled for the model's
concrete hazards and calls them directly.
"""
_hazards(model::MGPModel) = Tuple(ev.hazard for ev in model.events)

"""
    _alphas(hz, x, θ) -> NTuple{K,Float64}

`kli_hazard` of every event, in event order.
"""
_alphas(hz::Tuple, x, θ) = map(h -> Float64(h(x, θ)), hz)

function _kli_rates!(alpha, pi, slots, cols::Coloring{D,N}, x::NamedTuple, αs::Tuple, model::MGPModel;
                     lweights = nothing, proposal::Symbol = :guided,
                     pos::NTuple{N,Int} = _demepos(model, x, Val(N))) where {D,N}
    _check_demes(cols, model)
    ℓ = ell(cols)
    k = 1
    while k <= length(slots)
        i = slots[k][1]
        ev = model.events[i]
        α = αs[i]
        p0, p1, ns = _select(ev, ℓ, x, pos, model, lweights, proposal)
        for j in 1:ns
            slots[k][1] == i || error("kli_rates!: slots do not match kli_slots(model)")
            alpha[k] = α
            pi[k] = j == 1 ? p0 : p1
            k += 1
        end
    end
    _kli_decay(alpha, pi, slots, cols, x, αs, model, pos)
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
kli_decay(alpha, pi, slots, cols::Coloring, x::NamedTuple, θ, model::MGPModel) =
    _kli_decay(alpha, pi, slots, cols, x, _alphas(_hazards(model), x, θ), model)

function _kli_decay(alpha, pi, slots, cols::Coloring{D,N}, x::NamedTuple, αs::Tuple, model::MGPModel,
                    pos::NTuple{N,Int} = _demepos(model, x, Val(N))) where {D,N}
    ℓ = ell(cols)
    n = _counts(x, pos)
    ## total_decay (M07), summed in event order as there: SAMPLE hazards, and DEATH hazards when n_a ≤ ℓ_a
    λ = 0.0
    for (i, ev) in enumerate(model.events)
        if ev.type == SAMPLE
            λ += αs[i]
        elseif ev.type == DEATH
            λ += n[ev.from] <= ℓ[ev.from] ? αs[i] : 0.0
        end
    end
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
    apply_move!(cols, ev, slot, x1, model; rng = Random.default_rng(), lweights = nothing) -> Float64

Apply the coloring move of `slot` (see `kli_slots`) of regular event `ev`; `x1` is the post-event state.
Returns `log Φ - log q`, `q` the probability of the lineage choice; the regular step charges `-log π` separately.

- DEATH, NEUTRAL: no coloring change, returns 0 (the DEATH weight `1/π` comes from `-log π`).
- Slot 0: no change; `log(Φ_id + ℓ_a·Φ_inl)` (`_no_move_target`) at the current ℓ.
- Slot `d`: draw `b` from deme `a = ev.from`, `swap!` it to `d`, and return `log Φ_cr - log q` with `Φ_cr` at ℓ
  after the swap. Without `lweights`, `b` is uniform (`q = 1/ℓ_a`). With `lweights`, for every `proposal`, `b` has
  probability `q = h_b / Σh` (`rcateg` over `lweights(a, d)`, in the BitSet's order of `cols[a]`).

Without `lweights`, same RNG call as the cross slots of `mers_compiled_regular_part!` and `compiled_regular_part!`.
"""
function apply_move!(cols::Coloring{D,N}, ev::Event, slot::Integer, x1::NamedTuple, model::MGPModel;
                     rng::AbstractRNG = default_rng(), lweights = nothing,
                     pos::NTuple{N,Int} = _demepos(model, x1, Val(N))) where {D,N}
    _check_demes(cols, model)
    ev.type in (DEATH, NEUTRAL) && return 0.0
    a = ev.from
    n = _counts(x1, pos)
    if slot == 0
        return log(_no_move_target(ev, ell(cols), n))
    end
    if lweights === nothing
        ℓa = ell(cols, D(a))
        b = rand(rng, cols[D(a)])
        q = 1 / ℓa
    else
        ## guided: lineage b of deme a with probability ∝ its weight (`relhaz`), in the BitSet's order
        b, _, q = rcateg(lweights(a, slot), cols[D(a)], true; rng = rng)
    end
    swap!(cols, D(a), D(slot), b)
    log(_cross_phi(ev, slot, ell(cols), n)) - log(q)
end

"""
    singular_update!(cols, geneal, node, x, θ, model; guide = nothing) -> (Δll, x′)

Singular update at observed genealogy node `node`; `cols` is mutated.
`x` and `θ` are NamedTuples of compartments and parameters, in the names used by `model`.
Returns `-Inf` for `Δll` when the node is incompatible with the current coloring or state.
The state is then still advanced so that the coloring invariants hold.

With a `guide`, the root weight is multiplied by `guide[node].present[d, 1]` and the fork weight of an ordering `(o1, o2)`
by `present[o1, 1]·present[o2, 2]`; the sample step does not use it.

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
The sample φ is `kli_binomial_ratio` of the SAMPLE event; for SEIR it equals the rule of `NaiveSEIR.singular_part!`.
"""
function singular_update!(cols::Coloring{D}, geneal, node, x, θ, model::MGPModel;
                          guide::Union{Nothing,Guide} = nothing) where {D}
    _check_demes(cols, model)
    n = geneal[node]
    demes = model.demes
    bump(y, s, k) = merge(y, NamedTuple{(s,)}((getfield(y, s) + k,)))
    deme_of(b) = something(findfirst(d -> Int(b) ∈ cols[d], instances(D)), 0)
    if n.type == Root
        length(n.children) == 1 || error("root $(n.name) has $(length(n.children)) children, expected 1")
        w = Float64[getfield(x, demes[d]) - ell(cols, D(d)) for d in eachindex(demes)]
        guide === nothing || (w .*= guide[node].present[:, 1])
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
        ## φ is the production-slot ratio of the SAMPLE event at ℓ after the move: a tip leaves the slot of a
        ## host that stays (r = e_a) unoccupied, s = 0, giving (n_a − ℓ_a)/n_a; a sampled ancestor's child occupies
        ## it, s = r, giving 1/n_a; a sample that removes the host has r = 0 and φ = 1.
        if nkids == 0
            chop!(cols, D(a), n.lineage)
            s = zeros(Int, length(demes))
        else
            j = findfirst(==(1), ev.r)
            chop!(cols, D(a), n.lineage, D(j), geneal[only(n.children)].lineage)
            s = copy(ev.r)
        end
        phi = Float64(kli_binomial_ratio(ev.r, s, collect(Int, ell(cols)), Int[getfield(x, d) for d in demes]))
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
                gw = guide === nothing ? 1.0 : guide[node].present[o[1], 1] * guide[node].present[o[2], 2]
                push!(w, alpha / length(orders) * gw)
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
    regular_step!(cols, ll, t, dt, x, model, θ; rng = Random.default_rng(), guide = nothing, node = 0,
                  proposal = :guided, maxpop = typemax(Int)) -> (ll, x, t)

Advance the filter across [t, t+dt) with no observed genealogy events; `cols` is mutated and `t` is returned as `t + dt`.
Slots are `kli_slots(model)`, with rates `α_u·π` (`kli_rates!`) and decay `kli_decay`.
Each iteration draws the slot `k` (`rcateg`), then the waiting time `-log(rand())/Σα·π`, always, even when the
total rate is 0. If the jump falls before `t + dt`, charge `decay·step + log π_k`, apply `Δ`, and add `apply_move!`;
otherwise charge `decay` over the rest of the interval.
With `guide` (a `Guide` of the genealogy) and `node` (the guide node whose interval contains `[t, t+dt)`), the move
slots use `proposal` (`:soft`, `:guided` or `:hard`; see `kli_select`) and the moving lineage is drawn by `relhaz`.
Throws `ArgumentError` for any other `proposal`.
Without a guide, same RNG calls in the same order as `mers_compiled_regular_part!` and `compiled_regular_part!`.

`maxpop` caps the hosts of `model.demes` (summed): once an unobserved event takes them above `maxpop`, the step
returns `ll = -Inf` at once. Models without a susceptible pool (LBDP, BDEI, BDSS) otherwise let a fast-growing
particle simulate about exp(r·dt) hidden hosts, which may never finish. The estimate is then the likelihood of the
genealogy with the hosts never above `maxpop` between genealogy events: choose a cap far above any plausible count.
The default, `typemax(Int)`, is no cap.
"""
function regular_step!(cols::Coloring, ll::Real, t::Real, dt::Real, x::NamedTuple,
                       model::MGPModel, θ; rng::AbstractRNG = default_rng(),
                       guide::Union{Nothing,Guide} = nothing, node::Integer = 0, proposal::Symbol = :guided,
                       maxpop::Integer = typemax(Int))
    _check_demes(cols, model)
    proposal in (:soft, :guided, :hard) || throw(ArgumentError("regular_step!: proposal must be :soft, :guided or :hard"))
    _regular_step!(_hazards(model), cols, ll, t, dt, x, model, θ, rng, guide, node, proposal, Int(maxpop))
end

_enough_hosts(n::NTuple{N,Int}, ℓ::NTuple{N,Int}) where {N} = all(ntuple(k -> n[k] >= ℓ[k], Val(N)))

## Function barrier: compiled once per model, with its concrete hazard functions in `hz`.
## With a guide, the lineage that moves is drawn with weights `relhaz` at the time of the previous event
## (as in `choose_move`); for `:guided` and `:hard` the share of the move slot uses the same weights, while `:soft`
## keeps the naive share. Without a guide, all weights are 1.
function _regular_step!(hz::Tuple, cols::Coloring{D,N}, ll::Real, t::Real, dt::Real, x::NamedTuple,
                        model::MGPModel, θ, rng::AbstractRNG, gd, node::Integer, proposal::Symbol = :guided,
                        maxpop::Int = typemax(Int)) where {D,N}
    tf = t + dt
    slots = kli_slots(model)
    alpha = Vector{Float64}(undef, length(slots))
    pi = Vector{Float64}(undef, length(slots))
    w = Vector{Float64}(undef, length(slots))
    pos = _demepos(model, x, Val(N))
    tnow = Ref(Float64(t))
    lw = gd === nothing ? nothing :
        (a, d) -> relhaz(tnow[], gd, node, D(a), D(d), Int[gd[node].linmap[k] for k in cols[D(a)]])
    while t < tf
        _enough_hosts(_counts(x, pos), ell(cols)) || error("regular_step!: fewer hosts than tracked lineages")
        tnow[] = t
        decay = _kli_rates!(alpha, pi, slots, cols, x, _alphas(hz, x, θ), model; lweights = lw, proposal = proposal,
                            pos = pos)
        w .= alpha .* pi
        k, s = rcateg(w; rng = rng)
        step = -log(rand(rng)) / s
        if k > 0 && t + step < tf
            ll -= decay * step + log(pi[k])
            i, slot = slots[k]
            ev = model.events[i]
            x = apply_pop(x, ev)
            sum(_counts(x, pos)) > maxpop && return -Inf, x, tf
            ll += apply_move!(cols, ev, slot, x, model; rng = rng, lweights = lw, pos = pos)
            t += step
        else
            ll -= decay * (tf - t)
            break
        end
    end
    ll, x, tf
end
