# Filter pomp object for any MGPModel, built from the generic singular and regular steps.

export mgp_filter_pomp, generic_demeset

@demes MGPDemes2 d1 d2
@demes MGPDemes3 d1 d2 d3

"""
    generic_demeset(model) -> Module

A demeset with one deme per entry of `model.demes`, in the same order: `Unstructured` for one deme,
`MGPDemes2` or `MGPDemes3` for two or three. Use it for genealogies whose sample demes are not
observed. A genealogy that records sample demes needs the demeset it was parsed with.
"""
function generic_demeset(model::MGPModel)
    n = length(model.demes)
    n == 1 && return Unstructured
    n == 2 && return MGPDemes2
    n == 3 && return MGPDemes3
    throw(ArgumentError("generic_demeset: model `$(model.name)` has $n demes; pass a demeset"))
end

"""
    mgp_filter_pomp(gen, model; θ, x0, demeset = generic_demeset(model), guide = nothing, proposal = :guided,
                    maxpop = nothing)

Particle-filter `pomp` object for genealogy `gen` under `model`. Without a guide the proposal is the naive one.
Each step applies `singular_update!` at the next genealogy node, then `regular_step!` up to the next node time.

- `θ`: the model's parameters, a `NamedTuple` with the names the hazards use (`N` for SIR's population size).
  They become the pomp parameters.
- `x0`: the initial count of every compartment of `model`, as non-negative integers.
- `demeset`: the coloring's demeset; its demes must be listed in the order of `model.demes`.
- `guide`: a `FSMarkovProc` over `demeset`'s demes (built with `fsmarkov`), or `nothing` for the naive
  proposal. With a guide, the root deme, the demes of a fork's children and the lineage that moves at a
  regular event are drawn with weights from `Guide(gen, guide, _model_knowledge(model, D))` (`present`,
  `relhaz`), as in the hand-written guided filters. The target, and so the expected likelihood, is the same.
  The guide fixes the deme of a node that the genealogy records; for one it does not record, it uses the table:
  samples are in the deme of the SAMPLE events and branch points in the parent deme of the two-product BIRTH
  events, when that deme is the same for all of them (SEIR and BDEI: I). Without that, a genealogy with no demes
  (`Unstructured`) gives every lineage the same weight and the guide does nothing.
- `proposal` (with a guide): how the guide sets the share of the slot in which a lineage moves, as in the
  hand-written kernels of the same names. `:soft` keeps the naive share; `:guided` (default) weights the move
  by `Φ_cr·Σ relhaz` for a birth and `Σ relhaz` against `n_a − ℓ_a` for a migration; `:hard` does the same for a
  migration and weighs `Σ relhaz` against `n_a` times the naive no-move share for a birth. All three draw the
  moving lineage by `relhaz` and use the guide at nodes. Ignored without a guide.
- `maxpop` (default `nothing`, no cap): the most hosts of `model.demes`, summed, a particle may hold between
  genealogy events (see `regular_step!`). Needed for models without a susceptible pool (LBDP, BDEI, BDSS), where a
  particle with fast growth may otherwise never finish a step. The estimate is then the likelihood with the hosts
  never above the cap between genealogy events, so choose it far above any plausible count.

A particle whose step returns `-Inf` stays at `-Inf` for the rest of the genealogy.
Throws `ArgumentError` when `x0` does not name every compartment exactly once or holds a negative or
non-integer count.
"""
function mgp_filter_pomp(gen::Genealogy, model::MGPModel; θ::NamedTuple, x0::NamedTuple,
                         demeset::Module = generic_demeset(model),
                         guide::Union{Nothing,FSMarkovProc} = nothing, proposal::Symbol = :guided,
                         maxpop::Union{Nothing,Integer} = nothing)
    proposal in (:soft, :guided, :hard) ||
        throw(ArgumentError("mgp_filter_pomp: proposal must be :soft, :guided or :hard"))
    comps = Tuple(model.compartments)
    Set(keys(x0)) == Set(comps) && length(x0) == length(comps) ||
        throw(ArgumentError("mgp_filter_pomp: x0 must name each of $(join(comps, ", ")) once; got $(keys(x0))"))
    all(v -> isinteger(v) && v >= 0, values(x0)) ||
        throw(ArgumentError("mgp_filter_pomp: x0 counts must be non-negative integers; got $x0"))
    xinit = NamedTuple{comps}(Tuple(Int(getfield(x0, c)) for c in comps))
    pnames = keys(θ)
    vp = Val(pnames)
    cap = maxpop === nothing ? typemax(Int) : Int(maxpop)
    cap >= sum(Int(getfield(xinit, d)) for d in model.demes) ||
        throw(ArgumentError("mgp_filter_pomp: maxpop = $cap is below the initial number of hosts"))
    pos = _demepos(model, xinit, Val(length(model.demes)))   # deme positions in the state, computed once
    gd = guide === nothing ? nothing : Guide(gen, guide, _model_knowledge(model, demeset.DemeSet))
    pomp(
        params = map(Float64, θ),
        t0 = timezero(gen),
        times = times(gen),
        rinit = function (; _...)
            (node = one(Name), ll = zero(Prob), cols = Coloring(demeset), x = xinit, live = true)
        end,
        rprocess = onestep(
            function (; node, ll, cols, x, live, geneal, t, dt, args...)
                cols = copy(cols)
                ll = zero(Prob)
                th = _float_params(vp, values(args))
                if live
                    Δ, x = singular_update!(cols, geneal, node, x, th, model; guide = gd)
                    ll += Δ
                    live = isfinite(Δ) && _enough_hosts(_counts(x, pos), ell(cols))
                end
                if live && dt > 0
                    ll, x, _ = regular_step!(cols, ll, t, dt, x, model, th; guide = gd, node = node, proposal = proposal,
                                             maxpop = cap)
                    isfinite(ll) || (live = false)
                end
                live || (ll = Prob(-Inf))
                (; node = node + one(Name), ll, cols, x, live)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,),
    )
end

## The model parameters from the keyword arguments of `rprocess`, as Float64. With the names as a type parameter
## the selection is type-stable; a generator over a tuple of names held in a variable is not.
_float_params(::Val{names}, nt::NamedTuple) where {names} = NamedTuple{names}(map(Float64, Tuple(NamedTuple{names}(nt))))

## The guide's `knowledge!` from the model table: a deme recorded in the genealogy is used as it is; otherwise a
## sample is in the deme of the SAMPLE events and a branch point in the parent deme of the two-product BIRTH events,
## when that deme is the same for all of them. (Aaron's GuidedSEIR states the SEIR case by hand: Sample and Node in I.)
function _model_knowledge(model::MGPModel, ::Type{D}) where {D}
    sd = unique(ev.from for ev in model.events if ev.type == SAMPLE)
    fd = unique(ev.from for ev in model.events if ev.type == BIRTH && sum(ev.r) == 2)
    sdeme = length(sd) == 1 ? instances(D)[only(sd)] : nothing
    fdeme = length(fd) == 1 ? instances(D)[only(fd)] : nothing
    function (v; deme, type, time)
        if !ismissing(deme)
            demekron!(v, deme)
            true
        elseif type == Sample && sdeme !== nothing
            demekron!(v, sdeme)
            true
        elseif type == Node && fdeme !== nothing
            demekron!(v, fdeme)
            true
        else
            false
        end
    end
end
