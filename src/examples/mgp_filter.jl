# Generic KLI filter over an MGPModel (see mgp.jl).
#WIP
# Outline of a generic filter. Selection, decay, the recoloring weight and the singular update are not written yet and throw;
# the working filters are in mgp_seir_filter.jl and mgp_mers_filter.jl

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
    kli_select(ev, cols, x) -> Float64

Selection factor πᵤ for regular event `ev`; the driver is αᵤ·πᵤ.
Not implemented; throws.
"""
function kli_select(ev::Event, cols, x)
    error("kli_select: not implemented")
end

"""
    kli_decay(alpha, pi, cols, x, model) -> Float64

Decay rate λ(t,x) between genealogy events.
Not implemented; throws.
"""
function kli_decay(alpha, pi, cols, x, model::MGPModel)
    error("kli_decay: not implemented")
end

"""
    apply_move!(cols, ev, x) -> Float64

Apply the coloring move for `ev`; returns `log(phi_u) - log(q)`.
Not implemented; throws.
"""
function apply_move!(cols, ev::Event, x)
    error("apply_move!: not implemented")
end

"""
    singular_update!(cols, node, x, θ, model) -> (Δll, x′)

Singular update at an observed genealogy event:
fork κ, chop χ or swap σ by node type.
Roots enter through the initial condition.
Not implemented; throws.
"""
function singular_update!(cols, node, x, θ, model::MGPModel)
    error("singular_update!: not implemented")
end

"""
    regular_step!(cols, ll, t, dt, x, model, θ) -> (ll, x, t)

Advance the filter across [t, t+dt) with no observed genealogy events.
Proposals differ only in π (`kli_select`) and the coloring proposal (`apply_move!`).
"""
function regular_step!(cols, ll::Float64, t::Float64, dt::Float64, x, model::MGPModel, θ)
    tf = t + dt
    while t < tf
        alpha = Float64[ev.regular ? kli_hazard(ev, x, θ) : 0.0 for ev in model.events]
        pi = Float64[ev.regular ? kli_select(ev, cols, x) : 0.0 for ev in model.events]
        rates = alpha .* pi
        decay = kli_decay(alpha, pi, cols, x, model)
        total = sum(rates)
        if total <= 0
            ll -= decay * (tf - t)
            break
        end
        step = -log(rand()) / total
        if t + step < tf
            k, _ = rcateg(rates)
            ll -= decay*step + log(pi[k])
            x = apply_pop(x, model.events[k])
            ll += apply_move!(cols, model.events[k], x)
            t += step
        else
            ll -= decay * (tf - t)
            break
        end
    end
    ll, x, t
end
