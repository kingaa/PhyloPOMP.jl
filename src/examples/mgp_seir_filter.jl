# SEIR filter built from the compiler IR. The regular part is compiled.
# The singular part reuses NaiveSEIR.singular_part!, or the generic singular_update! with generic_singular = true.

export compiled_decay, compiled_event_rates!, compiled_regular_part!,
       compiled_filter_pomp

"""
    compiled_decay(model::MGPModel, x, θ, ℓ::AbstractVector{<:Integer},
                    n::AbstractVector{<:Integer}) -> Float64

Return `total_decay` plus, for each DEATH event, the leftover `α_full·1{n>ℓ} − α_reduced`,
where `α_reduced = α_full·(n−ℓ)/n`.
For a hazard of `γ·n` the sum is `γ·ℓ`, the decay of the naive filters.
`ℓ` and `n` are per-deme post-event counts in `model.demes` order.
"""
function compiled_decay(model::MGPModel, x, θ,
                         ℓ::AbstractVector{<:Integer},
                         n::AbstractVector{<:Integer})
    total = Float64(total_decay(model, x, θ, ℓ, n))
    for event in model.events
        event.type == DEATH || continue
        d = event.from
        αfull = Float64(event.hazard(x, θ))
        nd, ℓd = n[d], ℓ[d]
        above = nd > ℓd
        αreduced = above ? αfull * (nd - ℓd) / nd : 0.0
        leftover = (above ? αfull : 0.0) - αreduced
        total += leftover
    end
    total
end

"""
    compiled_event_rates!(alpha, pi, cols, S, E, I, R; β, σ, γ, ω, ψ, χ, pop, model, _...) -> decay

Fill `alpha`/`pi` for the 6 SEIR events and return the decay.
Slots: 1/2 infection no-move/cross, 3/4 progression identity/cross, 5 recovery, 6 waning.
`χ` is the culling rate. Culling is singular, so it enters only the decay.
"""
function compiled_event_rates!(
    alpha, pi_, cols,
    S, E, I, R;
    β, σ, γ, ω, ψ, χ, pop, model::MGPModel,
    _...,
)
    ellE, ellI = ell(cols)
    @assert I ≥ ellI && E ≥ ellE

    infection   = model.events[findfirst(e -> e.name == :infection, model.events)]
    progression = model.events[findfirst(e -> e.name == :progression, model.events)]
    recovery    = model.events[findfirst(e -> e.name == :recovery, model.events)]
    waning      = model.events[findfirst(e -> e.name == :waning, model.events)]

    x = (S = S, E = E, I = I, R = R)
    θ = (β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, N = pop)

    α_inf = Float64(infection.hazard(x, θ))
    alpha[2] = alpha[1] = α_inf
    α_prog = Float64(progression.hazard(x, θ))
    alpha[4] = alpha[3] = α_prog
    alpha[5] = @indicator(I > ellI, γ * (I - ellI))
    alpha[6] = Float64(waning.hazard(x, θ))

    # NaiveProposal: infection splits no-move : move by no_move_share; progression keeps the source-deme split.
    pi_[1] = no_move_share(ellI, ellE, I, E)   # infection, no lineage moves
    pi_[2] = 1 - pi_[1]                         # infection, tracked-I lineage moves to E
    pi_[3] = @indicator(E > 0, 1 - ellE / E)   # progression, untracked-E parent (identity)
    pi_[4] = @indicator(E > 0, ellE / E)       # progression, tracked-E parent (cross)
    pi_[6] = pi_[5] = 1.0

    ℓvec = [ellE, ellI]
    nvec = [E, I]
    compiled_decay(model, x, θ, ℓvec, nvec)
end

"""
    compiled_regular_part!(cols, ll, t, dt, S, E, I, R; model, kwargs...) -> (ll, S, E, I, R)

Compiled counterpart of `NaiveSEIR.regular_part!`.
Same control flow and RNG call order, so a fixed seed gives the same trajectory.
"""
function compiled_regular_part!(
    cols, ll,
    t, dt,
    S, E, I, R;
    model::MGPModel,
    kwargs...,
)
    tf = t + dt
    if t < tf
        alpha = similar(Vector{Prob}, 6)
        pi_ = similar(Vector{Prob}, 6)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellE, ellI = ell(cols)
        while t < tf
            decay = compiled_event_rates!(
                alpha, pi_, cols,
                S, E, I, R;
                model = model, kwargs...,
            )
            k, s = rcateg(alpha .* pi_)
            step = -log(rand()) / s
            if k > 0 && t + step < tf
                ll -= decay * step + log(pi_[k])
                infection = model.events[findfirst(e -> e.name == :infection, model.events)]
                progression = model.events[findfirst(e -> e.name == :progression, model.events)]
                if k == 1
                    S -= 1
                    E += 1
                    ℓpost = [ellE, ellI]
                    npost = [E, I]
                    # No move: Φ_id + ℓ_I·Φ_inl = 1 − ℓ_E/E'. Divided by pi_[1] in the ll update above.
                    ts = full_transitions(infection, ℓpost, npost)
                    Φid = sum((t.phi for t in ts if t isa IdentityTransition); init = 0 // 1)
                    Φinl = sum((t.phi for t in ts if t isa InlineSameDemeTransition); init = 0 // 1)
                    ll += log(Float64(Φid) + ellI * Float64(Φinl))
                elseif k == 2
                    ellI_pre = ellI
                    b = rand(cols[NaiveSEIR.Infec])
                    ellE, ellI = swap!(cols, NaiveSEIR.Infec, NaiveSEIR.Expos, b)
                    S -= 1
                    E += 1
                    ℓpost = [ellE, ellI]
                    npost = [E, I]
                    # CrossDemeTransition is a singleton per event and state.
                    ts = full_transitions(infection, ℓpost, npost)
                    Φcr = only(filter(t -> t isa CrossDemeTransition, ts)).phi
                    ll += log(Float64(Φcr)) - log(1 / ellI_pre)
                elseif k == 3
                    E -= 1
                    I += 1
                    ℓpost = [ellE, ellI]
                    npost = [E, I]
                    # progression has no slot in its own deme, so only IdentityTransition occurs.
                    ts = full_transitions(progression, ℓpost, npost)
                    Φid = only(filter(t -> t isa IdentityTransition, ts)).phi
                    ll += log(Float64(Φid))
                elseif k == 4
                    ellE_pre = ellE
                    b = rand(cols[NaiveSEIR.Expos])
                    ellE, ellI = swap!(cols, NaiveSEIR.Expos, NaiveSEIR.Infec, b)
                    E -= 1
                    I += 1
                    ℓpost = [ellE, ellI]
                    npost = [E, I]
                    ts = full_transitions(progression, ℓpost, npost)
                    Φcr = only(filter(t -> t isa CrossDemeTransition, ts)).phi
                    ll += log(Float64(Φcr)) - log(1 / ellE_pre)
                elseif k == 5
                    ll -= log(1 - ellI / I)
                    I -= 1
                    R += 1
                elseif k == 6
                    R -= 1
                    S += 1
                end
                t += step
            else
                step = tf - t
                ll -= decay * step
                break
            end
        end
        @assert I ≥ ellI && E ≥ ellE
    end
    ll, S, E, I, R
end

"""
    compiled_filter_pomp(gen; β, σ, γ, ω, ψ, χ, pop, S0, E0, I0, R0)

Build the SEIR filter POMP for genealogy `gen`.
Regular part is `compiled_regular_part!` with `model = SEIR`, or the generic `regular_step!` with `generic_regular = true`.
Singular part is `NaiveSEIR.singular_part!`, or the generic `singular_update!` with `generic_singular = true`.
"""
compiled_filter_pomp(
    gen::Genealogy;
    generic_singular = false,
    generic_regular = false,
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
    pop = 100,
    S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    pomp(
        params = (
            β = Float64(β), σ = Float64(σ), γ = Float64(γ),
            ω = Float64(ω), ψ = Float64(ψ), χ = Float64(χ),
            pop = Float64(pop),
            S0 = Float64(S0), E0 = Float64(E0),
            I0 = Float64(I0), R0 = Float64(R0),
        ),
        t0 = timezero(gen),
        times = times(gen),
        rinit = function (; S0, E0, I0, R0, pop, _...)
            m = pop/(S0+E0+I0+R0)
            (
                node = one(Name),
                ll = zero(Prob),
                cols = Coloring(NaiveSEIR.Demes),
                S = round(Int64, m*Float64(S0)),
                E = round(Int64, m*Float64(E0)),
                I = round(Int64, m*Float64(I0)),
                R = round(Int64, m*Float64(R0)),
                live = true,
            )
        end,
        rprocess = onestep(
            function (
                ; node, ll, cols, geneal,
                t, dt,
                S, E, I, R, live,
                args...,
                )
                cols = copy(cols)
                ll = zero(Prob)
                if generic_singular
                    ellE, ellI = ell(cols)
                    if live && (I < ellI || E < ellE)
                        live = false
                    end
                    if live
                        x = (S = S, E = E, I = I, R = R)
                        θ = (β = args[:β], σ = args[:σ], γ = args[:γ], ω = args[:ω],
                             ψ = args[:ψ], χ = args[:χ], N = args[:pop])
                        Δ, x = singular_update!(cols, geneal, node, x, θ, SEIR)
                        ll += Δ
                        S, E, I, R = x.S, x.E, x.I, x.R
                        isfinite(Δ) || (live = false)
                    end
                    live || (ll = Prob(-Inf))
                else
                    ll, S, E, I, R, live = NaiveSEIR.singular_part!(
                        cols, geneal, node, ll, live,
                        S, E, I, R;
                        args...,
                    )
                end
                if live && dt > 0 && isfinite(ll) && generic_regular
                    θ = (β = args[:β], σ = args[:σ], γ = args[:γ], ω = args[:ω],
                         ψ = args[:ψ], χ = args[:χ], N = args[:pop])
                    ll, x, _ = regular_step!(cols, ll, t, dt, (S = S, E = E, I = I, R = R), SEIR, θ)
                    S, E, I, R = x.S, x.E, x.I, x.R
                elseif live && dt > 0 && isfinite(ll)
                    ll, S, E, I, R = compiled_regular_part!(
                        cols, ll, t, dt,
                        S, E, I, R;
                        model = SEIR, args...,
                    )
                end
                (; node = node+1, ll = ll, cols = cols,
                 S = S, E = E, I = I, R = R, live = live)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,),
    )
end
