# MERS filter built from the compiler IR. The regular part is compiled.
# The singular part reuses NaiveMERS.singular_part!, or the generic singular_update! with generic_singular = true.
# TCC/THH: no-move weight is Φ_id + ℓ_d·Φ_inl (C(ℓ,s)-weighted sum).
# THC/TCH cross: boost(Φ_cross, 1/ℓ).

export mers_compiled_event_rates!, mers_compiled_regular_part!, mers_compiled_filter_pomp

"""
    _phi_of(ts::Vector{KLITransition}, ::Type{T}) -> Rational{Int}

Return the `.phi` of the `T` transition in `ts`, or `0//1` if absent.
A transition that needs a tracked lineage in deme d is absent when ℓ_d == 0.
"""
function _phi_of(ts::AbstractVector{<:KLITransition}, ::Type{T}) where {T<:KLITransition}
    i = findfirst(t -> t isa T, ts)
    isnothing(i) ? zero(Rational{Int}) : ts[i].phi
end

"""
    mers_compiled_event_rates!(alpha, pi_, cols, Sc, Ic, Sh, Ih; Beta_cc, Beta_ch,
                           Beta_hc, Beta_hh, gamma_c, gamma_h, chi_c, chi_h,
                           Bc, Bh, Nc, Nh, model, _...) -> decay

Fill `alpha`/`pi_` for the 12 MERS events and return the total decay (`compiled_decay`).
Slots: 1 TCC, 2 THH, 3/4 THC no-move/cross, 5/6 TCH no-move/cross, 7/8 removal, 9/10 birth, 11/12 death.
"""
function mers_compiled_event_rates!(
    alpha, pi_, cols,
    Sc, Ic, Sh, Ih;
    Beta_cc, Beta_ch, Beta_hc, Beta_hh,
    gamma_c, gamma_h, chi_c, chi_h, Bc, Bh, Nc, Nh,
    model::MGPModel,
    _...,
)
    ellc, ellh = ell(cols)
    @assert Ic ≥ ellc && Ih ≥ ellh

    tcc = model.events[findfirst(e -> e.name == :transmission_cc, model.events)]
    thh = model.events[findfirst(e -> e.name == :transmission_hh, model.events)]
    thc = model.events[findfirst(e -> e.name == :transmission_hc, model.events)]
    tch = model.events[findfirst(e -> e.name == :transmission_ch, model.events)]
    birth_c = model.events[findfirst(e -> e.name == :birth_c, model.events)]
    birth_h = model.events[findfirst(e -> e.name == :birth_h, model.events)]
    death_c = model.events[findfirst(e -> e.name == :death_c, model.events)]
    death_h = model.events[findfirst(e -> e.name == :death_h, model.events)]

    x = (S_c = Sc, I_c = Ic, S_h = Sh, I_h = Ih)
    θ = (β_cc = Beta_cc, β_ch = Beta_ch, β_hc = Beta_hc, β_hh = Beta_hh,
         γ_c = gamma_c, γ_h = gamma_h, χ_c = chi_c, χ_h = chi_h,
         B_c = Bc, B_h = Bh, N_c = Nc, N_h = Nh)

    alpha[1] = Float64(tcc.hazard(x, θ))
    alpha[2] = Float64(thh.hazard(x, θ))
    alpha[4] = alpha[3] = Float64(thc.hazard(x, θ))
    alpha[6] = alpha[5] = Float64(tch.hazard(x, θ))
    alpha[7] = @indicator(Ic > ellc, gamma_c * (Ic - ellc))
    alpha[8] = @indicator(Ih > ellh, gamma_h * (Ih - ellh))
    alpha[9]  = Float64(birth_c.hazard(x, θ))
    alpha[10] = Float64(birth_h.hazard(x, θ))
    alpha[11] = Float64(death_c.hazard(x, θ))
    alpha[12] = Float64(death_h.hazard(x, θ))

    pi_[1] = one(Prob)
    pi_[2] = one(Prob)
    pi_[3] = no_move_share(ellc, ellh, Ic, Ih)
    pi_[4] = 1 - pi_[3]
    pi_[5] = no_move_share(ellh, ellc, Ih, Ic)
    pi_[6] = 1 - pi_[5]
    for k in 7:12
        pi_[k] = one(Prob)
    end

    ℓvec = [ellc, ellh]
    nvec = [Ic, Ih]
    compiled_decay(model, x, θ, ℓvec, nvec)
end

"""
    mers_compiled_regular_part!(cols, ll, t, dt, Sc, Ic, Sh, Ih; model, kwargs...)
        -> (ll, Sc, Ic, Sh, Ih)

Compiled counterpart of `NaiveMERS.regular_part!`.
Same control flow and RNG call order, so a fixed seed gives the same trajectory.
"""
function mers_compiled_regular_part!(
    cols, ll,
    t, dt,
    Sc, Ic, Sh, Ih;
    model::MGPModel,
    kwargs...,
)
    tf = t + dt
    if t < tf
        alpha = similar(Vector{Prob}, 12)
        pi_ = similar(Vector{Prob}, 12)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellc, ellh = ell(cols)

        tcc = model.events[findfirst(e -> e.name == :transmission_cc, model.events)]
        thh = model.events[findfirst(e -> e.name == :transmission_hh, model.events)]
        thc = model.events[findfirst(e -> e.name == :transmission_hc, model.events)]
        tch = model.events[findfirst(e -> e.name == :transmission_ch, model.events)]

        while t < tf
            decay = mers_compiled_event_rates!(
                alpha, pi_, cols,
                Sc, Ic, Sh, Ih;
                model = model, kwargs...,
            )
            k, s = rcateg(alpha .* pi_)
            step = -log(rand()) / s
            if t + step < tf
                ll -= decay * step + log(pi_[k])
                if k == 1
                    Sc -= 1
                    Ic += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(tcc, ℓpost, npost)
                    Φid = _phi_of(ts, IdentityTransition)
                    Φinl = _phi_of(ts, InlineSameDemeTransition)
                    # Φ_id + ℓ_c·Φ_inl
                    ll += log(Float64(Φid) + ellc * Float64(Φinl))
                elseif k == 2
                    Sh -= 1
                    Ih += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(thh, ℓpost, npost)
                    Φid = _phi_of(ts, IdentityTransition)
                    Φinl = _phi_of(ts, InlineSameDemeTransition)
                    ll += log(Float64(Φid) + ellh * Float64(Φinl))
                elseif k == 3
                    Sh -= 1
                    Ih += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(thc, ℓpost, npost)
                    Φid = _phi_of(ts, IdentityTransition)
                    Φinl = _phi_of(ts, InlineSameDemeTransition)
                    ll += log(Float64(Φid) + ellc * Float64(Φinl))
                elseif k == 4
                    ellc_pre = ellc
                    b = rand(cols[NaiveMERS.Camel])
                    ellc, ellh = swap!(cols, NaiveMERS.Camel, NaiveMERS.Human, b)
                    Sh -= 1
                    Ih += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(thc, ℓpost, npost)
                    Φcr = _phi_of(ts, CrossDemeTransition)
                    ll += log(Float64(Φcr)) - log(1 / ellc_pre)
                elseif k == 5
                    Sc -= 1
                    Ic += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(tch, ℓpost, npost)
                    Φid = _phi_of(ts, IdentityTransition)
                    Φinl = _phi_of(ts, InlineSameDemeTransition)
                    ll += log(Float64(Φid) + ellh * Float64(Φinl))
                elseif k == 6
                    ellh_pre = ellh
                    b = rand(cols[NaiveMERS.Human])
                    ellc, ellh = swap!(cols, NaiveMERS.Human, NaiveMERS.Camel, b)
                    Sc -= 1
                    Ic += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(tch, ℓpost, npost)
                    Φcr = _phi_of(ts, CrossDemeTransition)
                    ll += log(Float64(Φcr)) - log(1 / ellh_pre)
                elseif k == 7
                    ll -= log(1 - ellc / Ic)
                    Ic -= 1
                elseif k == 8
                    ll -= log(1 - ellh / Ih)
                    Ih -= 1
                elseif k == 9
                    Sc += 1
                elseif k == 10
                    Sh += 1
                elseif k == 11
                    Sc -= 1
                elseif k == 12
                    Sh -= 1
                else
                    @assert false "impossible event" # COV_EXCL_LINE
                end
                t += step
            else
                step = tf - t
                ll -= decay * step
                break
            end
        end
        @assert Ic ≥ ellc && Ih ≥ ellh
    end
    ll, Sc, Ic, Sh, Ih
end

"""
    _mers_theta(args) -> NamedTuple

Parameters of `MERS`, from the keyword parameters of `mers_compiled_filter_pomp`.
"""
_mers_theta(args) = (β_cc = args[:Beta_cc], β_ch = args[:Beta_ch],
                     β_hc = args[:Beta_hc], β_hh = args[:Beta_hh],
                     γ_c = args[:gamma_c], γ_h = args[:gamma_h],
                     χ_c = args[:chi_c], χ_h = args[:chi_h],
                     B_c = args[:Bc], B_h = args[:Bh],
                     N_c = args[:Nc], N_h = args[:Nh])

"""
    mers_compiled_filter_pomp(; Beta_cc, Beta_ch, Beta_hc, Beta_hh, gamma_c,
                          gamma_h, chi_c, chi_h, Bc, Bh, Sc0, Sh0, Ic0, Ih0,
                          Nc, Nh)

Build the MERS filter POMP for genealogy `gen`.
Regular part is `mers_compiled_regular_part!` with `model = MERS`, or the generic `regular_step!` with `generic_regular = true`.
Singular part is `NaiveMERS.singular_part!`, or the generic `singular_update!` with `generic_singular = true`.
"""
mers_compiled_filter_pomp(
    gen::Genealogy;
    generic_singular = false,
    generic_regular = false,
    Beta_cc = 4.0, Beta_ch = 0.0, Beta_hc = 1.0, Beta_hh = 4.0,
    gamma_c = 1.0, gamma_h = 1.0,
    chi_c = 1.0, chi_h = 0.0,
    Bc = 0.1, Bh = 0.03,
    Sc0 = 1.0, Sh0 = 1.0,
    Ic0 = 0.01, Ih0 = 0.0,
    Nc = 10000, Nh = 10000,
) = begin
    pomp(
        params = (
            Beta_cc = Float64(Beta_cc), Beta_ch = Float64(Beta_ch),
            Beta_hc = Float64(Beta_hc), Beta_hh = Float64(Beta_hh),
            gamma_c = Float64(gamma_c), gamma_h = Float64(gamma_h),
            chi_c = Float64(chi_c), chi_h = Float64(chi_h),
            Bc = Float64(Bc), Bh = Float64(Bh),
            Sc0 = Float64(Sc0), Sh0 = Float64(Sh0),
            Ic0 = Float64(Ic0), Ih0 = Float64(Ih0),
            Nc = Float64(Nc), Nh = Float64(Nh),
        ),
        t0 = timezero(gen),
        times = times(gen),
        rinit = function (; Sc0, Sh0, Ic0, Ih0, Nc, Nh, _...)
            fc = Nc / (Sc0 + Ic0)
            fh = Nh / (Sh0 + Ih0)
            (
                node = one(Name),
                ll = zero(Prob),
                cols = Coloring(NaiveMERS.Demes),
                Sc = round(Int64, fc * Sc0),
                Ic = round(Int64, fc * Ic0),
                Sh = round(Int64, fh * Sh0),
                Ih = round(Int64, fh * Ih0),
            )
        end,
        rprocess = onestep(
            function (
                ; node, ll, cols, geneal,
                Sc, Ic, Sh, Ih,
                t, dt,
                args...,
                )
                cols = copy(cols)
                ll = zero(Prob)
                if generic_singular
                    x = (S_c = Sc, I_c = Ic, S_h = Sh, I_h = Ih)
                    θ = (β_cc = args[:Beta_cc], β_ch = args[:Beta_ch],
                         β_hc = args[:Beta_hc], β_hh = args[:Beta_hh],
                         γ_c = args[:gamma_c], γ_h = args[:gamma_h],
                         χ_c = args[:chi_c], χ_h = args[:chi_h],
                         B_c = args[:Bc], B_h = args[:Bh],
                         N_c = args[:Nc], N_h = args[:Nh])
                    Δ, x = singular_update!(cols, geneal, node, x, θ, MERS)
                    ll += Δ
                    Sc, Ic, Sh, Ih = x.S_c, x.I_c, x.S_h, x.I_h
                else
                    ll, Sc, Ic, Sh, Ih = NaiveMERS.singular_part!(
                        cols, ll, geneal, node,
                        Sc, Ic, Sh, Ih;
                        args...,
                    )
                end
                if isfinite(ll) && generic_regular
                    x = (S_c = Sc, I_c = Ic, S_h = Sh, I_h = Ih)
                    ll, x, _ = regular_step!(cols, ll, t, dt, x, MERS, _mers_theta(args))
                    Sc, Ic, Sh, Ih = x.S_c, x.I_c, x.S_h, x.I_h
                elseif isfinite(ll)
                    ll, Sc, Ic, Sh, Ih = mers_compiled_regular_part!(
                        cols, ll, t, dt,
                        Sc, Ic, Sh, Ih;
                        model = MERS, args...,
                    )
                end
                (; node = node + one(Name), ll, cols, Sc, Ic, Sh, Ih)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,),
    )
end
