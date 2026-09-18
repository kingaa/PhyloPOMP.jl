"""
    GuidedMERS

A module containing an implementation of the phylodynamic filter for
the two-host MERS model, using a filter guide and "semisoft"
proposals.  That is, when population-process events occur, the guide
can steer those events preferentially onto (or away from) particular
branches, but the overall event rate remains equal to that of the
underlying population process.
"""
module GuidedMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

"""
    knowledge!(v; deme, type, time)

Return `true` if the guide probabilities are fixed at this node (and,
if so, fill `v` with the appropriate probability vector); return
`false` otherwise.

For MERS the deme is fixed exactly when the genealogy supplies host
metadata, i.e. at the sample tips. This is identical to the default
`known_deme!`.
"""
knowledge!(
    v; deme, type, time,
) = begin
    if ismissing(deme)
        false
    else
        demekron!(v, deme)
        true
    end
end

check(
    gen::Genealogy,
) = begin
    for node ∈ eachindex(gen)
        n = gen[node]
        if n.type == Root
            @assert length(n.children)==1 "wrong number of children ($(length(n.children)) != 1) at root $(n.name), t=$(n.slate)"
        elseif n.type == Sample
            @assert length(n.children)<1 "too many children ($(length(n.children)) > 0) at sample $(n.name), t=$(n.time)"
            @assert n.deme in [Human, Camel] "sample node with unknown deme"
        elseif n.type == Node
            @assert length(n.children)==2 "wrong number of children ($(length(n.children)) ≠ 2) at node $(n.name), t=$(n.time)"
        end
    end
    nothing
end

transmission!(
    alpha, (;S_c, I_c, S_h, I_h), cols;
    β_cc, β_hh, β_ch, β_hc, N_c, N_h,
    _...,
) = begin
    alpha[1] = β_cc*S_c*I_c/N_c
    alpha[2] = β_hh*S_h*I_h/N_h
    alpha[3] = β_hc*S_h*I_c/N_c
    alpha[4] = β_ch*S_c*I_h/N_h
    zero(Prob)
end

removal!(
    alpha, (;S_c, I_c, S_h, I_h), cols;
    γ_c, γ_h,
    _...,
) = begin
    rate_c = γ_c*I_c
    rate_h = γ_h*I_h
    ellC = ell(cols,Camel)
    ellH = ell(cols,Human)
    alpha[1] = @indicator(I_c > ellC, rate_c*(1-ellC/I_c))
    alpha[2] = @indicator(I_h > ellH, rate_h*(1-ellH/I_h))
    rate_c + rate_h - sum(alpha)
end

demography!(
    alpha, (;S_c, I_c, S_h, I_h), cols;
    B_c, B_h, N_c, N_h, _...
        ) = begin
            alpha[1] = B_c
            alpha[2] = B_h
            alpha[3] = B_c*S_c/N_c
            alpha[4] = B_h*S_h/N_h
            zero(Prob)
        end

sampling(
    (;S_c, I_c, S_h, I_h), cols;
    χ_c, χ_h, _...,
) = begin
    χ_c*I_c + χ_h*I_h
end

event_rates!(alpha, state, cols; kwargs...) = begin
    decay = zero(Prob)
    decay += transmission!(@view(alpha[1:4]), state, cols; kwargs...,)
    decay += removal!(@view(alpha[5:6]), state, cols; kwargs...,)
    decay += demography!(@view(alpha[7:10]), state, cols; kwargs...,)
    decay += sampling(state, cols; kwargs...,)
    decay
end

deme_occupancy(;I_c, I_h, _...,) = begin
    ;I_c, I_h
end

live_condition((;S_c, I_c, S_h, I_h), cols) = begin
    ellC, ellH = ell(cols)
    I_c ≥ ellC && I_h ≥ ellH
end

singular_root!(
    gennode, guidenode::GuideNode{F,N,D}, cols, state;
    kwargs...,
) where {F, N, D <: Enum} = begin
    n = deme_occupancy(;state...)
    ells = ell(cols)
    i, _, p = rcateg(guidenode.present[:,1].*(n.-ells), true)
    if i > 0
        live = true
        ll = -log(p)
        plant!(cols,D(i),guidenode.chillins[1])
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, state
end

terminal_sample!(
    gennode, guidenode, cols, (;S_c, I_c, S_h, I_h); χ_c, χ_h, _...,
) = begin
    deme = gennode.deme
    if guidenode.parlin ∈ cols[deme]
        live = true
        chop!(cols, deme, guidenode.parlin)
        if deme == Camel
            ll = log(χ_c*I_c)
            I_c -= 1
        elseif deme == Human
            ll = log(χ_h*I_h)
            I_h -= 1
        else
            @assert false "impossible!"
        end
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, (;S_c, I_c, S_h, I_h)
end

singular_sample!(gennode, guidenode, cols, state; kwargs...) = begin
    terminal_sample!(gennode, guidenode, cols, state; kwargs...)
end

singular_branch!(
    gennode, guidenode, cols, (;S_c, I_c, S_h, I_h);
    β_cc, β_hh, β_hc, β_ch, N_c, N_h,
    _...,
) = begin
    if guidenode.parlin ∈ cols[Camel]
        live = true
        k, _, p = rcateg(
            [
                β_cc*S_c*I_c * guidenode.present[1,1]*guidenode.present[1,2],
                β_hc*S_h*I_c * guidenode.present[1,1]*guidenode.present[2,2],
                β_hc*S_h*I_c * guidenode.present[2,1]*guidenode.present[1,2],
            ],
            true,
        )
        if k==1
            ll = log(β_cc*S_c*I_c/N_c)-log(p)
            fork!(cols,Camel,guidenode.parlin,(Camel,Camel),guidenode.chillins)
            S_c -= 1
            I_c += 1
            ll -= log(I_c*(I_c-1))
        elseif k==2
            ll = log(β_hc*S_h*I_c/N_c)-log(p)
            fork!(cols,Camel,guidenode.parlin,(Camel,Human),guidenode.chillins)
            S_h -= 1
            I_h += 1
            ll -= log(I_c*I_h)
        elseif k==3
            ll = log(β_hc*S_h*I_c/N_c)-log(p)
            fork!(cols,Camel,guidenode.parlin,(Human,Camel),guidenode.chillins)
            S_h -= 1
            I_h += 1
            ll -= log(I_c*I_h)
        else
            live = false
            ll = Prob(-Inf)
        end
    elseif guidenode.parlin ∈ cols[Human]
        live = true
        k, _, p = rcateg(
            [
                β_hh*S_h*I_h * guidenode.present[2,1]*guidenode.present[2,2],
                β_ch*S_c*I_h * guidenode.present[1,1]*guidenode.present[2,2],
                β_ch*S_c*I_h * guidenode.present[2,1]*guidenode.present[1,2],
            ],
            true,
        )
        if k==1
            ll = log(β_hh*S_h*I_h/N_h)-log(p)
            fork!(cols,Human,guidenode.parlin,(Human,Human),guidenode.chillins)
            S_h -= 1
            I_h += 1
            ll -= log(I_h*(I_h-1))
        elseif k==2
            ll = log(β_ch*S_c*I_h/N_h)-log(p)
            fork!(cols,Human,guidenode.parlin,(Camel,Human),guidenode.chillins)
            S_c -= 1
            I_c += 1
            ll -= log(I_c*I_h)
        elseif k==3
            ll = log(β_ch*S_c*I_h/N_h)-log(p)
            fork!(cols,Human,guidenode.parlin,(Human,Camel),guidenode.chillins)
            S_c -= 1
            I_c += 1
            ll -= log(I_c*I_h)
        else
            live = false
            ll = Prob(-Inf)
        end
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, (;S_c, I_c, S_h, I_h)
end

singular_part!(
    cols, state, genealogy, guide, n;
    kwargs...,
) = begin
    ells = ell(cols)
    live = live_condition(state, cols)
    if live
        gennode = genealogy[n]
        guidenode = guide[n]
        if guidenode.type==Root
            singular_root!(gennode, guidenode, cols, state; kwargs...)
        elseif guidenode.type==Sample
            singular_sample!(gennode, guidenode, cols, state; kwargs...)
        elseif guidenode.type==Node
            singular_branch!(gennode, guidenode, cols, state; kwargs...)
        end
    else
        ll = Prob(-Inf)
        ;live, ll, state
    end
end

regular_transmission_cc!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    S_c -= 1
    I_c += 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_transmission_hh!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    S_h -= 1
    I_h += 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_transmission_hc!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    b,p = choose_branch(t,guide,node,I_c,cols,Camel,Human)
    ll = -log(p)
    S_h -= 1
    I_h += 1
    if b == 0
        ellH = ell(cols,Human)
        ll += log(1-ellH/I_h)
    else
        ellC, ellH = swap!(cols,Camel,Human,b)
        ll += log(1-ellC/I_c) - log(I_h)
    end
    ll, (;S_c, I_c, S_h, I_h)
end

regular_transmission_ch!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    b,p = choose_branch(t,guide,node,I_h,cols,Human,Camel)
    ll = -log(p)
    S_c -= 1
    I_c += 1
    if b == 0
        ellC = ell(cols,Camel)
        ll += log(1-ellC/I_c)
    else
        ellC, ellH = swap!(cols,Human,Camel,b)
        ll += log(1-ellH/I_h) - log(I_c)
    end
    ll, (;S_c, I_c, S_h, I_h)
end

regular_removal_c!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    ellC = ell(cols,Camel)
    ll = -log(1-ellC/I_c)       # preboost
    I_c -= 1
    ll, (;S_c, I_c, S_h, I_h)
end

regular_removal_h!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    ellH = ell(cols,Human)
    ll = -log(1-ellH/I_h)       # preboost
    I_h -= 1
    ll, (;S_c, I_c, S_h, I_h)
end

regular_birth_c!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    S_c += 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_birth_h!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    S_h += 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_death_c!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    S_c -= 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_death_h!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h),
    kwargs...,
) = begin
    S_h -= 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_part!(
    cols, state, guide, n, t, tf;
    kwargs...,
) = begin
    alpha = similar(Vector{Prob},10)
    step::Time = zero(Time)
    decay::Prob = zero(Prob)
    ll::Prob = zero(Prob)
    while t < tf
        decay = event_rates!(alpha, state, cols; kwargs...,)
        k, s = rcateg(alpha)
        step = -log(rand())/s
        if k > 0 && t+step < tf
            if k==1
                ll1, state = regular_transmission_cc!(t,guide,n,cols,state)
            elseif k==2
                ll1, state = regular_transmission_hh!(t,guide,n,cols,state)
            elseif k==3
                ll1, state = regular_transmission_hc!(t,guide,n,cols,state)
            elseif k==4
                ll1, state = regular_transmission_ch!(t,guide,n,cols,state)
            elseif k==5
                ll1, state = regular_removal_c!(t,guide,n,cols,state)
            elseif k==6
                ll1, state = regular_removal_h!(t,guide,n,cols,state)
            elseif k==7
                ll1, state = regular_birth_c!(t,guide,n,cols,state)
            elseif k==8
                ll1, state = regular_birth_h!(t,guide,n,cols,state)
            elseif k==9
                ll1, state = regular_death_c!(t,guide,n,cols,state)
            elseif k==10
                ll1, state = regular_death_h!(t,guide,n,cols,state)
            end
            ll += ll1 - decay*step
            t += step
        else
            step = tf - t
            ll -= decay*step
            break
        end
    end
    ll, state
end

mers_rinit(; S_c0, S_h0, I_c0, I_h0, N_c, N_h, _...,) = begin
    m_c = N_c / (S_c0 + I_c0)
    m_h = N_h / (S_h0 + I_h0)
    (
        S_c = round(Int64, m_c*Float64(S_c0)),
        I_c = round(Int64, m_c*Float64(I_c0)),
        S_h = round(Int64, m_h*Float64(S_h0)),
        I_h = round(Int64, m_h*Float64(I_h0)),
    )
end

"""
    filter_pomp(gen, m; β_cc = 4.0, ...)

Constructs a pomp object for the MERS genealogy-conditioned filter, based on
the filter guide built from genealogy `gen` and the guiding finite-state
Markov process `m` (construct `m` with `fsmarkov`).

Pass the built-in `mers_tree` as `gen` to reproduce the default data set; a call
`filter_pomp(mers_tree, fsmarkov(...))` constructs the guided filter.
"""
filter_pomp(
    gen::Genealogy,
    m::FSMarkovProc;
    β_cc = 4.0, β_ch = 1.0, β_hc = 1.0, β_hh = 4.0,
    γ_c = 1.0, γ_h = 1.0,
    χ_c = 1.0, χ_h = 1.0,
    B_c = 10.0, B_h = 10.0,
    S_c0 = 1.0, S_h0 = 1.0,
    I_c0 = 0.01, I_h0 = 0.0,
    N_c = 10000, N_h = 10000,
) = begin
    check(gen)
    guidegen = guide(gen,m,knowledge!)
    pomp(
        params = (
            β_cc = Float64(β_cc), β_ch = Float64(β_ch),
            β_hc = Float64(β_hc), β_hh = Float64(β_hh),
            γ_c = Float64(γ_c), γ_h = Float64(γ_h),
            χ_c = Float64(χ_c), χ_h = Float64(χ_h),
            B_c = Float64(B_c), B_h = Float64(B_h),
            S_c0 = Float64(S_c0), S_h0 = Float64(S_h0),
            I_c0 = Float64(I_c0), I_h0 = Float64(I_h0),
            N_c = Float64(N_c), N_h = Float64(N_h),
        ),
        t0 = timezero(guidegen),
        times = times(guidegen),
        userdata = (genealogy = gen, guide = guidegen,),
        init_state = (
            node=one(Name),
            ll=zero(Prob),
            live=true,
            cols=Coloring(Demes),
            state=(S_c=Int64(0), I_c=Int64(0), S_h=Int64(0), I_h=Int64(0)),
        ),
        rinit = function (; kwargs...)
            (
                node = one(Name),
                ll = zero(Prob),
                live = true,
                cols = Coloring(Demes),
                state=mers_rinit(;kwargs...),
            )
        end,
        rprocess = onestep(
            function (;t, dt, genealogy, guide, node, ll, live, cols, state, kwargs...,)
                if live
                    cols = copy(cols)
                    live, ll, state = singular_part!(
                        cols, state, genealogy, guide, node;
                        kwargs...,
                    )
                else
                    ll = Prob(-Inf)
                end
                if live && dt > 0
                    ll1, state = regular_part!(
                        cols, state, guide, node, t, t+dt;
                        kwargs...,
                    )
                    ll += ll1
                end
                (;node = node+1, ll, live, cols, state)
            end,
        ),
        logdmeasure = function (; live, ll, _...)
            live ? ll : Prob(-Inf)
        end,
    )
end

end
