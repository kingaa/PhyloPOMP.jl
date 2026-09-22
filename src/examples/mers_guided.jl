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
    alpha, (;Sc, Ic, Sh, Ih), cols;
    βcc, βhh, βch, βhc, Nc, Nh,
    _...,
) = begin
    alpha[1] = βcc*Sc*Ic/Nc
    alpha[2] = βhh*Sh*Ih/Nh
    alpha[3] = βhc*Sh*Ic/Nc
    alpha[4] = βch*Sc*Ih/Nh
    zero(Prob)
end

removal!(
    alpha, (;Sc, Ic, Sh, Ih), cols;
    γc, γh,
    _...,
) = begin
    ratec = γc*Ic
    rateh = γh*Ih
    ellC = ell(cols,Camel)
    ellH = ell(cols,Human)
    alpha[1] = @indicator(Ic > ellC, ratec*(1-ellC/Ic))
    alpha[2] = @indicator(Ih > ellH, rateh*(1-ellH/Ih))
    ratec + rateh - sum(alpha)
end

demography!(
    alpha, (;Sc, Ic, Sh, Ih), cols;
    Bc, Bh, Nc, Nh, _...
        ) = begin
            alpha[1] = Bc
            alpha[2] = Bh
            alpha[3] = Bc*Sc/Nc
            alpha[4] = Bh*Sh/Nh
            zero(Prob)
        end

sampling(
    (;Sc, Ic, Sh, Ih), cols;
    χc, χh, _...,
) = begin
    χc*Ic + χh*Ih
end

event_rates!(alpha, state, cols; kwargs...) = begin
    decay = zero(Prob)
    decay += transmission!(@view(alpha[1:4]), state, cols; kwargs...,)
    decay += removal!(@view(alpha[5:6]), state, cols; kwargs...,)
    decay += demography!(@view(alpha[7:10]), state, cols; kwargs...,)
    decay += sampling(state, cols; kwargs...,)
    decay
end

deme_occupancy(;Ic, Ih, _...,) = begin
    ;Ic, Ih
end

livecondition((;Sc, Ic, Sh, Ih), cols) = begin
    ellC, ellH = ell(cols)
    Ic ≥ ellC && Ih ≥ ellH
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
    gennode, guidenode, cols, (;Sc, Ic, Sh, Ih); χc, χh, _...,
) = begin
    deme = gennode.deme
    if guidenode.parlin ∈ cols[deme]
        live = true
        chop!(cols, deme, guidenode.parlin)
        if deme == Camel
            ll = log(χc*Ic)
            Ic -= 1
        elseif deme == Human
            ll = log(χh*Ih)
            Ih -= 1
        else
            @assert false "impossible!"
        end
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, (;Sc, Ic, Sh, Ih)
end

singular_sample!(gennode, guidenode, cols, state; kwargs...) = begin
    terminal_sample!(gennode, guidenode, cols, state; kwargs...)
end

singular_branch!(
    gennode, guidenode, cols, (;Sc, Ic, Sh, Ih);
    βcc, βhh, βhc, βch, Nc, Nh,
    _...,
) = begin
    if guidenode.parlin ∈ cols[Camel]
        live = true
        k, _, p = rcateg(
            [
                βcc*Sc*Ic * guidenode.present[1,1]*guidenode.present[1,2],
                βhc*Sh*Ic * guidenode.present[1,1]*guidenode.present[2,2],
                βhc*Sh*Ic * guidenode.present[2,1]*guidenode.present[1,2],
            ],
            true,
        )
        if k==1
            ll = log(βcc*Sc*Ic/Nc)-log(p)
            fork!(cols,Camel,guidenode.parlin,(Camel,Camel),guidenode.chillins)
            Sc -= 1
            Ic += 1
            ll -= log(Ic*(Ic-1)/2)
        elseif k==2
            ll = log(βhc*Sh*Ic/Nc)-log(p)
            fork!(cols,Camel,guidenode.parlin,(Camel,Human),guidenode.chillins)
            Sh -= 1
            Ih += 1
            ll -= log(Ic*Ih)
        elseif k==3
            ll = log(βhc*Sh*Ic/Nc)-log(p)
            fork!(cols,Camel,guidenode.parlin,(Human,Camel),guidenode.chillins)
            Sh -= 1
            Ih += 1
            ll -= log(Ic*Ih)
        else
            live = false
            ll = Prob(-Inf)
        end
    elseif guidenode.parlin ∈ cols[Human]
        live = true
        k, _, p = rcateg(
            [
                βhh*Sh*Ih * guidenode.present[2,1]*guidenode.present[2,2],
                βch*Sc*Ih * guidenode.present[1,1]*guidenode.present[2,2],
                βch*Sc*Ih * guidenode.present[2,1]*guidenode.present[1,2],
            ],
            true,
        )
        if k==1
            ll = log(βhh*Sh*Ih/Nh)-log(p)
            fork!(cols,Human,guidenode.parlin,(Human,Human),guidenode.chillins)
            Sh -= 1
            Ih += 1
            ll -= log(Ih*(Ih-1)/2)
        elseif k==2
            ll = log(βch*Sc*Ih/Nh)-log(p)
            fork!(cols,Human,guidenode.parlin,(Camel,Human),guidenode.chillins)
            Sc -= 1
            Ic += 1
            ll -= log(Ic*Ih)
        elseif k==3
            ll = log(βch*Sc*Ih/Nh)-log(p)
            fork!(cols,Human,guidenode.parlin,(Human,Camel),guidenode.chillins)
            Sc -= 1
            Ic += 1
            ll -= log(Ic*Ih)
        else
            live = false
            ll = Prob(-Inf)
        end
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, (;Sc, Ic, Sh, Ih)
end

singular_part!(
    cols, state, genealogy, guide, n;
    kwargs...,
) = begin
    ells = ell(cols)
    live = livecondition(state, cols)
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
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    Sc -= 1
    Ic += 1
    zero(Prob), (;Sc, Ic, Sh, Ih)
end

regular_transmission_hh!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    Sh -= 1
    Ih += 1
    zero(Prob), (;Sc, Ic, Sh, Ih)
end

regular_transmission_hc!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    b,p = choose_branch(t,guide,node,Ic,cols,Camel,Human)
    ll = -log(p)
    Sh -= 1
    Ih += 1
    if b == 0
        ellH = ell(cols,Human)
        ll += log(1-ellH/Ih)
    else
        ellC, ellH = swap!(cols,Camel,Human,b)
        ll += log(1-ellC/Ic) - log(Ih)
    end
    ll, (;Sc, Ic, Sh, Ih)
end

regular_transmission_ch!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    b,p = choose_branch(t,guide,node,Ih,cols,Human,Camel)
    ll = -log(p)
    Sc -= 1
    Ic += 1
    if b == 0
        ellC = ell(cols,Camel)
        ll += log(1-ellC/Ic)
    else
        ellC, ellH = swap!(cols,Human,Camel,b)
        ll += log(1-ellH/Ih) - log(Ic)
    end
    ll, (;Sc, Ic, Sh, Ih)
end

regular_removal_c!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    ellC = ell(cols,Camel)
    ll = -log(1-ellC/Ic)       # preboost
    Ic -= 1
    ll, (;Sc, Ic, Sh, Ih)
end

regular_removal_h!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    ellH = ell(cols,Human)
    ll = -log(1-ellH/Ih)       # preboost
    Ih -= 1
    ll, (;Sc, Ic, Sh, Ih)
end

regular_birth_c!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    Sc += 1
    zero(Prob), (;Sc, Ic, Sh, Ih)
end

regular_birth_h!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    Sh += 1
    zero(Prob), (;Sc, Ic, Sh, Ih)
end

regular_death_c!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    Sc -= 1
    zero(Prob), (;Sc, Ic, Sh, Ih)
end

regular_death_h!(
    t, guide, node, cols, (;Sc, Ic, Sh, Ih),
    kwargs...,
) = begin
    Sh -= 1
    zero(Prob), (;Sc, Ic, Sh, Ih)
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

mers_rinit(; Sc0, Sh0, Ic0, Ih0, Nc, Nh, _...,) = begin
    mc = Nc / (Sc0 + Ic0)
    mh = Nh / (Sh0 + Ih0)
    (
        Sc = round(Int64, mc*Float64(Sc0)),
        Ic = round(Int64, mc*Float64(Ic0)),
        Sh = round(Int64, mh*Float64(Sh0)),
        Ih = round(Int64, mh*Float64(Ih0)),
    )
end

"""
    filter_pomp(gen, m; ...)

Constructs a pomp object for the MERS genealogy-conditioned filter, based on
the filter guide built from genealogy `gen` and the guiding finite-state
Markov process `m` (construct `m` with `fsmarkov`).

Pass the built-in `mers_tree` as `gen` to reproduce the default data set; a call
`filter_pomp(mers_tree, fsmarkov(...))` constructs the guided filter.
"""
filter_pomp(
    gen::Genealogy,
    m::FSMarkovProc;
    βcc = 4.0, βch = 1.0, βhc = 1.0, βhh = 4.0,
    γc = 1.0, γh = 1.0,
    χc = 1.0, χh = 1.0,
    Bc = 10.0, Bh = 10.0,
    Sc0 = 1.0, Sh0 = 1.0,
    Ic0 = 0.01, Ih0 = 0.0,
    Nc = 10000, Nh = 10000,
) = begin
    check(gen)
    guidegen = guide(gen,m,knowledge!)
    pomp(
        params = map(
            Float64,
            (;βcc, βch, βhc, βhh,
             γc, γh, χc, χh, Bc, Bh,
             Sc0, Sh0, Ic0, Ih0, Nc, Nh,
             )
        ),
        t0 = timezero(guidegen),
        times = times(guidegen),
        userdata = (genealogy = gen, guide = guidegen,),
        init_state = (
            node=one(Name),
            ll=zero(Prob),
            live=true,
            cols=Coloring(Demes),
            state=(Sc=Int64(0), Ic=Int64(0), Sh=Int64(0), Ih=Int64(0)),
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
