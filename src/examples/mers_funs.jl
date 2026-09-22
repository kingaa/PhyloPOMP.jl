## This file contains definitions that are used in the guided filters
## for the MERS-CoV two-host (Camel/Human) model.
##
## It is the MERS analogue of seir_funs.jl and is `include`d by both
## SoftMERS (mers_soft.jl) and HardMERS (mers_hard.jl).  Everything
## except the `regular_part!` kernel lives here; each of those two
## modules supplies only its own `regular_part!`, exactly as
## seir_soft.jl / seir_hard.jl do on top of seir_funs.jl.
##
## The singular (genealogical-event) machinery below --- `knowledge!`,
## `check`, `deme_occupancy`, `live_condition`, `singular_root!`,
## `terminal_sample!`, `singular_sample!`, `singular_branch!`,
## `singular_part!`, `mers_rinit` and `filter_pomp` --- is copied
## verbatim from GuidedMERS (mers_guided.jl), so that all three MERS
## kernels share one, identical, singular part and differ *only* in
## how they propose population events on the open intervals.
##
## The two structural differences from seir_funs.jl are:
##
##   1. A guide node (`GuideNode`) carries no `.deme` field, but a MERS
##      Sample's host species (Camel or Human) is *essential* data (it
##      selects χ_c vs χ_h and which color to chop). We therefore pass
##      the underlying `Genealogy` alongside the `Guide` and read the
##      sampled deme from `genealogy[node].deme`. `guide[node]` and
##      `genealogy[node]` share the same node index.
##
##   2. An internal Node (coalescence) is NOT deme-fixed in MERS: the
##      parent can transmit within-host (cc, hh) or across-host (hc, ch).
##      Hence `knowledge!` fixes the guide only at Samples --- exactly the
##      behavior of PhyloPOMP's default `known_deme!`.

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

## ---------------------------------------------------------------------
## Population-process rates.
##
## These are written in seir_funs.jl's `on*`/`off*` keyword style so that
## SoftMERS and HardMERS can share a single `event_rates!`:
##
##   Soft passes  onC = ellC                          (offC defaults to I_c-onC)
##                onH = ellH                          (offH defaults to I_h-onH)
##   Hard passes  onC = sum_relhaz(rh,n,cols,Camel,Human), offC = I_c-ellC
##                onH = sum_relhaz(rh,n,cols,Human,Camel), offH = I_h-ellH
##
## In both cases `pi[k]` is the proposal fraction of the underlying
## population rate that is directed at slot `k`, and `alpha[k]` is the
## resulting proposal intensity `rate*pi[k]`.  Whatever population rate
## is left unclaimed is returned as "decay" and charged continuously
## (`ll -= decay*step`); whichever slot fires is charged `-log(pi[k])`.
##
## Slot layout (12 slots):
##   1  camel → camel transmission
##   2  human → human transmission
##   3  camel → human transmission, identity (no tracked-branch swap)
##   4  camel → human transmission, tracked-branch swap
##   5  human → camel transmission, identity (no tracked-branch swap)
##   6  human → camel transmission, tracked-branch swap
##   7  camel removal
##   8  human removal
##   9  S_c birth
##  10  S_h birth
##  11  S_c death
##  12  S_h death
## ---------------------------------------------------------------------

transmission!(
    alpha, pi, state, cols;
    β_cc, β_hh, β_hc, β_ch, N_c, N_h,
    onC, offC = state.I_c - onC,
    onH, offH = state.I_h - onH,
    _...,
) = begin
    (;S_c, I_c, S_h, I_h) = state
    rate_cc = β_cc*S_c*I_c/N_c
    rate_hh = β_hh*S_h*I_h/N_h
    rate_hc = β_hc*S_h*I_c/N_c   # a camel infects a human
    rate_ch = β_ch*S_c*I_h/N_h   # a human infects a camel
    ## Within-host births cannot move a tracked lineage between demes,
    ## so they are proposed at the full population rate (pi = 1).
    pi[1] = one(Prob)
    pi[2] = one(Prob)
    pi[3] = @indicator(I_c > 0, offC/I_c)
    pi[4] = @indicator(I_c > 0, onC/I_c)
    pi[5] = @indicator(I_h > 0, offH/I_h)
    pi[6] = @indicator(I_h > 0, onH/I_h)
    alpha[1] = rate_cc*pi[1]
    alpha[2] = rate_hh*pi[2]
    alpha[3] = rate_hc*pi[3]
    alpha[4] = rate_hc*pi[4]
    alpha[5] = rate_ch*pi[5]
    alpha[6] = rate_ch*pi[6]
    ## Soft has pi[3]+pi[4] = pi[5]+pi[6] = 1, so this leftover is zero;
    ## Hard's unnormalized onC/onH make it nonzero in either direction.
    rate_hc - alpha[3] - alpha[4] + rate_ch - alpha[5] - alpha[6]
end

removal!(
    alpha, pi, state, cols;
    γ_c, γ_h,
    _...,
) = begin
    (;S_c, I_c, S_h, I_h) = state
    ellC, ellH = ell(cols)
    rate_c = γ_c*I_c
    rate_h = γ_h*I_h
    ## A tracked lineage cannot be removed, so only the (1-ell/I) share
    ## of the removal rate is proposable.
    pi[1] = @indicator(I_c > 0, 1-ellC/I_c)
    pi[2] = @indicator(I_h > 0, 1-ellH/I_h)
    alpha[1] = @indicator(I_c > ellC, rate_c*pi[1])
    alpha[2] = @indicator(I_h > ellH, rate_h*pi[2])
    rate_c + rate_h - sum(alpha)
end

demography!(
    alpha, pi, state, cols;
    B_c, B_h, N_c, N_h,
    _...,
) = begin
    (;S_c, I_c, S_h, I_h) = state
    pi[1] = pi[2] = pi[3] = pi[4] = one(Prob)
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

event_rates!(alpha, pi, state, cols; kwargs...) = begin
    decay = zero(Prob)
    decay += transmission!(
        @view(alpha[1:6]), @view(pi[1:6]), state, cols; kwargs...,
    )
    decay += removal!(
        @view(alpha[7:8]), @view(pi[7:8]), state, cols; kwargs...,
    )
    decay += demography!(
        @view(alpha[9:12]), @view(pi[9:12]), state, cols; kwargs...,
    )
    decay += sampling(state, cols; kwargs...,)
    decay
end

## ---------------------------------------------------------------------
## Singular part: verbatim from GuidedMERS (mers_guided.jl).
## ---------------------------------------------------------------------

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
            ll -= log(I_c*(I_c-1)/2)  # 1/binomial(I_c,2): unordered pair of camel lineages
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
            ll -= log(I_h*(I_h-1)/2)  # 1/binomial(I_h,2): unordered pair of human lineages
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
