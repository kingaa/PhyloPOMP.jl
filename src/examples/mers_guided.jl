"""
    GuidedMERS

A module containing an implementation of the phylodynamic filter for the
MERS-CoV two-host (Camel/Human) model, using a filter guide and "guided"
proposals.  That is, when population-process events occur, the guide can
steer those events preferentially onto (or away from) particular
branches, but the overall event rate remains equal to that of the
underlying population process.

Unlike `SoftMERS`, which preserves the naive identity/tracked-branch mass
split and only guides the choice *within* the tracked-branch group,
`GuidedMERS` jointly normalizes the identity weight (the untracked
source-host count) together with every tracked branch's relative hazard,
then draws from a single categorical over all of them at once.  The
identity/tracked split itself is therefore guide-driven here, not fixed
to the naive ratio --- this is what makes `GuidedMERS` mathematically
distinct from `SoftMERS`.

This is the MERS analogue of `GuidedSEIR` (seir_guided.jl); the joint
chooser uses `choose_branch` from `guide.jl`.  It is self-contained: it
does *not* `include` mers_funs.jl, which serves `SoftMERS` and
`HardMERS`.

Two structural differences from the SEIR case shape the code below:

 1. A `GuideNode` carries no `.deme` field, but a MERS `Sample`'s host
    species (Camel or Human) is *essential* data (it selects χ_c vs χ_h
    and which color to chop).  We therefore carry the underlying
    `Genealogy` alongside the `Guide` in the pomp object's `userdata` and
    read the sampled deme from `geneal[node].deme`.  `guide[node]` and
    `geneal[node]` share the same node index.

 2. An internal `Node` (coalescence) is not deme-fixed in MERS: the
    parent can transmit within-host (cc, hh) or across-host (hc, ch).
    Hence `knowledge!` fixes the guide only at `Sample`s, and a branch
    point must be handled separately according to the parent's deme.
"""
module GuidedMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

const mers_tree = parse_newick(first(mers_trees), t0=0, demes=Demes)

"""
    knowledge!(v; deme, type, time)

This function should return true if the guide probabilities will be fixed
at this node and false otherwise. If the former, it should fill the vector
`v` with an appropriate probability vector.

For MERS the deme is fixed exactly when the genealogy supplies host
metadata, i.e., at the sample tips.  Internal nodes are latent, since a
branch point may be a within-host or a cross-host transmission.
"""
knowledge!(
    v; deme, type, time,
) = begin
    if ismissing(deme)
        false
    else
        demekron!(v,deme)
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
            @assert length(n.children)==0 "MERS sampling is destructive: sample $(n.name), t=$(n.slate), has $(length(n.children)) children"
            @assert !ismissing(n.deme) "MERS samples must carry deme metadata (Camel or Human) at sample $(n.name), t=$(n.slate)"
        elseif n.type == Node
            @assert length(n.children)==2 "wrong number of children ($(length(n.children)) ≠ 2) at node $(n.name), t=$(n.slate)"
        end
    end
    nothing
end

transmission_cc!(alpha, (;S_c, I_c, S_h, I_h), cols; β_cc, N_c, _...,) = begin
    alpha[1] = β_cc*S_c*I_c/N_c
    zero(Prob)
end

transmission_hh!(alpha, (;S_c, I_c, S_h, I_h), cols; β_hh, N_h, _...,) = begin
    alpha[1] = β_hh*S_h*I_h/N_h
    zero(Prob)
end

## THC: camel → human spillover (camel parent, newly infected human)
spillover_hc!(alpha, (;S_c, I_c, S_h, I_h), cols; β_hc, N_c, _...,) = begin
    alpha[1] = β_hc*S_h*I_c/N_c
    zero(Prob)
end

## TCH: human → camel spillover (human parent, newly infected camel)
spillover_ch!(alpha, (;S_c, I_c, S_h, I_h), cols; β_ch, N_h, _...,) = begin
    alpha[1] = β_ch*S_c*I_h/N_h
    zero(Prob)
end

removal_c!(alpha, (;S_c, I_c, S_h, I_h), cols; γ_c, _...,) = begin
    rate = γ_c*I_c
    ellC = ell(cols,Camel)
    alpha[1] = @indicator(I_c > ellC, rate*(1-ellC/I_c))
    rate-sum(alpha)
end

removal_h!(alpha, (;S_c, I_c, S_h, I_h), cols; γ_h, _...,) = begin
    rate = γ_h*I_h
    ellH = ell(cols,Human)
    alpha[1] = @indicator(I_h > ellH, rate*(1-ellH/I_h))
    rate-sum(alpha)
end

demography!(alpha, (;S_c, I_c, S_h, I_h), cols; B_c, B_h, N_c, N_h, _...,) = begin
    alpha[1] = B_c                                  # S_c birth
    alpha[2] = B_h                                  # S_h birth
    alpha[3] = B_c*S_c/N_c                          # S_c death
    alpha[4] = B_h*S_h/N_h                          # S_h death
    zero(Prob)
end

sampling((;S_c, I_c, S_h, I_h), cols; χ_c, χ_h, _...,) = begin
    χ_c*I_c + χ_h*I_h
end

event_rates!(alpha, state, cols; kwargs...) = begin
    decay = zero(Prob)
    decay += transmission_cc!(@view(alpha[1]), state, cols; kwargs...,)
    decay += transmission_hh!(@view(alpha[2]), state, cols; kwargs...,)
    decay += spillover_hc!(@view(alpha[3]), state, cols; kwargs...,)
    decay += spillover_ch!(@view(alpha[4]), state, cols; kwargs...,)
    decay += removal_c!(@view(alpha[5]), state, cols; kwargs...,)
    decay += removal_h!(@view(alpha[6]), state, cols; kwargs...,)
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
    node::GuideNode{F,N,D}, cols, state;
    kwargs...,
) where {F, N, D <: Enum} = begin
    n = deme_occupancy(;state...)
    ells = ell(cols)
    i, _, p = rcateg(node.present[:,1].*(n.-ells), true)
    if i > 0
        live = true
        ll = -log(p)
        plant!(cols,D(i),node.chillins[1])
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, state
end

terminal_sample!(node, deme, cols, (;S_c, I_c, S_h, I_h); χ_c, χ_h, _...) = begin
    χ, I = (deme==Camel) ? (χ_c, I_c) : (χ_h, I_h)
    if node.parlin ∈ cols[deme] && χ*I > 0
        live = true
        ll = log(χ*I)               # α_{S_d} = χ_d·I_d(x'), the pre-removal count
        chop!(cols,deme,node.parlin)
        if deme==Camel
            I_c -= 1
        else
            I_h -= 1
        end
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, (;S_c, I_c, S_h, I_h)
end

singular_sample!(node, deme, cols, state; kwargs...) = begin
    if length(node.chillins) == 0
        terminal_sample!(node, deme, cols, state; kwargs...)
    else                            # MERS sampling is destructive
        live = false
        ll = Prob(-Inf)
        ;live, ll, state
    end
end

## Weights for the three ways a camel parent can explain a branch point:
## a within-camel birth, or either orientation of a camel→human spillover.
camel_branch_options(x::AbstractArray{Prob,2}, rate_cc, rate_hc) = begin
    [
        rate_cc*x[1,1]*x[1,2],      # TCC: (Camel,Camel)
        rate_hc*x[1,1]*x[2,2],      # THC: (Camel,Human)
        rate_hc*x[2,1]*x[1,2],      # THC: (Human,Camel)
    ]
end

## The mirror, for a human parent.
human_branch_options(x::AbstractArray{Prob,2}, rate_hh, rate_ch) = begin
    [
        rate_hh*x[2,1]*x[2,2],      # THH: (Human,Human)
        rate_ch*x[1,1]*x[2,2],      # TCH: (Camel,Human)
        rate_ch*x[2,1]*x[1,2],      # TCH: (Human,Camel)
    ]
end

camel_branch!(
    node, cols, (;S_c, I_c, S_h, I_h);
    β_cc, β_hc, N_c, _...,
) = begin
    rate_cc = @indicator(N_c > 0, β_cc*S_c*I_c/N_c)
    rate_hc = @indicator(N_c > 0, β_hc*S_h*I_c/N_c)
    k, _, p = rcateg(camel_branch_options(node.present,rate_cc,rate_hc), true)
    if k == 0
        live = false
        ll = Prob(-Inf)
    else
        live = true
        ll = -log(p)
        if k == 1                   # within-camel birth
            fork!(cols,Camel,node.parlin,(Camel,Camel),node.chillins)
            S_c -= 1
            I_c += 1
            ll += log(rate_cc)-log(I_c*(I_c-1)/2)
        else                        # camel → human spillover
            if k == 2
                fork!(cols,Camel,node.parlin,(Camel,Human),node.chillins)
            else
                fork!(cols,Camel,node.parlin,(Human,Camel),node.chillins)
            end
            S_h -= 1
            I_h += 1
            ll += log(rate_hc)-log(I_c*I_h)
        end
    end
    ;live, ll, (;S_c, I_c, S_h, I_h)
end

human_branch!(
    node, cols, (;S_c, I_c, S_h, I_h);
    β_hh, β_ch, N_h, _...,
) = begin
    rate_hh = @indicator(N_h > 0, β_hh*S_h*I_h/N_h)
    rate_ch = @indicator(N_h > 0, β_ch*S_c*I_h/N_h)
    k, _, p = rcateg(human_branch_options(node.present,rate_hh,rate_ch), true)
    if k == 0
        live = false
        ll = Prob(-Inf)
    else
        live = true
        ll = -log(p)
        if k == 1                   # within-human birth
            fork!(cols,Human,node.parlin,(Human,Human),node.chillins)
            S_h -= 1
            I_h += 1
            ll += log(rate_hh)-log(I_h*(I_h-1)/2)
        else                        # human → camel spillover
            if k == 2
                fork!(cols,Human,node.parlin,(Camel,Human),node.chillins)
            else
                fork!(cols,Human,node.parlin,(Human,Camel),node.chillins)
            end
            S_c -= 1
            I_c += 1
            ll += log(rate_ch)-log(I_c*I_h)
        end
    end
    ;live, ll, (;S_c, I_c, S_h, I_h)
end

singular_branch!(node, cols, state; kwargs...) = begin
    if node.parlin ∈ cols[Camel]
        camel_branch!(node, cols, state; kwargs...)
    elseif node.parlin ∈ cols[Human]
        human_branch!(node, cols, state; kwargs...)
    else
        live = false
        ll = Prob(-Inf)
        ;live, ll, state
    end
end

singular_part!(
    cols, state, guide, geneal, n;
    kwargs...,
) = begin
    live = live_condition(state, cols)
    if live
        node = guide[n]
        if node.type==Root
            singular_root!(node, cols, state; kwargs...)
        elseif node.type==Sample
            singular_sample!(node, geneal[n].deme, cols, state; kwargs...)
        elseif node.type==Node
            singular_branch!(node, cols, state; kwargs...)
        end
    else
        ll = Prob(-Inf)
        ;live, ll, state
    end
end

regular_transmission_cc!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
    kwargs...,
) = begin
    ellC = ell(cols,Camel)
    S_c -= 1
    I_c += 1
    ll = log(1-ellC*(ellC-1)/I_c/(I_c-1))
    ll, (;S_c, I_c, S_h, I_h)
end

regular_transmission_hh!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
    kwargs...,
) = begin
    ellH = ell(cols,Human)
    S_h -= 1
    I_h += 1
    ll = log(1-ellH*(ellH-1)/I_h/(I_h-1))
    ll, (;S_c, I_c, S_h, I_h)
end

regular_spillover_hc!(          # THC: camel → human
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
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
        ll += log(1-ellC/I_c)-log(I_h)
    end
    ll, (;S_c, I_c, S_h, I_h)
end

regular_spillover_ch!(          # TCH: human → camel
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
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
        ll += log(1-ellH/I_h)-log(I_c)
    end
    ll, (;S_c, I_c, S_h, I_h)
end

regular_removal_c!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
    kwargs...,
) = begin
    ellC = ell(cols,Camel)
    ll = -log(1-ellC/I_c)
    I_c -= 1
    ll, (;S_c, I_c, S_h, I_h)
end

regular_removal_h!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
    kwargs...,
) = begin
    ellH = ell(cols,Human)
    ll = -log(1-ellH/I_h)
    I_h -= 1
    ll, (;S_c, I_c, S_h, I_h)
end

regular_birth_c!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
    kwargs...,
) = begin
    S_c += 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_birth_h!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
    kwargs...,
) = begin
    S_h += 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_death_c!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
    kwargs...,
) = begin
    S_c -= 1
    zero(Prob), (;S_c, I_c, S_h, I_h)
end

regular_death_h!(
    t, guide, node, cols, (;S_c, I_c, S_h, I_h);
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
                ll1, state = regular_spillover_hc!(t,guide,n,cols,state)
            elseif k==4
                ll1, state = regular_spillover_ch!(t,guide,n,cols,state)
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

mers_rinit(; S_c0, I_c0, S_h0, I_h0, N_c, N_h, _...) = begin
    m_c = N_c/(S_c0+I_c0)
    m_h = N_h/(S_h0+I_h0)
    (
        S_c=round(Int64, m_c*Float64(S_c0)),
        I_c=round(Int64, m_c*Float64(I_c0)),
        S_h=round(Int64, m_h*Float64(S_h0)),
        I_h=round(Int64, m_h*Float64(I_h0)),
    )
end

"""
    filter_pomp(g, m; β_cc = 4.0, β_ch = 0.0, β_hc = 0.0, β_hh = 4.0,
         γ_c = 1.0, γ_h = 1.0, χ_c = 1.0, χ_h = 0.0, B_c = 0.0, B_h = 0.0,
         S_c0 = 1.0, S_h0 = 1.0, I_c0 = 0.01, I_h0 = 0.0,
         N_c = 10000, N_h = 10000)

Constructs a pomp object based on the genealogy `g` and finite-state
Markov guiding process `m`.

Pass the built-in `mers_tree` as `g` to reproduce the default data set; a
call `filter_pomp(mers_tree, fsmarkov(...))` constructs the guided filter.
"""
filter_pomp(
    gen::Genealogy,
    m::FSMarkovProc;
    β_cc = 4.0, β_ch = 0.0, β_hc = 0.0, β_hh = 4.0,
    γ_c = 1.0, γ_h = 1.0, χ_c = 1.0, χ_h = 0.0,
    B_c = 0.0, B_h = 0.0,
    S_c0 = 1.0, S_h0 = 1.0, I_c0 = 0.01, I_h0 = 0.0,
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
        userdata = (guide = guidegen, geneal = gen,),
        init_state = (
            node=one(Name),
            ll=zero(Prob),
            live=true,
            cols=Coloring(Demes),
            state=(S_c=Float64(0), I_c=Float64(0), S_h=Float64(0), I_h=Float64(0)),
        ),
        rinit = function (; kwargs...)
            (
                node = one(Name),
                ll = zero(Prob),
                live = true,
                cols = Coloring(Demes),
                state = mers_rinit(;kwargs...),
            )
        end,
        rprocess = onestep(
            function (;t, dt, guide, geneal, node, ll, live, cols, state, kwargs...,)
                if live
                    cols = copy(cols)
                    live, ll, state = singular_part!(
                        cols, state, guide, geneal, node;
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
