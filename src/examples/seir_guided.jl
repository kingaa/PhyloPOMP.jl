"""
    GuidedSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using a filter guide and "soft" proposals.  That is, when
population-process events occur, the guide can steer those events
preferentially onto (or away from) particular branches, but the overall
event rate remains equal to that of the underlying population process.
"""
module GuidedSEIR

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Expos Infec
using .Demes: Expos, Infec, DemeSet

include("seir_trees.jl")

"""
    knowledge!(v; deme, type, time)

This function should return true if the guide probabilities will be fixed
at this node and false otherwise. If the former, it should fill the vector
`v` with an appropriate probability vector.
"""
knowledge!(
    v; deme, type, time,
) = begin
    if type==Sample || type==Node
        demekron!(v,Infec)
        true
    else
        false
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
            @assert length(n.children)<2 "too many children ($(length(n.children)) > 1) at sample $(n.name), t=$(n.time)"
        elseif n.type == Node
            @assert length(n.children)==2 "wrong number of children ($(length(n.children)) ≠ 2) at node $(n.name), t=$(n.time)"
        end
    end
    nothing
end

transmission!(alpha, (;S, E, I, R), cols; β, pop, _...,) = begin
    alpha[1] = β*S*I/pop
    zero(Prob)
end

progression!(alpha, (;S, E, I, R), cols; σ, _...,) = begin
    alpha[1] = σ*E
    zero(Prob)
end

recovery!(alpha, (;S, E, I, R), cols; γ, _...,) = begin
    rate = γ*I
    ellI = ell(cols,Infec)
    alpha[1] = @indicator(I > ellI, rate*(1-ellI/I))
    rate-sum(alpha)
end

waning!(alpha, (;S, E, I, R), cols; ω, _...,) = begin
    alpha[1] = ω*R
    zero(Prob)
end

sampling((;S, E, I, R), cols; ψ, χ, _...,) = begin
    (ψ+χ)*I
end

event_rates!(alpha, state, cols; kwargs...) = begin
    decay = zero(Prob)
    decay += transmission!(@view(alpha[1]), state, cols; kwargs...,)
    decay += progression!(@view(alpha[2]), state, cols; kwargs...,)
    decay += recovery!(@view(alpha[3]), state, cols; kwargs...,)
    decay += waning!(@view(alpha[4]), state, cols; kwargs...,)
    decay += sampling(state, cols; kwargs...,)
    decay
end

deme_occupancy(;E, I, _...,) = begin
    ;E, I
end

live_condition((;S, E, I, R), cols) = begin
    ellE, ellI = ell(cols)
    I ≥ ellI && E ≥ ellE
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
        ll = log(p)
        plant!(cols,D(i),node.chillins[1])
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, state
end

terminal_sample!(node, cols, (;S, E, I, R); ψ, χ, _...) = begin
    if node.parlin ∈ cols[Infec]
        k,_,p = rcateg([ψ, χ],true)
        ll = -log(p)
        ellE, ellI = chop!(cols,Infec,node.parlin)
        if k==0                     # both ψ and χ are zero
            live = false
            ll = Prob(-Inf)
        elseif k==1                 # non-destructive sample
            live = true
            ll += log(ψ*(I-ellI))
        elseif k==2                 # destructive sample
            live = true
            ll += log(χ*I)
            I -= 1
        end
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, (;S, E, I, R)
end

inline_sample!(node, cols, state; ψ, _...) = begin
    if node.parlin ∈ cols[Infec]
        live = true
        ll = log(ψ)
        chop!(cols,Infec,node.parlin,Infec,node.chillins[1])
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, state
end

branch_options(x::AbstractArray{Prob,2}) = begin
    [x[1,1]*x[2,2], x[1,2]*x[2,1]]
end


singular_sample!(node, cols, state; kwargs...) = begin
    if length(node.chillins) == 0
        terminal_sample!(node, cols, state; kwargs...)
    elseif length(node.chillins) == 1
        inline_sample!(node, cols, state; kwargs...)
    end
end

singular_branch!(node, cols, (;S, E, I, R); β, pop, _...) = begin
    if node.parlin ∈ cols[Infec]
        live = true
        ll = log(β*S*I/pop)
        k, _, p = rcateg(branch_options(node.present), true)
        ll -= log(p)
        @assert k ≠ 0
        if k==1
            fork!(cols,Infec,node.parlin,(Expos,Infec),node.chillins)
        else
            fork!(cols,Infec,node.parlin,(Infec,Expos),node.chillins)
        end
        if S > 0
            S -= 1
        end
        E += 1
        ll -= log(E*I)
    else
        live = false
        ll = Prob(-Inf)
    end
    ;live, ll, (;S, E, I, R)
end

singular_part!(
    cols, state, guide, n;
    kwargs...,
) = begin
    ells = ell(cols)
    live = live_condition(state, cols)
    if live
        node = guide[n]
        if node.type==Root
            singular_root!(node, cols, state; kwargs...)
        elseif node.type==Sample
            singular_sample!(node, cols, state; kwargs...)
        elseif node.type==Node
            singular_branch!(node, cols, state; kwargs...)
        end
    else
        ll = Prob(-Inf)
        ;live, ll, state
    end
end

regular_transmission!(
    t, guide, node, cols, (;S, E, I, R);
    kwargs...,
) = begin
    b,p = choose_branch(t,guide,node,I,cols,Infec,Expos)
    ll = -log(p)
    S -= 1
    E += 1
    if b == 0
        ellE = ell(cols,Expos)
        ll += log(1-ellE/E)
    else
        ellE, ellI = swap!(cols,Infec,Expos,b)
        ll += log(1-ellI/I)-log(E)
    end
    ll, (;S, E, I, R)
end

regular_progression!(
    t, guide, node, cols, (;S, E, I, R);
    kwargs...,
) = begin
    b,p = choose_branch(t,guide,node,E,cols,Expos,Infec)
    ll = -log(p)
    E -= 1
    I += 1
    if b == 0
        ellI = ell(cols,Infec)
        ll += log(1-ellI/I)
    else
        swap!(cols,Expos,Infec,b)
        ll -= log(I)
    end
    ll, (;S, E, I, R)
end

regular_recovery!(
    t, guide, node, cols, (;S, E, I, R);
    kwargs...,
) = begin
    ellI = ell(cols,Infec)
    ll = -log(1-ellI/I)
    I -= 1
    R += 1
    ll, (;S, E, I, R)
end

regular_waning!(
    t, guide, node, cols, (;S, E, I, R);
    kwargs...,
) = begin
    R -= 1
    S += 1
    zero(Prob), (;S, E, I, R)
end

regular_part!(
    cols, state, guide, n, t, tf;
    kwargs...,
) = begin
    node = guide[n]
    alpha = similar(Vector{Prob},4)
    step::Time = zero(Time)
    decay::Prob = zero(Prob)
    ll::Prob = zero(Prob)
    while t < tf
        decay = event_rates!(alpha, state, cols; kwargs...,)
        k, s = rcateg(alpha)
        step = -log(rand())/s
        if k > 0 && t+step < tf
            if k==1
                ll1, state = regular_transmission!(t,guide,n,cols,state)
            elseif k==2
                ll1, state = regular_progression!(t,guide,n,cols,state)
            elseif k==3
                ll1, state = regular_recovery!(t,guide,n,cols,state)
            elseif k==4
                ll1, state = regular_waning!(t,guide,n,cols,state)
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

seir_rinit(; S0, E0, I0, R0, pop, _...) = begin
    m = pop/(S0+E0+I0+R0)
    (
        S=round(Int64, m*Float64(S0)),
        E=round(Int64, m*Float64(E0)),
        I=round(Int64, m*Float64(I0)),
        R=round(Int64, m*Float64(R0)),
    )
end

"""
    filter_pomp(g, m; β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02,
         χ = 0.0, pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08)

Constructs a pomp object based on the genealogy `g` and finite-state
Markov guiding process `m`.
"""
filter_pomp(
    gen::Genealogy,
    m::FSMarkovProc;
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
    pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    check(gen)
    guidegen = guide(gen,m,knowledge!)
    pomp(
        params = (
            β = Float64(β), σ = Float64(σ), γ = Float64(γ),
            ω = Float64(ω), ψ = Float64(ψ), χ = Float64(χ),
            pop = Float64(pop),
            S0 = Float64(S0), E0 = Float64(E0),
            I0 = Float64(I0), R0 = Float64(R0),
        ),
        t0 = timezero(guidegen),
        times = times(guidegen),
        userdata = (guide = guidegen,),
        init_state = (
            node=one(Name),
            ll=zero(Prob),
            live=true,
            cols=Coloring(Demes),
            state=(S=Float64(0), E=Float64(0), I=Float64(0), R=Float64(0)),
        ),
        rinit = function (; kwargs...)
            (
                node = one(Name),
                ll = zero(Prob),
                live = true,
                cols = Coloring(Demes),
                state=seir_rinit(;kwargs...),
            )
        end,
        rprocess = onestep(
            function (;t, dt, guide, node, ll, live, cols, state, kwargs...,)
                if live
                    cols = copy(cols)
                    live, ll, state = singular_part!(
                        cols, state, guide, node;
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
