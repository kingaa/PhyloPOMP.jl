"""
    NaiveSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using so-called naïve proposals. These proposals are
non-anticipatory. In particular, population events occurring at any
time result in color changes with a probability proportional to the
fraction of appropriate lineages represented in the genealogy at that
time.
"""
module NaiveSEIR

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time

@demes Demes Expos Infec
using .Demes: Expos, Infec, DemeSet

include("seir_trees.jl")

## gives the occupancy of each deme
deme_occupancy(;E, I, _...,) = begin
    ;E, I
end

## Is the population state compatible with the genealogy?
live_condition((;S, E, I, R), cols) = begin
    ellE, ellI = ell(cols)
    I ≥ ellI && E ≥ ellE
end

seir_rinit(; S0, E0, I0, R0, pop, _...) =
    barycentric(Int64,(S=S0,E=E0,I=I0,R=R0),pop)

singular_root!(node, cols, state, t; kwargs...,) = begin
    n = deme_occupancy(;state...)
    ells = ell(cols)
    i, _, p = rcateg(n.-ells, true)
    ll = -log(p)
    if i > 0
        live = true
        plant!(cols, DemeSet(i), node.lineage)
    else
        live = false
        ll = Prob(-Inf)
    end
    live, ll, state
end

terminal_sample!(node, cols, (;S, E, I, R), t; ψ, χ, _...) = begin
    boost::Prob = 1
    live = true
    if node.lineage ∈ cols[Infec]
        k,_,p = rcateg([ψ, χ], true)
        _, ellI = chop!(cols, Infec, node.lineage)
        if k==0
            live = false
        elseif k==1         # non-destructive sample
            boost = ψ*(I-ellI)/p
        elseif k==2         # destructive sample
            boost = χ*I/p
            I -= 1
        end
    else
        live = false
        boost = 0
    end
    live, log(boost), (;S, E, I, R)
end

inline_sample!(geneal, node, cols, state, t; ψ, _...) = begin
    boost::Prob = 1
    if node.lineage ∈ cols[Infec]
        live = true
        boost = ψ
        chillin = geneal[node.children[1]].lineage
        chop!(cols, Infec, node.lineage, Infec, chillin)
    else
        live = false
        boost = 0
    end
    live, log(boost), state
end

singular_sample!(geneal, node, cols, state, t; kwargs...) = begin
    if length(node.children) == 0
        terminal_sample!(node, cols, state, t; kwargs...)
    elseif length(node.children) == 1
        inline_sample!(geneal, node, cols, state, t; kwargs...)
    else
        @assert false "impossible" # COV_EXCL_LINE
    end
end

singular_branch!(geneal, node, cols, (;S, E, I, R), t; β, pop, _...) = begin
    boost::Prob = 1
    if node.lineage ∈ cols[Infec] && S > 0
        live = true
        chillins = map(node.children) do i
            geneal[i].lineage
        end
        k, _, p = rcateg([1, 1], true)
        boost = β*S*I/pop/p
        @assert k ≠ 0
        if k==1
            ellE, ellI = fork!(cols, Infec, node.lineage, (Expos,Infec), chillins)
        else
            ellE, ellI = fork!(cols, Infec, node.lineage, (Infec,Expos), chillins)
        end
        S -= 1
        E += 1
        boost *= @binomratio((E,I),(ellE,ellI),(1,1),(1,1))
    else
        live = false
        boost = 0
    end
    live, log(boost), (;S, E, I, R)
end

singular_part!(
    cols, state, t, geneal, n;
    kwargs...,
) = begin
    live = live_condition(state, cols)
    ll::Prob = 0
    if live
        node = geneal[n]
        if node.type==Root
            live, ll, state = singular_root!(node, cols, state, t; kwargs...)
        elseif node.type==Sample
            live, ll, state = singular_sample!(geneal, node, cols, state, t; kwargs...)
        elseif node.type==Node
            live, ll, state = singular_branch!(geneal, node, cols, state, t; kwargs...)
        else
            @assert false # COV_EXCL_LINE
        end
    else
        ll = Prob(-Inf)
    end
    live, ll, state
end

rate_transmission!(alpha, pi, t, (;S, E, I, R), cols; β, pop, _...,) = begin
    ellI = ell(cols,Infec)
    rate = β*S*I/pop
    pi[1] = @indicator(I > 0, (I+1-ellI)/(I+1))
    pi[2] = @indicator(I > 0, ellI/(I+1))
    alpha .= rate.*pi
    rate-sum(alpha)
end

rate_progression!(alpha, pi, t, (;S, E, I, R), cols; σ, _...,) = begin
    rate = σ*E
    ellE = ell(cols,Expos)
    pi[1] = @indicator(E > 0, 1-ellE/E)
    pi[2] = @indicator(E > 0, ellE/E)
    alpha .= rate.*pi
    rate-sum(alpha)
end

rate_recovery!(alpha, pi, t, (;S, E, I, R), cols; γ, _...,) = begin
    rate = γ*I
    ellI = ell(cols,Infec)
    pi[1] = @indicator(I > ellI, 1-ellI/I)
    alpha .= rate.*pi
    rate-sum(alpha)
end

rate_waning!(alpha, pi, t, (;S, E, I, R), cols; ω, _...,) = begin
    rate = ω*R
    pi[1] = 1
    alpha .= rate.*pi
    rate-sum(alpha)
end

rate_sampling(t, (;S, E, I, R), cols; ψ, χ, _...,) = begin
    rate = (ψ+χ)*I
    rate
end

event_rates!(alpha, pi, t, state, cols; kwargs...) = begin
    decay = zero(Prob)
    decay += rate_transmission!(@view(alpha[1:2]), @view(pi[1:2]), t, state, cols; kwargs...,)
    decay += rate_progression!(@view(alpha[3:4]), @view(pi[3:4]), t, state, cols; kwargs...,)
    decay += rate_recovery!(@view(alpha[5]), @view(pi[5]), t, state, cols; kwargs...,)
    decay += rate_waning!(@view(alpha[6]), @view(pi[6]), t, state, cols; kwargs...,)
    decay += rate_sampling(t, state, cols; kwargs...,)
    decay
end

move_transmission_off!(t, (;S, E, I, R), cols) = begin
    S -= 1
    E += 1
    boost = @binomratio(E,ell(cols,Expos),1,0)
    log(boost), (;S, E, I, R)
end

move_transmission_on!(t, (;S, E, I, R), cols) = begin
    b, p = rand(cols[Infec]), inv(ell(cols,Infec))
    S -= 1
    E += 1
    ellE, ellI = swap!(cols,Infec,Expos,b)
    boost = @binomratio((E,I),(ellE,ellI),(1,1),(1,0))/p
    log(boost), (;S, E, I, R)
end

move_progression_off!(t, (;S, E, I, R), cols) = begin
    E -= 1
    I += 1
    boost = @binomratio(I,ell(cols,Infec),1,0)
    log(boost), (;S, E, I, R)
end

move_progression_on!(t, (;S, E, I, R), cols) = begin
    b, p = rand(cols[Expos]), inv(ell(cols,Expos))
    ellE, ellI = swap!(cols,Expos,Infec,b)
    E -= 1
    I += 1
    boost = @binomratio(I,ellI,1,1)/p
    log(boost), (;S, E, I, R)
end

move_recovery!(t, (;S, E, I, R), cols) = begin
    I -= 1
    R += 1
    boost = 1.0
    log(boost), (;S, E, I, R)
end

move_waning!(t, (;S, E, I, R), cols) = begin
    R -= 1
    S += 1
    boost = 1.0
    log(boost), (;S, E, I, R)
end

regular_move!(k, t, state, cols) = begin
    if k==1
        move_transmission_off!(t, state, cols)
    elseif k==2
        move_transmission_on!(t, state, cols)
    elseif k==3
        move_progression_off!(t, state, cols)
    elseif k==4
        move_progression_on!(t, state, cols)
    elseif k==5
        move_recovery!(t, state, cols)
    elseif k==6
        move_waning!(t, state, cols)
    else
        @assert false "impossible" # COV_EXCL_LINE
    end
end

regular_part!(
    cols, state, t, tf;
    kwargs...,
) = begin
    if t < tf
        alpha = similar(Vector{Prob}, 6)
        pi = similar(Vector{Prob}, 6)
        step::Time = 0
        decay::Prob = 0
        ll::Prob = 0
        while t < tf
            decay = event_rates!(alpha, pi, t, state, cols; kwargs...,)
            k, s = rcateg(alpha)
            step = -log(rand())/s
            if k > 0 && t+step < tf
                ll1, state = regular_move!(k, t, state, cols)
                ll += ll1 - decay*step - log(pi[k])
                t += step
            else
                step = tf - t
                ll -= decay*step
                break
            end
        end
    end
    ll, state
end

"""
    filter_pomp(gen; β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
         pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08)

Constructs a pomp object based on the genealogy `gen` and implementing
the naïve genealogical filter.
"""
filter_pomp(
    gen::Genealogy;
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
    pop = 100,
    S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    check(gen)
    pomp(
        params = map(
            Float64,
            (;β,σ,γ,ω,ψ,χ,pop,S0,E0,I0,R0)
        ),
        t0 = timezero(gen),
        times = times(gen),
        userdata = (geneal = gen,),
        init_state = (
            node=one(Name),
            ll=zero(Prob),
            live=true,
            cols=Coloring(Demes),
            state=(S=Int64(0), E=Int64(0), I=Int64(0), R=Int64(0)),
        ),
        rinit = function (;kwargs...)
            (
                node = one(Name),
                ll = zero(Prob),
                live = true,
                cols = Coloring(Demes),
                state = seir_rinit(;kwargs...),
            )
        end,
        rprocess = onestep(
            function (;node, ll, live, cols, state, geneal, t, dt, kwargs...)
                if live
                    cols = copy(cols)
                    live, ll, state = singular_part!(cols, state, t, geneal, node; kwargs...)
                else
                    ll = Prob(-Inf)
                end
                if live && dt > 0
                    ll1, state = regular_part!(cols, state, t, t+dt; kwargs...)
                    ll += ll1
                end
                (;node = node+1, ll, live, cols, state)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
    )
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

end
