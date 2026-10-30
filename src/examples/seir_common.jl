using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

## Construct an enumeration of the demes for this model.
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

## Singular-part section.
##
## The following functions deal with the cases of having a root,
## sample, or branch at an observed event time. They are called in
## turn by the `singular_part!` function defined at the bottom of this
## section.

singular_root!(
    node::GuideNode{F,N,D}, cols, state, t;
    kwargs...,
) where {F, N, D <: Enum} = begin
    n = deme_occupancy(;state...)
    ells = ell(cols)
    i, _, p = rcateg(node.present[:,1].*(n.-ells), true)
    q::Prob = 0
    if i > 0
        live = true
        q = inv(p)
        plant!(cols, D(i), node.chillins[1])
    else
        live = false
        q = 0
    end
    live, log(q), state
end

terminal_sample!(node, cols, (;S, E, I, R), t; ψ, χ, _...) = begin
    boost::Prob = 1
    if node.parlin ∈ cols[Infec]
        k,_,p = rcateg([ψ, χ],true)
        _, ellI = chop!(cols, Infec, node.parlin)
        if k==0                     # both ψ and χ are zero
            live = false
            boost = 0
        elseif k==1                 # non-destructive sample
            live = true
            boost = ψ*(I-ellI)/p
        elseif k==2                 # destructive sample
            live = true
            boost = χ*I/p
            I -= 1
        end
    else
        live = false
        boost = 0
    end
    live, log(boost), (;S, E, I, R)
end

inline_sample!(node, cols, state, t; ψ, _...) = begin
    boost::Prob = 1
    if node.parlin ∈ cols[Infec]
        live = true
        boost = ψ
        chop!(cols, Infec, node.parlin, Infec, node.chillins[1])
    else
        live = false
        boost = 0
    end
    live, log(boost), state
end

singular_sample!(node, cols, state, t; kwargs...) = begin
    if length(node.chillins) == 0
        terminal_sample!(node, cols, state, t; kwargs...)
    elseif length(node.chillins) == 1
        inline_sample!(node, cols, state, t; kwargs...)
    else
        @assert false "impossible" # COV_EXCL_LINE
    end
end

branch_options(x::AbstractArray{Prob,2}) = begin
    [x[1,1]*x[2,2], x[1,2]*x[2,1]]
end

singular_branch!(node, cols, (;S, E, I, R), t; β, pop, _...) = begin
    boost::Prob = 1
    if node.parlin ∈ cols[Infec] && S > 0
        live = true
        k, _, p = rcateg(branch_options(node.present), true)
        boost = β*S*I/pop/p
        @assert k ≠ 0
        if k==1
            ellE, ellI = fork!(cols,Infec,node.parlin,(Expos,Infec),node.chillins)
        else
            ellE, ellI = fork!(cols,Infec,node.parlin,(Infec,Expos),node.chillins)
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
    cols, state, guide, n, t;
    kwargs...,
) = begin
    live = live_condition(state, cols)
    if live
        node = guide[n]
        if node.type==Root
            singular_root!(node, cols, state, t; kwargs...)
        elseif node.type==Sample
            singular_sample!(node, cols, state, t; kwargs...)
        elseif node.type==Node
            singular_branch!(node, cols, state, t; kwargs...)
        else
            @assert false "impossible" # COV_EXCL_LINE
        end
    else
        ll = Prob(-Inf)
        live, ll, state
    end
end

"""
    seir_rinit(;S0, E0, I0, R0, pop, _...)

Initializes the state of the population process.  The parameters `S0`,
`E0`, `I0`, `R0` are taken to be barycentric coordinates for a
position on the lattice simplex with sum `pop`.  That is, the returned
population state is a tuple `(S,E,I,R)` of integers, where `S + E + I
+ R = pop`.
"""
seir_rinit(; S0, E0, I0, R0, pop, _...) =
    barycentric(Int64,(S=S0,E=E0,I=I0,R=R0),pop)

"""
    filter_pomp(g, m; β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02,
         χ = 0.0, pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08)

Constructs a pomp object based on the genealogy `g` and finite-state
Markov guiding process `m`. The keyword arguments are model
parameters.
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
        params = map(
            Float64,
            (;β,σ,γ,ω,ψ,χ,pop,S0,E0,I0,R0)
        ),
        t0 = timezero(guidegen),
        times = times(guidegen),
        userdata = (guide = guidegen,),
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
            function (;t, dt, guide, node, ll, live, cols, state, kwargs...,)
                if live
                    cols = copy(cols)
                    live, ll, state = singular_part!(
                        cols, state, guide, node, t;
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

## Check the genealogy for incompatibility with the model.
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
