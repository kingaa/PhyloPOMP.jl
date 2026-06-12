import Base: eachindex, length, getindex, eachindex, show

"""
    GuideNode{F,D}

A `GuideNode` contains the information needed to guide a phylodynamic
Monte Carlo filter on a particular interval between genealogical events.
`F` is an `AbstractFloat` type and `D` is the enumeration of the demes
(see [`@demes`](@ref)).

A `GuideNode` contains:
- the time interval to which it pertains
- the type of event (`NodeType`) at the left endpoint of the interval
- the name of the parent lineage (ancestral to the genealogical node at the left endpoint)
- the names of the child lineages (immediately descended from the parent node)
- the names of all lineages present on the interval
- the guide probability for each child lineage at the left interval endpoint
- the target probabilities for all lineages at the right interval endpoint.
"""
struct GuideNode{F<:AbstractFloat,D<:Enum}
    "type of node at left endpoint of the interval"
    type::NodeType
    "left end of the time interval"
    tbeg::Time
    "right end of the time interval"
    tend::Time
    "identity of the parent lineage"
    parlin::Name
    "identities of child lineage(s)"
    chillins::Vector{Name}
    "identities of all lineages present on the interval"
    alllins::Vector{Name}
    "Dict that maps lineages to indices of the various matrices"
    linmap::Dict{Name,Name}
    "Matrix giving the probabilities of finding a lineage in a deme"
    present::Matrix{F}
    """
    Matrix giving the target probabilities (of finding a lineage in
    a given deme at the end of the interval).
    """
    target::Matrix{F}
    "`target` rexpressed in the eigenbasis of the guide process"
    dtarget::Matrix{F} ## FIXME: inelegant to store redundant information
end

"""
    Guide{F,D}

This type holds a "filter-guide" based on a reverse-time finite-state Markov process.
Construct it with a call to [`guide`](@ref).
`F` is an `AbstractFloat` type and `D` is the enumeration of the demes.
It contains the [`FSMarkovProc`](@ref) object encoding the Markov process and a sequence of [`GuideNode`](@ref) objects.
"""
struct Guide{F<:AbstractFloat,D<:Enum}
    "finite-state Markov process"
    process::FSMarkovProc{F,D}
    "vector of nodes, one for each interval"
    nodes::Vector{GuideNode{F,D}}
    "total number of sample lineages"
    nsample::Size
    Guide(
        g::Genealogy{E},
        m::FSMarkovProc{F,D},
        knowledge!::Function = known_deme!,
    ) where {F<:AbstractFloat,D<:Enum,E} = begin
        probs = zeros(F,length(instances(D)),g.nsample)
        nodes = Vector{GuideNode{F,D}}(undef,length(g))
        tend::Time = g.time
        lins::Set{Name} = Set{Name}()
        for n ∈ reverse(eachindex(g))
            ells = collect(lins)
            parlin = g[n].lineage
            chillins = map(x->g[x].lineage, g[n].children)
            target = probs[:,ells]
            ## FIXME: inelegant to store redundant information
            dtarget = m.right_trans * target
            probs[:,ells] = forward_action(m,tend-g[n].slate,probs[:,ells])
            present = probs[:,chillins]
            linmap = Dict(zip(ells,eachindex(ells)))
            nodes[n] = GuideNode{F,D}(
                g[n].type,
                g[n].slate,tend,
                parlin,chillins,ells,linmap,
                present,target,dtarget,
            )
            known = knowledge!(
                @view(probs[:,parlin]);
                deme=g[n].deme,
                type=g[n].type,
                time=g[n].slate,
            )
            if known
                @inbounds probs[:,parlin] = probs[:,parlin]./sum(probs[:,parlin])
            else
                @inbounds assimil!(
                    @view(probs[:,parlin]),
                    @view(probs[:,chillins]),
                    statdist(m)
                )
            end
            tend = g[n].slate
            union!(setdiff!(lins,chillins),parlin)
        end
        new{F,D}(m,nodes,g.nsample)
    end
end

"""
    known_deme!(v; deme, type, time)

The default behavior of [`guide`](@ref) is to assume the deme is known
if and only if it is fixed in the genealogy.  This function implements
that default.
"""
known_deme!(
    v::AbstractVector{F};
    deme, type, time,
) where F = begin
    if (ismissing(deme))
        false
    else
        demekron!(v,deme)
        true
    end
end

"""
    guide(g, m, knowledge!)

- `g` is a [`Genealogy`](@ref).
- `m`: is an [`FSMarkovProc`](@ref) encoding the guiding finite-state Markov process.
  Construct it with a call to [`fsmarkov`](@ref).
- `knowledge!` is a function with signature `knowledge!(v; deme, type, time)`.
  It should overwrite the vector `v` with probabilities of being in each
  possible deme, according to the given `deme`, `type`, and `time`.
  By default, `knowledge! = known_deme!` (see [`known_deme!`](@ref)).
"""
guide(args...) = Guide(args...)

"""
    demekron([F::Type], d)

Kronecker delta for deme `d`. By default, `F = Float64`.
"""
demekron(F::Type, d::D,) where {D <: Enum} = begin
    [i == d ? one(F) : zero(F) for i ∈ instances(D)]
end

demekron(d::D,) where {D <: Enum} = demekron(Prob,d)

"""
    demekron!(v, d)

Kronecker delta for deme `d`.  Replaces the entries of `v` with zeros in all places except the one corresponding to deme `d`.
"""
demekron!(v::AbstractVector{F}, d::D,) where {F,D <: Enum} = begin
    @assert length(v)==length(instances(D))
    for k ∈ eachindex(v)
        v[k] = k == Int(d) ? one(F) : zero(F)
    end
    nothing
end

"""
    relhaz(t, g, n, i, j, lins)

Compute the relative reverse-time hazards associated with filter guide `g` at time `t` relative
to the target at guide-node `n`.  Only columns denoted by the indices in `lins` are considered.

This returns a vector of hazards for the `i ⟶ j` transition relative to the hazard at stationarity.
"""
relhaz(
    t::Time,
    g::Guide{F,D},
    n::Integer,
    i::D,
    j::D,
    lins::Union{Integer,AbstractVector{<:Integer}},
) where {F,D} = begin
    @assert g[n].tbeg ≤ t < g[n].tend "t ∉ [$(g[n].tbeg),$(g[n].tend))"
    p = relhaz_action(
        g.process,
        g[n].tend-t,
        @view(g[n].dtarget[:,lins]),
    )
    @views p[Int(j),:]./p[Int(i),:]
end

"""
    relhaz(t, g, n, i, j, cols)

In this form, `cols` is a `Coloring`.
"""
relhaz(
    t::Time,
    g::Guide{F,D},
    n::Integer,
    i::D,
    j::D,
    cols::Coloring{D}
) where {F,D} = begin
    lins = [g[n].linmap[k] for k ∈ cols[i]]
    relhaz(t,g,n,i,j,lins)
end

choose_branch(
    t::Time,
    guide::Guide{F,D},
    node::Integer,
    n,
    cols::Coloring{D},
    i::D, j::D,
) where {F,D} = begin
    lins = [guide[node].linmap[k] for k ∈ cols[i]]
    p = relhaz(t,guide,node,i,j,lins)
    s1 = F(n-ell(cols,i))
    s2 = sum(p)
    s = s1+s2
    r = s*rand(F)
    if r > s1
        r -= s1
        k::Name = 1
        while r > p[k]
            r -= p[k]
            k += 1
        end
        guide[node].alllins[lins[k]], p[k]/s
    else
        zero(Name), s1/s
    end
end

assimil!(
    parprob::AbstractVector{F},
    chilprobs::AbstractMatrix{F},
    pi::AbstractVector{F},
) where F = begin
    @assert size(parprob,1)==size(chilprobs,1)
    if isempty(chilprobs)
        parprob[:] = pi
    else
        for i ∈ axes(chilprobs,1)
            for j ∈ axes(chilprobs,2)
                @inbounds chilprobs[i,j] /= pi[i]
            end
        end
        s = zero(F)
        for i ∈ axes(chilprobs,1)
            @inbounds parprob[i] = pi[i]*prod(chilprobs[i,:])
            @inbounds s += parprob[i]
        end
        for i ∈ eachindex(parprob)
            @inbounds parprob[i] /= s
        end
    end
    nothing
end

Base.getindex(g::Guide, i::Integer) = g.nodes[i]
Base.getindex(g::Guide, i::AbstractVector{<:Integer}) = g.nodes[i]
Base.length(g::Guide) = length(g.nodes)
Base.eachindex(g::Guide) = eachindex(g.nodes)

timezero(g::Guide) = g.nodes[1].tbeg
times(g::Guide) = map(n->n.tend,g.nodes)

Base.show(
    io::IO,
    g::Union{Guide,GuideNode};
    kwargs...,
) = begin
    print(io,pretty_string(g;kwargs...))
end

pretty_string(g::Guide; sigdigits = 4) = begin
    nodes = map(eachindex(g)) do i
        "  $i: "*pretty_string(g[i],sigdigits=sigdigits)
    end
    "<guide:\n" * join(nodes,'\n') * ">"
end

pretty_string(n::GuideNode; sigdigits = 4) = begin
    t1 = round(n.tbeg,sigdigits=sigdigits)
    t2 = round(n.tend,sigdigits=sigdigits)
    chillins = map(eachindex(n.chillins)) do i
        ell = n.chillins[i]
        prob = round.(n.present[:,i],sigdigits=sigdigits)
        "$ell=>$prob"
    end
    targs = map(eachindex(n.alllins)) do i
        ell = n.alllins[i]
        prob = round.(n.target[:,i],sigdigits=sigdigits)
        "$ell=>$prob"
    end
    "<$(lowercase(String(Symbol(n.type)))): t ∈ [$t1,$t2] " *
        "parlin=$(n.parlin)" *
        " children:{" * join(chillins,',') *
        "} targets:{" * join(targs,',') * "}>"
end
