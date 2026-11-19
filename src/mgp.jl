using PartiallyObservedMarkovProcesses: pomp

const MgpDict = Dict{Symbol,Any}
const FuncMap = Dict{Tuple{Int,Int},Vector{Symbol}}

include("extract_symbols.jl")
include("parse_mgp.jl")
include("switchblock.jl")
include("filter_singular.jl")
include("filter_regular.jl")

macro mgp_soft_filter(args...)
    expr = soft_filter_def(args...)
    esc(:(Base.@__doc__ $expr))
end

soft_filter_def(args::Expr...) = begin
    dict = MgpDict(
        :demes => :((default,)),
        :state => :(),
        :eventmap => FuncMap(),
        :movemap => FuncMap(),
        :logdmeasure => :((;_...) -> zero(Prob)),
        :knowledge => quote
            function(v; deme, type, time)
                if ismissing(deme)
                    false
                else
                    demekron!(v,deme)
                    true
                end
            end
        end,
        :decls => Expr[
            :(using PhyloPOMP),
            :(using PhyloPOMP: Name, Time, Prob, Root, Node, Sample),
        ],
    )
    parse_mgp_spec!(dict, args...)
    modname = pop!(dict,:name,Symbol(:MGP_,gensym()))
    process_demes!(dict)
    process_jumps!(dict)
    make_knowledge!(dict)
    make_rinitfn!(dict)
    make_init_state!(dict)
    make_logdmeasure!(dict)
    make_rprocessfn!(dict)
    decls = peelblock.(dict[:decls])
    mdoc = Expr(:block, Expr(:meta, :doc), modname)
    Expr(:toplevel,
         Expr(:module, true, modname,
              Expr(:block, decls...)),
         mdoc, modname)
end

process_demes!(dict) = begin
    @assert haskey(dict,:demes) "demes not specified"
    @assert dict[:demes].head == :tuple "demes should be specified as a tuple"
    dict[:demes] = dict[:demes].args
    state = dict[:state]
    dict[:demeoccfn] = fname = :deme_occupancy
    @assert haskey(dict,:deme_occupancy) "deme_occupancy not specified"
    demeocc = pop!(dict,:deme_occupancy)
    defn = quote
        const $fname = function($state,)
            $demeocc
        end
    end
    push!(dict[:decls],:(@eval @demes(DemeSet, $(dict[:demes]...))))
    push!(dict[:decls],defn)
    nothing
end

process_jumps!(dict) = begin
    jumps = get(dict, :jumps, MgpDict())
    names = keys(jumps)
    for name ∈ names
        jump = jumps[name]
        jump[:name] = name
        process_jump!(jump, dict)
    end
    nothing
end

process_jump!(jump, dict) = begin
    make_ratefn!(jump, dict)
    make_popmovefn!(jump, dict)
    process_moves!(jump, dict)
end

## construct one rate-function per jump and store in the dictionary
make_ratefn!(jump, dict) = begin
    @assert haskey(jump,:rate) "no rate specified for jump $name"
    rate = pop!(jump,:rate)
    state = dict[:state]
    vars = setdiff(extract_symbols(rate),extract_symbols(state),(:t,))
    fname = Symbol(:rate_,jump[:name])
    jump[:ratefn] = fname
    defn = quote
        const $fname = function(t, $state; $(vars...), _...)
            $rate
        end
    end
    push!(dict[:decls],defn)
    nothing
end

## construct one popmove-function per jump and store in the dictionary
make_popmovefn!(jump, dict) = begin
    fname = Symbol(:popmove_,jump[:name])
    if haskey(jump,:popmove)
        move = (pop!(jump,:popmove),)
    else
        move = ()
    end
    state = dict[:state]
    vars = setdiff(extract_symbols(move),extract_symbols(state),(:t,))
    fname = Symbol(:popmove_,jump[:name])
    jump[:popmovefn] = fname
    defn = quote
        const $fname = function(t, $state; $(vars...), _...)
            $(move...)
            $state
        end
    end
    push!(dict[:decls],defn)
    nothing
end

process_moves!(jump, dict) = begin
    @assert haskey(jump,:moves)
    moves = jump[:moves]
    names = keys(moves)
    for name ∈ names
        move = moves[name]
        move[:name] = Symbol(jump[:name],:_,name)
        process_move!(move, jump, dict)
    end
    nothing
end

process_move!(move, jump, dict) = begin
    @assert haskey(move,:type)
    make_propfn!(move,jump,dict)
    type = pop!(move,:type)
    if type == :id
        jump[:has_reg] = :true
        jump[:regpropfn] = pop!(move,:propfn)
        make_id_move!(move, jump, dict)
    elseif type == :swap
        jump[:has_reg] = :true
        jump[:regpropfn] = pop!(move,:propfn)
        make_swap_move!(move, jump, dict)
    elseif type == :fork
        jump[:has_sing] = :true
        make_fork_propnmove!(move, jump, dict)
    elseif type == :chop
        jump[:has_sing] = :true
        make_chop_propnmove!(move, jump, dict)
    else
        error("unknown type $type in move $(move[:name]) for jump $(jump[:name]).")
    end
    nothing
end

make_propfn!(move, jump, dict) = begin
    if haskey(move,:pi)
        prop = pop!(move,:pi)
    else
        prop = 1
    end
    state = dict[:state]
    enns = dict[:enns]
    ells = dict[:ells]
    vars = setdiff(
        extract_symbols(prop),
        extract_symbols(state),extract_symbols(enns),extract_symbols(ells),
    )
    fname = Symbol(:pi_,move[:name])
    move[:propfn] = fname
    defn = quote
        const $fname = function($state, $enns, $ells; $(vars...), _...)
            $prop
        end
    end
    push!(dict[:decls],defn)
    nothing
end

## construct the boost function for a move
make_boostfn!(boost, move, dict, name) = begin
    fname = Symbol(:boost_,move[:name],:_,name)
    enns = dict[:enns]
    ells = dict[:ells]
    defn = quote
        const $fname = function($enns, $ells)
            $boost
        end
    end
    push!(dict[:decls],defn)
    fname
end

make_init_state!(dict) = begin
    @assert haskey(dict,:init_state) "init_state not specified"
    init_state = eval(pop!(dict,:init_state))
    defn = quote
        const init_state = (
            node=one(Name),
            ll=zero(Prob),
            live=true,
            cols=Coloring(DemeSet),
            state=$init_state,
        )
    end
    push!(dict[:decls],defn)
    nothing
end

make_rinitfn!(dict) = begin
    @assert haskey(dict,:rinit) "rinit not specified"
    rinit = pop!(dict,:rinit)
    push!(dict[:decls],:(const user_rinit = $rinit))
    push!(
        dict[:decls],
        quote
            const rinitfn = function(;kwargs...)
                (
                    node = one(Name),
                    ll = zero(Prob),
                    live = true,
                    cols = Coloring(DemeSet),
                    state = user_rinit(;kwargs...),
                )
            end
        end
    )
    nothing
end

make_rprocessfn!(dict) = begin
    make_regular_part!(dict)
    make_singular_part!(dict)
    delete!(dict,:jumps)
    defn = quote
        const rprocessfn = function(;t, dt, guide, node, ll, live, cols, state, kwargs...,)
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
                local ll1
                ll1, state = regular_part!(
                    cols, state, guide, node, t, t+dt;
                    kwargs...,
                )
                ll += ll1
            end
            (;node = node+1, ll, live, cols, state)
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_logdmeasure!(dict) = begin
    @assert haskey(dict,:logdmeasure) "logdmeasure not specified"
    f = pop!(dict,:logdmeasure)
    defn = :(const user_logdmeas = $f)
    push!(dict[:decls], defn)
    push!(
        dict[:decls],
        quote
            const logdmeasurefn = function(; node, ll, live, cols, state, kwargs...)
                if live
                    ll + user_logdmeas(;state..., kwargs...)
                else
                    Prob(-Inf)
                end
            end
        end
    )
    nothing
end

make_knowledge!(dict) = begin
    @assert haskey(dict,:knowledge) "knowledge not specified"
    know = pop!(dict,:knowledge)
    push!(dict[:decls],:(const knowledge! = $know))
    nothing
end

filter_pomp(
    x::Module,
    g::Guide;
    kwargs...,
) = begin
    pomp(
        ;t0 = timezero(g),
        times = times(g),
        userdata = (guide = g,),
        init_state = x.init_state,
        rinit = x.rinitfn,
        rprocess = onestep(x.rprocessfn),
        logdmeasure = x.logdmeasurefn,
        kwargs...,
    )
end

macro guide(x, g, arg)
    x = __module__.eval(x)
    @assert x isa Module && isdefined(x,:knowledge!) "improper call to @guide"
    m = x.eval(arg)
    k = x.knowledge!
    esc(:(PhyloPOMP.guide($g,$m,$k)))
end
