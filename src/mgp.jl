include("extract_symbols.jl")
include("parse_mgp.jl")
include("switchblock.jl")
include("filter_singular.jl")
include("filter_regular.jl")

using PartiallyObservedMarkovProcesses: pomp

struct GuidedMGP
    demes::Module
    decls::Vector{Expr}
    init_state::NamedTuple
    rinit::Function
    rprocess::Function
    logdmeasure::Function
    knowledge::Function
    fsmp::FSMarkovProc
end

macro mgp_soft_filter(args...)
    esc(soft_filter_def(args...))
end

soft_filter_def(args::Expr...) = begin
    dict = MgpDict(
        :deme_module => Symbol(:Demes_,gensym()),
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
        :types => quote
            Name = PhyloPOMP.Name
            Time = PhyloPOMP.Time
            Prob = PhyloPOMP.Prob
            Root = PhyloPOMP.Root
            Node = PhyloPOMP.Node
            Sample = PhyloPOMP.Sample
        end,
        :decls => Expr[]
    )
    parse_mgp_spec!(dict, args...)
    repair_mgp_spec!(dict)
    process_demes!(dict)
    process_jumps!(dict)
    rinit = make_rinitfn(dict)
    init_state = make_init_state(dict)
    rprocess = make_rprocessfn!(dict)
    logdmeas = make_logdmeasure!(dict)
    knowledge = make_knowledge!(dict)
    fsmp = make_fsmp!(dict)
    decls = Expr[
        :(knowledge! = $knowledge),
        :(fsmp = $fsmp),
        :(logdmeasure = $logdmeas),
        :(init_state = $init_state),
        :(rinit = $rinit),
        :(rprocess = $rprocess),
    ]
    quote
        @eval @demes($(dict[:deme_module]), $(dict[:demes]...))

        PhyloPOMP.GuidedMGP(
            $(dict[:deme_module]),
            $decls,
            $init_state,
            $rinit,
            $rprocess,
            $logdmeas,
            $knowledge,
            $fsmp,
        )
    end
end

showdefs(obj::GuidedMGP; _...) = begin
    decls = Base.remove_linenums!.(copy(obj.decls))
    strings = map(decls) do decl
        if decl.head == :block
            string(decl.args...)
        else
            string(decl)
        end
    end
    join(strings,"\n\n")
end

Base.show(
    io::IO,
    obj::GuidedMGP,
    kwargs...,
) = begin
    print(io, showdefs(obj; kwargs...))
end

process_demes!(dict) = begin
    @assert haskey(dict,:demes) "demes not specified"
    @assert dict[:demes].head == :tuple "demes should be specified as a tuple"
    dict[:demes] = dict[:demes].args
    nd = length(dict[:demes])
    dict[:nullproduc] = (collect(fill(0,nd))...,)
    demedecls = map(dict[:demes]) do d
        :($d = $(dict[:deme_module]).$d)
    end
    dict[:demedecls] = quote
        $(demedecls...)
    end
    state = dict[:state]
    dict[:demeoccfn] = fname = :deme_occupancy
    @assert haskey(dict,:deme_occupancy) "deme_occupancy not specified"
    demeocc = pop!(dict,:deme_occupancy)
    defn = quote
        $fname = function($state,)
            $demeocc
        end
    end
    push!(dict[:decls],defn)
    nothing
end

process_jumps!(dict) = begin
    @assert haskey(dict, :jumps) "no jumps specified: null process"
    jumps = dict[:jumps]
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
        $fname = function(t, $state; $(vars...), _...)
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
        $fname = function(t, $state; $(vars...), _...)
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
        $fname = function($state, $enns, $ells; $(vars...), _...)
            $prop
        end
    end
    push!(dict[:decls],defn)
    nothing
end

## construct the boost function for a move
make_boostfn!(production, saturation, move, dict, name) = begin
    fname = Symbol(:boost_,move[:name],:_,name)
    defn = :($fname = PhyloPOMP.binomratio($production,$saturation))
    push!(dict[:decls],defn)
    fname
end

make_init_state(dict) = begin
    @assert haskey(dict,:init_state) "init_state not specified"
    init_state = eval(pop!(dict,:init_state))
    types = dict[:types]
    demes = dict[:deme_module]
    quote
        let
            $types
            (
                node=one(Name),
                ll=zero(Prob),
                live=true,
                cols=Coloring($demes),
                state=$init_state,
            )
        end
    end
end

make_rinitfn(dict) = begin
    @assert haskey(dict,:rinit) "rinit not specified"
    fname = :rinitfn
    rinit = pop!(dict,:rinit)
    types = dict[:types]
    demes = dict[:deme_module]
    defn = :($fname = $rinit)
    quote
        let
            $types
            $defn
            function (;kwargs...)
                (
                    node = one(Name),
                    ll = zero(Prob),
                    live = true,
                    cols = Coloring($demes),
                    state = $fname(;kwargs...),
                )
            end
        end
    end
end

make_rprocessfn!(dict) = begin
    make_regular_part!(dict)
    make_singular_part!(dict)
    fname = :rprocessfn
    defn = quote
        $fname = function (;t, dt, guide, node, ll, live, cols, state, kwargs...,)
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
        end
    end
    push!(dict[:decls],defn)
    types = dict[:types]
    demedecls = dict[:demedecls]
    quote
        let
            $types
            $demedecls
            $(dict[:decls]...)
            $fname
        end
    end
end

make_logdmeasure!(dict) = begin
    @assert haskey(dict,:logdmeasure) "logdmeasure not specified"
    fname = :logdmeas
    f = pop!(dict,:logdmeasure)
    types = dict[:types]
    demes = dict[:deme_module]
    defn = :($fname = $f)
    quote
        let
            $types
            $defn
            function (; node, ll, live, cols, state, kwargs...)
                if live
                    ll + $fname(;state..., kwargs...)
                else
                    Prob(-Inf)
                end
            end
        end
    end
end

make_knowledge!(dict) = begin
    @assert haskey(dict,:knowledge) "knowledge not specified"
    know = pop!(dict,:knowledge)
    types = dict[:types]
    demes = dict[:deme_module]
    demedecls = dict[:demedecls]
    quote
        let
            $types
            $demedecls
            $know
        end
    end
end

make_fsmp!(dict) = begin
    @assert haskey(dict,:guide) "guide not specified"
    g = pop!(dict,:guide)
    types = dict[:types]
    demes = dict[:deme_module]
    demedecls = dict[:demedecls]
    quote
        let
            $types
            $demedecls
            $g
        end
    end
end

filter_pomp(
    x::GuidedMGP,
    g::Genealogy;
    kwargs...,
) = begin
    guidegen = guide(g, x.fsmp, x.knowledge)
    pomp(
        ;t0 = timezero(guidegen),
        times = times(guidegen),
        userdata = (guide = guidegen,),
        init_state = x.init_state,
        rinit = x.rinit,
        rprocess = onestep(x.rprocess),
        logdmeasure = x.logdmeasure,
        kwargs...,
    )
end
