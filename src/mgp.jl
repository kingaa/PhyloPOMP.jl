include("extract_symbols.jl")
include("parse_mgp.jl")
include("switchblock.jl")

struct GuidedMGP
    demes::Module
    decls::Vector{Expr}
    rinit::Function
    rprocess::Function
end

const FuncMap = Dict{Tuple{Int,Int},Vector{Symbol}}

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
        :decls => Expr[]
    )
    parse_mgp_spec!(dict, args...)
    repair_mgp_spec!(dict)
    process_demes!(dict)
    process_livecond!(dict)
    process_jumps!(dict)
    rinit = make_rinitfn(dict)
    rprocess = make_rprocessfn!(dict)
    quote
        @eval(@demes($(dict[:deme_module]), $(dict[:demes]...)))

        PhyloPOMP.GuidedMGP(
            $(dict[:deme_module]),
            $(dict[:decls]),
            $rinit,
            $rprocess
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

process_livecond!(dict) = begin
    enns = dict[:enns]
    ells = dict[:ells]
    live_condition = pop!(dict,:live_condition)
    fname = :isalive
    dict[:livecondfn] = fname
    defn = quote
        isalive = function($enns, $ells)
            $live_condition
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
        make_fork_move!(move, jump, dict)
    elseif type == :chop
        jump[:has_sing] = :true
        make_chop_move!(move, jump, dict)
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

make_id_move!(move, jump, dict) = begin
    demeocc = dict[:demeoccfn]
    popmove = jump[:popmovefn]
    production = pop!(move,:production)
    saturation = pop!(move,:saturation)
    boost = make_boostfn!(production,saturation,move,dict,:id)
    jump[:regmovefn] = fname = Symbol(:move_,move[:name])
    defn = quote
        $fname = function(t, guide, node, cols, state; kwargs...)
            state = $popmove(t, state; kwargs...)
            enns = $demeocc(state)
            ells = ell(cols)
            boost = $boost(enns, ells)
            log(boost), state
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_swap_move!(move, jump, dict) = begin
    demeocc = dict[:demeoccfn]
    popmove = jump[:popmovefn]
    from = pop!(move,:from)
    to = pop!(move,:to)
    production = pop!(move,:production)
    saturation = pop!(move,:saturation)
    onboost = make_boostfn!(production,saturation,move,dict,:swap)
    offboost = make_boostfn!(production,dict[:nullproduc],move,dict,:id)
    jump[:regmovefn] = fname = Symbol(:move_,move[:name])
    defn = quote
        $fname = function(t, guide, node, cols, state; kwargs...)
            enns = $demeocc(state)
            b, p = choose_branch(t, guide, node, enns, cols, $from, $to)
            state = $popmove(t, state; kwargs...)
            enns = $demeocc(state)
            if b == 0
                ells = ell(cols)
                boost = $offboost(n, ells)
            else
                ells = swap!(cols, $from, $to, b)
                boost = $onboost(n, ells)
            end
            log(boost/p), state
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_fork_move!(move, jump, dict) = begin
    signature = (0,2)
    fname = Symbol(:move_,move[:name])
    from = pop!(move,:from)
    to = pop!(move,:to).args
    production = pop!(move,:production)
    saturation = pop!(move,:saturation)
    boost = make_boostfn!(production,saturation,move,dict,:fork)
    rate = jump[:ratefn]
    demeocc = dict[:demeoccfn]
    popmove = jump[:popmovefn]
    defn = quote
        $fname = function(node, cols, state, t; kwargs...)
            if node.parlin ∈ cols[$from]
                k, _, p = rcateg(branch_options(node.present), true) # FIXME: care is needed in defn of 'branch_options'
                rate = $rate(t, state; kwargs...)
                if k==1
                    ells = fork!(cols,$from,node.parlin,($(to[1]),$(to[2])),node.chillins)
                else
                    ells = fork!(cols,$from,node.parlin,($(to[2]),$(to[1])),node.chillins)
                end
                state = $popmove(t, state; kwargs...)
                enns = $demeocc(state)
                boost = $boost(enns,ells)
                log(rate*boost/p), state
            else
                Prob(-Inf), state
            end
        end
    end
    push!(dict[:decls],defn)
    if !haskey(dict[:movemap],signature)
        dict[:movemap][signature] = Symbol[]
    end
    push!(dict[:movemap][signature], fname)
    fname = Symbol(:prop_,move[:name])
    prop = pop!(move,:propfn)
    defn = quote
        $fname = function(node, cols, state, t; kwargs...)
            if node.parlin ∈ cols[$from]
                rate = $rate(t, state; kwargs...)
                enns = $demeocc(state)
                ells = ell(cols)
                prop = $prop(state, enns, ells; kwargs...)
                rate*prop
            else
                0
            end
        end
    end
    push!(dict[:decls],defn)
    if !haskey(dict[:eventmap],signature)
        dict[:eventmap][signature] = Symbol[]
    end
    push!(dict[:eventmap][signature], fname)
    nothing
end

make_chop_move!(move, jump, dict) = begin
    signature = eval(move[:signature])
    fname = Symbol(:move_,move[:name])
    from = pop!(move,:from)
    to = if haskey(move,:to)
        pop!(move,:to)
    else
        nothing
    end
    production = pop!(move,:production)
    saturation = pop!(move,:saturation)
    boost = make_boostfn!(production,saturation,move,dict,:chop)
    rate = jump[:ratefn]
    demeocc = dict[:demeoccfn]
    popmove = jump[:popmovefn]
    op = if isnothing(to)
        :(chop!(cols, $from, node.parlin))
    else
        :(chop!(cols, $from, node.parlin, $to, node.chillins[1]))
    end
    defn = quote
        $fname = function(node, cols, state, t; kwargs...)
            if node.parlin ∈ cols[$from]
                rate = $rate(t, state; kwargs...)
                state = $popmove(t, state; kwargs...)
                ells = $op
                enns = $demeocc(state)
                boost = $boost(enns,ells)
                log(rate*boost), state
            else
                Prob(-Inf), state
            end
        end
    end
    push!(dict[:decls],defn)
    if !haskey(dict[:movemap],signature)
        dict[:movemap][signature] = Symbol[]
    end
    push!(dict[:movemap][signature], fname)
    fname = Symbol(:prop_,move[:name])
    prop = pop!(move,:propfn)
    defn = quote
        $fname = function(node, cols, state, t; kwargs...)
            if node.parlin ∈ cols[$from]
                rate = $rate(t, state; kwargs...)
                enns = $demeocc(state)
                ells = ell(cols)
                prop = $prop(state, enns, ells; kwargs...)
                rate*prop
            else
                0
            end
        end
    end
    push!(dict[:decls],defn)
    if !haskey(dict[:eventmap],signature)
        dict[:eventmap][signature] = Symbol[]
    end
    push!(dict[:eventmap][signature], fname)
    nothing
end

make_regmovefn!(dict) = begin
    @assert haskey(dict,:jumps) "no jumps specified: null process"
    jumps = filter(x -> haskey(x.second,:has_reg), dict[:jumps])
    dict[:nregjump] = length(jumps)
    names = keys(jumps)
    items = map(enumerate(names)) do (k,name)
        fun = pop!(jumps[name],:regmovefn)
        cond = :((k == $k))
        call = :((ll, state) = $fun(t, guide, node, cols, state; kwargs...))
        (cond, call)
    end
    body = switchblock(items..., :(error("unreachable condition in `regular_move!` function")))
    fname = :regular_move!
    dict[:regmovefn] = fname
    defn = quote
        $fname = function(k, t, guide, node, cols, state)
            $body
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_regeventratefn!(dict) = begin
    @assert haskey(dict,:jumps) "no jumps specified: null process"
    jumps = filter(x -> haskey(x.second,:has_reg), dict[:jumps])
    mpieces = map(enumerate(values(jumps))) do (i,jump)
        rate = pop!(jump,:ratefn)
        prop = pop!(jump,:regpropfn)
        quote
            rate += r = $rate(t, state; kwargs...)
            pi = $prop(state, ells; kwargs...)
            alpha[$i] = r*pi
        end
    end
    jumps = filter(x -> !haskey(x.second,:has_reg), dict[:jumps])
    dpieces = map(values(jumps)) do jump
        rate = jump[:ratefn]
        quote
            rate += $rate(t, state; kwargs...)
        end
    end
    fname = :reg_event_rates!
    dict[:regeventratesfn] = fname
    defn = quote
        $fname = function(alpha, pi, t, state, cols; kwargs...)
            rate::Prob = 0
            ells = ell(cols)
            $(mpieces...)
            $(dpieces...)
            rate-sum(alpha)
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_regular_part!(dict) = begin
    make_regeventratefn!(dict)
    make_regmovefn!(dict)
    nj = pop!(dict,:nregjump)
    eventratesfn = pop!(dict,:regeventratesfn)
    regmovefn = pop!(dict,:regmovefn)
    dict[:regularpartfn] = fname = :regular_part!
    defn = quote
        $fname = function(cols, state, guide, n, t, tf; kwargs...,)
            node = guide[n]
            alpha = similar(Vector{Prob},$nj)
            pi = similar(Vector{Prob},$nj)
            step::Time = 0
            decay::Prob = 0
            ll::Prob = 0
            while t < tf
                decay = $eventratesfn(alpha, pi, t, state, cols; kwargs...,)
                k, s = rcateg(alpha)
                step = -log(rand())/s
                if k > 0 && t+step < tf
                    ll1, state = $regmovefn(k, t, guide, n, cols, state)
                    ll += ll1 - decay*step - log(pi[k])
                    t += step
                else
                    step = tf - t
                    ll -= decay*step
                    break
                end
            end
            ll, state
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_singmovefn!(dict) = begin
    @assert haskey(dict,:eventmap) && haskey(dict,:movemap)
    eventmap = pop!(dict,:eventmap)
    movemap = pop!(dict,:movemap)
    sigs = collect(keys(eventmap))
    pieces = map(sigs) do (m,n)
        v = eventmap[(m,n)]
        w = movemap[(m,n)]
        nv = length(v)
        expr = if nv > 1
            quote
                ratefuns = [$(v...)]
                movefuns = [$(w...)]
                alpha = map(ratefuns) do f
                    f(node, cols, state, t; kwargs...)
                end
                k,_,p = rcateg(alpha, true)
                @assert k > 0
                ll, state = movefuns[k](node, cols, state, t; kwargs...)
                ll -= log(p)
            end
        else
            quote
                ll, state = $(w[1])(node, cols, state, t; kwargs...)
            end
        end
        (:(nsample == $m && nchild == $n), expr)
    end
    msg = "unrecognized node signature ($(Expr(:$,:nsample)), $(Expr(:$,:nchild))): model is incompatible with genealogy."
    body = switchblock(pieces...,:(error($msg)))
    fname = :singular_move!
    dict[:singmovefn] = fname
    defn = quote
        $fname = function(nsample, nchild, node, t, state, cols; kwargs...)
            $body
            live = ll > Prob(-Inf)
            live, ll, state
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_singular_rootfn!(dict) = begin
    demeocc = dict[:demeoccfn]
    livecond = dict[:livecondfn]
    dict[:singrootfn] = fname = :singular_root!
    defn = quote
        $fname = function(node::GuideNode{F,N,D}, cols, state, t; kwargs...) where {F, N, D <: Enum}
            ll::Prob = -Inf
            enns = $demeocc(state)
            ells = ell(cols)
            live = $livecond(enns, ells)
            if live
                ## FIXME: should this be proportional to (n-ℓ)?
                i, _, p = rcateg(node.present[:,1].*(enns .- ells), true)
                if i > 0
                    ll = -log(p)
                    plant!(cols, D(i), node.chillins[1])
                else
                    live = false
                end
            end
            live, ll, state
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_singular_part!(dict) = begin
    make_singmovefn!(dict)
    make_singular_rootfn!(dict)
    livecond = dict[:livecondfn]
    fname = :singular_part!
    rootfn = pop!(dict,:singrootfn)
    singmovefn = pop!(dict,:singmovefn)
    dict[:singularpartfn] = fname
    defn = quote
        $fname = function(cols, state, guide, n, t; kwargs...)
            ll::Prob = -Inf
            live = $livecond(state, ell(cols))
            if live
                node = guide[n]
                nchild = length(node.chillins)
                if node.type == PhyloPOMP.Root
                    @assert nchild==1 "a root should have exactly one child"
                    live, ll, state = $rootfn(node, cols, state, t; kwargs...)
                elseif node.type == PhyloPOMP.Node
                    live, ll, state = $singmovefn(0, nchild, node, cols, state, t; kwargs...)
                elseif node.type == PhyloPOMP.Sample
                    live, ll, state = $singmovefn(1, nchild, node, cols, state, t; kwargs...)
                end
            end
            live, ll, state
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_rinitfn(dict) = begin
    @assert haskey(dict,:rinit) "rinit not specified"
    fname = :rinitfn
    rinit = pop!(dict,:rinit)
    defn = :($fname = $rinit)
    quote
        let
            $defn
            function (;kwargs...)
                (
                    node = one(Name),
                    ll = zero(Prob),
                    live = true,
                    cols = Coloring(Demes),
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
    demedecls = pop!(dict,:demedecls)
    quote
        let
            $demedecls
            $(dict[:decls]...)
            $fname
        end
    end
end
