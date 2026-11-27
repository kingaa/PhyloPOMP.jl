isalive(enns, ells) = begin
    all(enns .≥ ells)
end

make_fork_propnmove!(move, jump, dict) = begin
    move[:to] = move[:to].args
    make_fork_prop!(move, jump, dict)
    make_fork_move!(move, jump, dict)
    delete!(move,:from)
    delete!(move,:to)
    nothing
end

make_fork_prop!(move, jump, dict) = begin
    signature = (0,2)
    from = move[:from]
    to = move[:to]
    rate = jump[:ratefn]
    demeocc = dict[:demeoccfn]
    prop = pop!(move,:propfn)
    fname = Symbol(:prop_,move[:name])
    defn = quote
        const $fname = function(
            node::GenealNode, cols::Coloring,
            state::NamedTuple, t::Time; kwargs...)
            if node.parlin ∈ cols[$from]
                local rate = $rate(t, state; kwargs...)::Prob
                local enns = $demeocc(state)
                local ells = ell(cols)
                local p1 = $prop(state, enns, ells; kwargs...)::Prob
                local p2 = branch_options_prop($(to[1]), $(to[2]), node.present)::Prob
                rate*p1*p2
            else
                zero(Prob)
            end
        end
    end
    push!(dict[:decls],defn)
    if !haskey(dict[:eventmap], signature)
        dict[:eventmap][signature] = Symbol[]
    end
    push!(dict[:eventmap][signature], fname)
    nothing
end

make_fork_move!(move, jump, dict) = begin
    signature = (0,2)
    fname = Symbol(:move_,move[:name])
    from = move[:from]
    to = move[:to]
    boost = pop!(move,:boost,one(Prob))
    boost = make_boostfn!(boost,move,dict,:fork)
    rate = jump[:ratefn]
    demeocc = dict[:demeoccfn]
    popmove = jump[:popmovefn]
    defn = quote
        const $fname = function(
            node::GuideNode, cols::Coloring, state::S, t::Time;
            kwargs...) where {S <: NamedTuple}
            if node.parlin ∈ cols[$from]
                local rate = $rate(t, state; kwargs...)::Prob
                local s = $popmove(t, state; kwargs...)::S
                local enns = $demeocc(s)
                local i, j, p = branch_options_move($(to[1]), $(to[2]), node.present)
                local ells = fork!(cols,$from,node.parlin,(i,j),node.chillins)
                local boost = $boost(enns, ells)::Prob
                local ll::Prob = log(rate*boost/p)
                return (ll, s)
            else
                return (Prob(-Inf), state)
            end
        end
    end
    push!(dict[:decls],defn)
    if !haskey(dict[:movemap], signature)
        dict[:movemap][signature] = Symbol[]
    end
    push!(dict[:movemap][signature], fname)
    nothing
end

make_chop_propnmove!(move, jump, dict) = begin
    move[:signature] = eval(move[:signature])
    make_chop_prop!(move, jump, dict)
    make_chop_move!(move, jump, dict)
    delete!(move, :signature)
    nothing
end

make_chop_prop!(move, jump, dict) = begin
    signature = move[:signature]
    from = move[:from]
    prop = pop!(move,:propfn)
    rate = jump[:ratefn]
    demeocc = dict[:demeoccfn]
    fname = Symbol(:prop_,move[:name])
    defn = quote
        const $fname = function(
            node::GuideNode, cols::Coloring,
            state::NamedTuple, t::Time; kwargs...)
            if node.parlin ∈ cols[$from]
                local rate = $rate(t, state; kwargs...)::Prob
                local enns = $demeocc(state)
                local ells = ell(cols)
                local prop = $prop(state, enns, ells; kwargs...)::Prob
                rate*prop
            else
                zero(Prob)
            end
        end
    end
    push!(dict[:decls], defn)
    if !haskey(dict[:eventmap], signature)
        dict[:eventmap][signature] = Symbol[]
    end
    push!(dict[:eventmap][signature], fname)
    nothing
end

make_chop_move!(move, jump, dict) = begin
    signature = move[:signature]
    from = pop!(move, :from)
    to = pop!(move, :to, nothing)
    boost = pop!(move,:boost,one(Prob))
    boost = make_boostfn!(boost,move,dict,:chop)
    rate = jump[:ratefn]
    popmove = jump[:popmovefn]
    demeocc = dict[:demeoccfn]
    op = if isnothing(to)
        :(chop!(cols, $from, node.parlin))
    else
        :(chop!(cols, $from, node.parlin, $to, node.chillins[1]))
    end
    fname = Symbol(:move_,move[:name])
    defn = quote
        const $fname = function(
            node::GuideNode, cols::Coloring,
            state::S, t::Time; kwargs...) where {S <: NamedTuple}
            if node.parlin ∈ cols[$from]
                local rate = $rate(t, state; kwargs...)::Prob
                local s = $popmove(t, state; kwargs...)::S
                local enns = $demeocc(s)
                local ells = $op
                local boost = $boost(enns,ells)::Prob
                return (log(rate*boost), s)
            else
                return (Prob(-Inf), state)
            end
        end
    end
    push!(dict[:decls],defn)
    if !haskey(dict[:movemap], signature)
        dict[:movemap][signature] = Symbol[]
    end
    push!(dict[:movemap][signature], fname)
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
                local movefuns = ($(w...),)
                local alpha = [f(node, cols, state, t; kwargs...)::Prob
                               for f ∈ ($(v...),)]
                local k,_,p = rcateg(alpha, true)
                if k > 0
                    ll, s = movefuns[k](node, cols, state, t; kwargs...)
                    ll -= log(p)
                else
                    ll = Prob(-Inf)
                    s = state
                end
            end
        else
            quote
                ll, s = $(w[1])(node, cols, state, t; kwargs...)
            end
        end
        (:(nsample == $m && nchild == $n), expr)
    end
    msg = :("unrecognized node signature (nsample = $nsample, nchild = $nchild): model is incompatible with genealogy.")
    body = switchblock(pieces...,:(error($msg)))
    fname = :singular_move!
    dict[:singmovefn] = fname
    defn = quote
        const $fname = function(
            nsample::Integer, nchild::Integer, node::GuideNode,
            cols::Coloring, state::S, t::Time; kwargs...) where {S <: NamedTuple}
            local live::Bool
            local ll::Prob
            local s::S
            $body
            live = ll > Prob(-Inf)
            live, ll, s
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_singular_rootfn!(dict) = begin
    demeocc = dict[:demeoccfn]
    dict[:singrootfn] = fname = :singular_root!
    defn = quote
        const $fname = function(node::GuideNode{F,N,D}, cols, state, t; kwargs...) where {F, N, D <: Enum}
            local i, _, p = rcateg(node.present[:,1], true)
            if i > 0
                local q = inv(p)
                plant!(cols, D(i), node.chillins[1])
                local enns = $demeocc(state)
                local ells = ell(cols)
                local live = isalive(enns, ells)
                return live, log(q), state
            else
                return false, Prob(-Inf), state
            end
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_singular_part!(dict) = begin
    make_singmovefn!(dict)
    make_singular_rootfn!(dict)
    demeocc = dict[:demeoccfn]
    fname = :singular_part!
    rootfn = pop!(dict,:singrootfn)
    singmovefn = pop!(dict,:singmovefn)
    defn = quote
        const $fname = function(
            cols::Coloring, state::NamedTuple, guide::Guide,
            n::Integer, t::Time; kwargs...)
            local ll::Prob = -Inf
            local enns = $demeocc(state)
            local ells = ell(cols)
            local live = isalive(enns, ells)
            if live
                local node = guide[n]
                local nchild = length(node.chillins)
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
