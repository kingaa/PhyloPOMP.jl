const FuncMap = Dict{Tuple{Int,Int},Vector{Symbol}}

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
        $fname = function(node, cols, state, t; kwargs...)
            if node.parlin ∈ cols[$from]
                rate = $rate(t, state; kwargs...)
                enns = $demeocc(state)
                ells = ell(cols)
                p1 = $prop(state, enns, ells; kwargs...)
                p2 = branch_options_prob($(to[1]), $(to[2]), node.present)
                rate*p1*p2
            else
                0
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
    production = pop!(move,:production)
    saturation = pop!(move,:saturation)
    boost = make_boostfn!(production,saturation,move,dict,:fork)
    rate = jump[:ratefn]
    demeocc = dict[:demeoccfn]
    popmove = jump[:popmovefn]
    defn = quote
        $fname = function(node, cols, state, t; kwargs...)
            if node.parlin ∈ cols[$from]
                i, j, p = branch_options_move($(to[1]), $(to[2]), node.present)
                ells = fork!(cols,$from,node.parlin,(i,j),node.chillins)
                rate = $rate(t, state; kwargs...)
                state = $popmove(t, state; kwargs...)
                enns = $demeocc(state)
                boost = $boost(enns, ells)
                log(rate*boost/p), state
            else
                Prob(-Inf), state
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
    production = pop!(move,:production)
    saturation = pop!(move,:saturation)
    boost = make_boostfn!(production,saturation,move,dict,:chop)
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
                ratefuns = [$(v...)]
                movefuns = [$(w...)]
                alpha = map(ratefuns) do f
                    f(node, cols, state, t; kwargs...)
                end
                k,_,p = rcateg(alpha, true)
                if k > 0
                    ll, state = movefuns[k](node, cols, state, t; kwargs...)
                    ll -= log(p)
                else
                    ll = Prob(-Inf)
                end
            end
        else
            quote
                ll, state = $(w[1])(node, cols, state, t; kwargs...)
            end
        end
        (:(nsample == $m && nchild == $n), expr)
    end
    msg = :("unrecognized node signature (nsample = $nsample, nchild = $nchild): model is incompatible with genealogy.")
    body = switchblock(pieces...,:(error($msg)))
    fname = :singular_move!
    dict[:singmovefn] = fname
    defn = quote
        $fname = function(nsample, nchild, node, cols, state, t; kwargs...)
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
    dict[:singrootfn] = fname = :singular_root!
    defn = quote
        $fname = function(node::GuideNode{F,N,D}, cols, state, t; kwargs...) where {F, N, D <: Enum}
            q::Prob = 0
            i, _, p = rcateg(node.present[:,1], true)
            if i > 0
                q = inv(p)
                plant!(cols, D(i), node.chillins[1])
                enns = $demeocc(state)
                ells = ell(cols)
                live = isalive(enns, ells)
            else
                live = false
            end
            live, log(q), state
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
    dict[:singularpartfn] = fname
    defn = quote
        $fname = function(cols, state, guide, n, t; kwargs...)
            ll::Prob = -Inf
            enns = $demeocc(state)
            ells = ell(cols)
            live = isalive(enns, ells)
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
