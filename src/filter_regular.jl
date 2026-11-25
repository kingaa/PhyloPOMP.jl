make_id_move!(move, jump, dict) = begin
    demeocc = dict[:demeoccfn]
    popmove = jump[:popmovefn]
    boost = pop!(move,:boost,one(Prob))
    boost = make_boostfn!(boost,move,dict,:id)
    jump[:regmovefn] = fname = Symbol(:move_,move[:name])
    defn = quote
        const $fname = function(t, guide, node, cols, state; kwargs...)
            local enns, ells, boost
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
    idboost = pop!(move,:idboost,one(Prob))
    swapboost = pop!(move,:swapboost,one(Prob))
    idboost = make_boostfn!(idboost,move,dict,:id)
    swapboost = make_boostfn!(swapboost,move,dict,:swap)
    jump[:regmovefn] = fname = Symbol(:move_,move[:name])
    defn = quote
        const $fname = function(t, guide, node, cols, state; kwargs...)
            local enns, ells, boost, b, p
            enns = $demeocc(state)
            b, p = choose_branch(t, guide, node, enns, cols, $from, $to)
            state = $popmove(t, state; kwargs...)
            enns = $demeocc(state)
            if b == 0
                ells = ell(cols)
                boost = $idboost(enns, ells)
            else
                ells = swap!(cols, $from, $to, b)
                boost = $swapboost(enns, ells)
            end
            log(boost/p), state
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_regmovefn!(dict) = begin
    jumps = filter(
        x -> haskey(x.second,:has_reg),
        get(dict, :jumps, MgpDict())
    )
    dict[:nregjump] = length(jumps)
    names = keys(jumps)
    items = map(enumerate(names)) do (k,name)
        fun = pop!(jumps[name],:regmovefn)
        cond = :((k == $k))
        call = :((ll, state) = $fun(t, guide, node, cols, state; kwargs...))
        (cond, call)
    end
    body = switchblock(
        items...,
        :(error("unreachable condition in `regular_move!` function"))
    )
    fname = :regular_move!
    dict[:regmovefn] = fname
    defn = quote
        const $fname = function(k, t, guide, node, cols, state; kwargs...)
            $body
        end
    end
    push!(dict[:decls],defn)
    nothing
end

make_regeventratefn!(dict) = begin
    jumps = filter(
        x -> haskey(x.second,:has_reg),
        get(dict, :jumps, MgpDict())
    )
    mpieces = map(enumerate(values(jumps))) do (i,jump)
        rate = pop!(jump,:ratefn)
        prop = pop!(jump,:regpropfn)
        quote
            rate += r = $rate(t, state; kwargs...)::Prob
            pi[$i] = $prop(state, enns, ells; kwargs...)::Prob
            alpha[$i] = r*pi[$i]
        end
    end
    jumps = filter(x -> !haskey(x.second,:has_reg), dict[:jumps])
    dpieces = map(values(jumps)) do jump
        rate = jump[:ratefn]
        quote
            rate += $rate(t, state; kwargs...)::Prob
        end
    end
    demeocc = dict[:demeoccfn]
    fname = :reg_event_rates!
    dict[:regeventratesfn] = fname
    defn = quote
        const $fname = function(alpha, pi, t, state, cols; kwargs...)
            local rate, enns, ells
            rate::Prob = 0
            enns = $demeocc(state)
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
    fname = :regular_part!
    defn = quote
        const $fname = function(
            cols::Coloring, state::NamedTuple, guide::Guide,
            n::Integer, t::Time, tf::Time; kwargs...,)
            local node, alpha, pi, step, decay, ll, ll1, k, s
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
                    ll += ll1 - (decay*step + log(pi[k]))
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
