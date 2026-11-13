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
                boost = $offboost(enns, ells)
            else
                ells = swap!(cols, $from, $to, b)
                boost = $onboost(enns, ells)
            end
            log(boost/p), state
        end
    end
    push!(dict[:decls],defn)
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
        $fname = function(k, t, guide, node, cols, state; kwargs...)
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
            pi[$i] = $prop(state, enns, ells; kwargs...)
            alpha[$i] = r*pi[$i]
        end
    end
    jumps = filter(x -> !haskey(x.second,:has_reg), dict[:jumps])
    dpieces = map(values(jumps)) do jump
        rate = jump[:ratefn]
        quote
            rate += $rate(t, state; kwargs...)
        end
    end
    demeocc = dict[:demeoccfn]
    fname = :reg_event_rates!
    dict[:regeventratesfn] = fname
    defn = quote
        $fname = function(alpha, pi, t, state, cols; kwargs...)
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
