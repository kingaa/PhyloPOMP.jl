"""
    HardSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using a filter guide and "hard" proposals.  That is, color
changes can be proposed on branches at rates that exceed the event rates
in the underlying population process.
"""
module HardSEIR

include("seir_common.jl")

rate_transmission!(
    alpha, pi, t, node, cols, rh, (;S, E, I, R);
    β, pop, _...,
) = begin
    rate = β*S*I/pop
    on = sum_relhaz(rh,node,cols,Infec,Expos)
    pi[1] = @indicator(I > 0, 1-ell(cols,Infec)/(I+1))
    pi[2] = @indicator(I > 0, on/(I+1))
    alpha .= rate.*pi
    rate-sum(alpha)
end

rate_progression!(
    alpha, pi, t, node, cols, rh, (;S, E, I, R);
    σ, _...,
) = begin
    rate = σ*E
    on = sum_relhaz(rh,node,cols,Expos,Infec)
    pi[1] = @indicator(E > 0, 1-ell(cols,Expos)/E)
    pi[2] = @indicator(E > 0, on/E)
    alpha .= rate.*pi
    rate-sum(alpha)
end

rate_recovery!(
    alpha, pi, t, node, cols, rh, (;S, E, I, R);
    γ, _...,
) = begin
    rate = γ*I
    on = ell(cols,Infec)/I
    pi[1] = @indicator(on < 1, 1-on)
    alpha .= rate.*pi
    rate-sum(alpha)
end

rate_waning!(
    alpha, pi, t, node, cols, rh, (;S, E, I, R);
    ω, _...,
) = begin
    rate = ω*R
    pi[1] = 1
    alpha .= rate.*pi
    rate-sum(alpha)
end

rate_sampling(
    t, node, cols, rh, (;S, E, I, R);
    ψ, χ,
    _...,
) = begin
    rate = (ψ+χ)*I
    rate
end

event_rates!(alpha, pi, t, node, cols, rh, state; kwargs...) = begin
    decay = zero(Prob)
    decay += rate_transmission!(@view(alpha[1:2]), @view(pi[1:2]), t, node, cols, rh, state; kwargs...,)
    decay += rate_progression!(@view(alpha[3:4]), @view(pi[3:4]), t, node, cols, rh, state; kwargs...,)
    decay += rate_recovery!(@view(alpha[5]), @view(pi[5]), t, node, cols, rh, state; kwargs...,)
    decay += rate_waning!(@view(alpha[6]), @view(pi[6]), t, node, cols, rh, state; kwargs...,)
    decay += rate_sampling(t, node, cols, rh, state; kwargs...,)
    decay
end

move_transmission_off!(t, (;S, E, I, R), node, cols, rh) = begin
    S -= 1
    E += 1
    boost = @binomratio(E,ell(cols,Expos),1,0)
    log(boost), (;S, E, I, R)
end

move_transmission_on!(t, (;S, E, I, R), node, cols, rh) = begin
    b, p = choose_branch(rh,node,cols,Infec,Expos)
    ellE, ellI = swap!(cols,Infec,Expos,b)
    S -= 1
    E += 1
    boost = @binomratio((E,I),(ellE,ellI),(1,1),(1,0))/p
    log(boost), (;S, E, I, R)
end

move_progression_off!(t, (;S, E, I, R), node, cols, rh) = begin
    E -= 1
    I += 1
    boost = @binomratio(I,ell(cols,Infec),1,0)
    log(boost), (;S, E, I, R)
end

move_progression_on!(t, (;S, E, I, R), node, cols, rh) = begin
    b, p = choose_branch(rh,node,cols,Expos,Infec)
    _, ellI = swap!(cols,Expos,Infec,b)
    E -= 1
    I += 1
    boost = @binomratio(I,ellI,1,1)/p
    log(boost), (;S, E, I, R)
end

move_recovery!(t, (;S, E, I, R), node, cols, rh) = begin
    I -= 1
    R += 1
    boost = 1.0
    log(boost), (;S, E, I, R)
end

move_waning!(t, (;S, E, I, R), node, cols, rh) = begin
    R -= 1
    S += 1
    boost = 1.0
    log(boost), (;S, E, I, R)
end

regular_move!(k, t, node, cols, state, rh) = begin
    if k==1
        move_transmission_off!(t, state, node, cols, rh)
    elseif k==2
        move_transmission_on!(t, state, node, cols, rh)
    elseif k==3
        move_progression_off!(t, state, node, cols, rh)
    elseif k==4
        move_progression_on!(t, state, node, cols, rh)
    elseif k==5
        move_recovery!(t, state, node, cols, rh)
    elseif k==6
        move_waning!(t, state, node, cols, rh)
    else
        @assert false "impossible" # COV_EXCL_LINE
    end
end

regular_part!(
    cols, state, guide, n, t, tf;
    kwargs...,
) = begin
    node = guide[n]
    alpha = similar(Vector{Prob},6)
    pi = similar(Vector{Prob},6)
    rh = relhaz_alloc(guide,n)
    step::Time = zero(Time)
    decay::Prob = zero(Prob)
    ll::Prob = zero(Prob)
    while t < tf
        relhaz!(rh,t,guide,n)
        decay = event_rates!(alpha, pi, t, node, cols, rh, state; kwargs...)
        k, s = rcateg(alpha)
        step = -log(rand())/s
        if k > 0 && t+step < tf
            ll1, state = regular_move!(k, t, node, cols, state, rh)
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
