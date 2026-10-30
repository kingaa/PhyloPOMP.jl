"""
    SoftSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using a filter guide and "soft" proposals.  That is, when
population-process events occur, the guide can steer those events
preferentially onto (or away from) particular branches, but the overall
event rate remains equal to that of the underlying population process.
"""
module SoftSEIR

include("seir_common.jl")

rate_transmission!(
    alpha, pi, t, cols, (;S, E, I, R);
    β, pop, _...,
) = begin
    rate = β*S*I/pop
    on = ell(cols,Infec)
    pi[1] = @indicator(I > 0, (I-on+1)/(I+1))
    pi[2] = @indicator(I > 0, on/(I+1))
    alpha[1] = rate*pi[1]
    alpha[2] = rate*pi[2]
    rate-sum(alpha)
end

rate_progression!(
    alpha, pi, t, cols, (;S, E, I, R);
    σ, _...,
) = begin
    rate = σ*E
    on = ell(cols,Expos)
    pi[1] = @indicator(E > 0, (E-on)/E)
    pi[2] = @indicator(E > 0, on/E)
    alpha[1] = rate*pi[1]
    alpha[2] = rate*pi[2]
    rate-sum(alpha)
end

rate_recovery!(
    alpha, pi, t, cols, (;S, E, I, R);
    γ, _...,
) = begin
    rate = γ*I
    ellI = ell(cols,Infec)
    pi[1] = @indicator(I > ellI, 1-ellI/I)
    alpha[1] = rate*pi[1]
    rate-sum(alpha)
end

rate_waning!(
    alpha, pi, t, cols, (;S, E, I, R);
    ω, _...,
) = begin
    rate = ω*R
    pi[1] = 1
    alpha[1] = rate*pi[1]
    rate-sum(alpha)
end

rate_sampling(
    t, cols, (;S, E, I, R);
    ψ, χ,
    _...,
) = begin
    rate = (ψ+χ)*I
    rate
end

event_rates!(alpha, pi, t, cols, state; kwargs...) = begin
    decay = zero(Prob)
    decay += rate_transmission!(@view(alpha[1:2]), @view(pi[1:2]), t, cols, state; kwargs...,)
    decay += rate_progression!(@view(alpha[3:4]), @view(pi[3:4]), t, cols, state; kwargs...,)
    decay += rate_recovery!(@view(alpha[5]), @view(pi[5]), t, cols, state; kwargs...,)
    decay += rate_waning!(@view(alpha[6]), @view(pi[6]), t, cols, state; kwargs...,)
    decay += rate_sampling(t, cols, state; kwargs...,)
    decay
end

move_transmission_off!(t, (;S, E, I, R), guide, n, cols, rh) = begin
    S -= 1
    E += 1
    boost = @binomratio(E,ell(cols,Expos),1,0)
    log(boost), (;S, E, I, R)
end

move_transmission_on!(t, (;S, E, I, R), guide, n, cols, rh) = begin
    relhaz!(rh,t,guide,n)
    b, p = choose_branch(rh,guide[n],cols,Infec,Expos)
    _, ellI = swap!(cols,Infec,Expos,b)
    S -= 1
    E += 1
    boost = @binomratio(E,1,1,1)*@binomratio(I,ellI,1,0)/p
    log(boost), (;S, E, I, R)
end

move_progression_off!(t, (;S, E, I, R), guide, n, cols, rh) = begin
    E -= 1
    I += 1
    boost = @binomratio(I,ell(cols,Infec),1,0)
    log(boost), (;S, E, I, R)
end

move_progression_on!(t, (;S, E, I, R), guide, n, cols, rh) = begin
    relhaz!(rh,t,guide,n)
    b, p = choose_branch(rh,guide[n],cols,Expos,Infec)
    ellE, ellI = swap!(cols,Expos,Infec,b)
    E -= 1
    I += 1
    boost = @binomratio(I,ellI,1,1)/p
    log(boost), (;S, E, I, R)
end

move_recovery!(t, (;S, E, I, R), guide, n, cols, rh) = begin
    I -= 1
    R += 1
    boost = 1.0
    log(boost), (;S, E, I, R)
end

move_waning!(t, (;S, E, I, R), guide, n, cols, rh) = begin
    R -= 1
    S += 1
    boost = 1.0
    log(boost), (;S, E, I, R)
end

regular_move!(k, t, guide, n, cols, state, rh) = begin
    if k==1
        move_transmission_off!(t, state, guide, n, cols, rh)
    elseif k==2
        move_transmission_on!(t, state, guide, n, cols, rh)
    elseif k==3
        move_progression_off!(t, state, guide, n, cols, rh)
    elseif k==4
        move_progression_on!(t, state, guide, n, cols, rh)
    elseif k==5
        move_recovery!(t, state, guide, n, cols, rh)
    elseif k==6
        move_waning!(t, state, guide, n, cols, rh)
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
        decay = event_rates!(alpha, pi, t, cols, state; kwargs...)
        k, s = rcateg(alpha)
        step = -log(rand())/s
        if k > 0 && t+step < tf
            ll1, state = regular_move!(k, t, guide, n, cols, state, rh)
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
