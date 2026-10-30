"""
    GuidedSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using a filter guide and "semisoft" proposals.  That
is, when population-process events occur, the guide can steer those
events preferentially onto (or away from) particular branches, but the
overall event rate remains equal to that of the underlying population
process.
"""
module GuidedSEIR

include("seir_common.jl")

## Event-rates section.
##
## Each of these functions computes one or more event-rates (α in the
## paper) and proposals (π in the paper) and fills the arrays `alpha`
## and `pi` with these. In addition, it returns any decay term
## associated.  These are only called by the function `event_rates!`
## defined at the bottom of this section.

rate_transmission!(alpha, pi, t, (;S, E, I, R), cols; β, pop, _...,) = begin
    rate = β*S*I/pop
    pi[1] = 1
    alpha[1] = rate*pi[1]
    rate-sum(alpha)
end

rate_progression!(alpha, pi, t, (;S, E, I, R), cols; σ, _...,) = begin
    rate = σ*E
    pi[1] = 1
    alpha[1] = rate*pi[1]
    rate-sum(alpha)
end

rate_recovery!(alpha, pi, t, (;S, E, I, R), cols; γ, _...,) = begin
    rate = γ*I
    ellI = ell(cols,Infec)
    pi[1] = @indicator(I > ellI, 1-ellI/I)
    alpha[1] = rate*pi[1]
    rate-sum(alpha)
end

rate_waning!(alpha, pi, t, (;S, E, I, R), cols; ω, _...,) = begin
    rate = ω*R
    pi[1] = 1
    alpha[1] = rate*pi[1]
    rate-sum(alpha)
end

rate_sampling(t, (;S, E, I, R), cols; ψ, χ, _...,) = begin
    rate = (ψ+χ)*I
    rate
end

event_rates!(alpha, pi, t, state, cols; kwargs...) = begin
    decay = zero(Prob)
    decay += rate_transmission!(@view(alpha[1]), @view(pi[1]), t, state, cols; kwargs...,)
    decay += rate_progression!(@view(alpha[2]), @view(pi[2]), t, state, cols; kwargs...,)
    decay += rate_recovery!(@view(alpha[3]), @view(pi[3]), t, state, cols; kwargs...,)
    decay += rate_waning!(@view(alpha[4]), @view(pi[4]), t, state, cols; kwargs...,)
    decay += rate_sampling(t, state, cols; kwargs...,)
    decay
end

## Regular-part section.
##
## The following functions execute the appropriate actions for each
## kind of event (see Event-rate section above) during the integration
## of the regular part of the filter equation. These are called in
## turn by the `regular_part!` function defined at the bottom of the
## section.

move_transmission!(
    t, guide, node, cols, (;S, E, I, R);
    kwargs...,
) = begin
    b,p = choose_branch(t,guide,node,I,cols,Infec,Expos)
    S -= 1
    E += 1
    if b == 0
        ellE = ell(cols,Expos)
        boost = (1-ellE/E)/p
    else
        ellE, ellI = swap!(cols,Infec,Expos,b)
        boost = (1-ellI/I)/E/p
    end
    log(boost), (;S, E, I, R)
end

move_progression!(
    t, guide, node, cols, (;S, E, I, R);
    kwargs...,
) = begin
    b,p = choose_branch(t,guide,node,E,cols,Expos,Infec)
    E -= 1
    I += 1
    if b == 0
        ellI = ell(cols,Infec)
        boost = (1-ellI/I)/p
    else
        swap!(cols,Expos,Infec,b)
        boost = inv(I*p)
    end
    log(boost), (;S, E, I, R)
end

move_recovery!(
    t, guide, node, cols, (;S, E, I, R);
    kwargs...,
) = begin
    I -= 1
    R += 1
    zero(Prob), (;S, E, I, R)
end

move_waning!(
    t, guide, node, cols, (;S, E, I, R);
    kwargs...,
) = begin
    R -= 1
    S += 1
    zero(Prob), (;S, E, I, R)
end

regular_move!(k, t, guide, node, cols, state) = begin
    if k==1
        move_transmission!(t,guide,node,cols,state)
    elseif k==2
        move_progression!(t,guide,node,cols,state)
    elseif k==3
        move_recovery!(t,guide,node,cols,state)
    elseif k==4
        move_waning!(t,guide,node,cols,state)
    else
        @assert false "unreachable" # COV_EXCL_LINE
    end
end

regular_part!(
    cols, state, guide, n, t, tf;
    kwargs...,
) = begin
    node = guide[n]
    alpha = similar(Vector{Prob},4)
    pi = similar(Vector{Prob},4)
    step::Time = zero(Time)
    decay::Prob = zero(Prob)
    ll::Prob = zero(Prob)
    while t < tf
        decay = event_rates!(alpha, pi, t, state, cols; kwargs...,)
        k, s = rcateg(alpha)
        step = -log(rand())/s
        if k > 0 && t+step < tf
            ll1, state = regular_move!(k, t, guide, n, cols, state)
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
