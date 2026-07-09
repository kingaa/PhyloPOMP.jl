"""
    HardSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using a filter guide and "hard" proposals.  That is, color
changes can be proposed on branches at rates that exceed the event rates
in the underlying population process.
"""
module HardSEIR

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Expos Infec
using .Demes: Expos, Infec, DemeSet

include("seir_trees.jl")
include("seir_funs.jl")

transmission!(
    alpha, pi, srh;
    β, S, I, pop, ellI,
    _...,
) = begin
    rate = β*S*I/pop
    pi[1] = @indicator(I > 0, 1-ellI/I)
    pi[2] = @indicator(I > 0, srh/I)
    alpha[1] = rate*pi[1]
    alpha[2] = rate*pi[2]
    rate-sum(alpha)
end

progression!(
    alpha, pi, srh;
    σ, E, ellE,
    _...,
) = begin
    rate = σ*E
    pi[1] = @indicator(E > 0, 1-ellE/E)
    pi[2] = @indicator(E > 0, srh/E)
    alpha[1] = rate*pi[1]
    alpha[2] = rate*pi[2]
    rate-sum(alpha)
end

recovery!(
    alpha, pi;
    γ, I, ellI,
    _...,
) = begin
    rate = γ*I
    pi[1] = @indicator(I > ellI, 1-ellI/I)
    alpha[1] = rate*pi[1]
    rate-sum(alpha)
end

waning!(
    alpha, pi;
    ω, R,
    _...,
) = begin
    rate = ω*R
    pi[1] = 1
    alpha[1] = rate
    zero(Float64)
end

sampling(
    ;ψ, χ, I,
    _...,
) = begin
    (ψ+χ)*I
end

event_rates!(
    alpha, pi, rh,
    t, guide, node, cols;
    kwargs...,
) = begin
    relhaz!(rh,t,guide,node)
    decay = zero(Prob)
    decay += transmission!(
        @view(alpha[1:2]),@view(pi[1:2]),
        sum_relhaz(rh,guide[node],cols,Infec,Expos);
        kwargs...,
    )
    decay += progression!(
        @view(alpha[3:4]),@view(pi[3:4]),
        sum_relhaz(rh,guide[node],cols,Expos,Infec);
        kwargs...,
    )
    decay += recovery!(
        @view(alpha[5]),@view(pi[5]);
        kwargs...,
    )
    decay += waning!(
        @view(alpha[6]),@view(pi[6]);
        kwargs...,
    )
    decay += sampling(;kwargs...,)
    decay
end

regular_part!(
    cols, guide, node, ll,
    t, tf, S, E, I, R;
    kwargs...,
) = begin
    n = guide[node]
    alpha = similar(Vector{Prob},6)
    pi = similar(Vector{Prob},6)
    rh = relhaz_alloc(guide,node)
    step::Time = zero(Time)
    decay::Prob = zero(Prob)
    ellE, ellI = ell(cols)
    while t < tf
        decay = event_rates!(
            alpha, pi, rh,
            t, guide, node, cols;
            ellE, ellI,
            S, E, I, R,
            kwargs...,
        )
        k, s = rcateg(alpha)
        step = -log(rand())/s
        if t+step < tf
            ll -= decay*step + log(pi[k])
            if k==1
                S -= 1
                E += 1
                ll += log(1-ellE/E)
            elseif k==2
                rates = relhaz(rh,n,cols,Infec,Expos)
                b,_,p = rcateg(rates,cols[Infec],true)
                ll -= log(p)
                S -= 1
                E += 1
                ellE, ellI = swap!(cols,Infec,Expos,b)
                ll += log(1-ellI/I)-log(E)
            elseif k==3
                E -= 1
                I += 1
                ll += log(1-ellI/I)
            elseif k==4
                rates = relhaz(rh,n,cols,Expos,Infec)
                b,_,p = rcateg(rates,cols[Expos],true)
                ll -= log(p)
                E -= 1
                I += 1
                ellE, ellI = swap!(cols,Expos,Infec,b)
                ll -= log(I)
            elseif k==5
                I -= 1
                R += 1
            elseif k==6
                R -= 1
                S += 1
            else
                @assert false "impossible error!" # COV_EXCL_LINE
            end
            t += step
        else
            step = tf - t
            ll -= decay*step
            break
        end
    end
    @assert I ≥ ellI && E ≥ ellE
    ll, S, E, I, R
end

end
