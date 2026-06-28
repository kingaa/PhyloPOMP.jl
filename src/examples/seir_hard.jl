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

event_rates!(
    alpha, preboost, t, guide, node,
    cols, ellE, ellI,
    S, E, I, R;
    β, σ, γ, ω, ψ, χ, pop,
    _...,
) = begin
    rhEI = relhaz(t,guide,node,Expos,Infec,cols)
    rhIE = relhaz(t,guide,node,Infec,Expos,cols)
    a = β*S*I/pop
    preboost[1] = I > 0 ? 1-ellI/I : 0
    alpha[1] = a*preboost[1]
    preboost[2] = I > 0 ? sum(rhIE)/I : 1
    alpha[2] = a*preboost[2]
    b = σ*E
    preboost[3] = E > 0 ? 1-ellE/E : 0
    alpha[3] = b*preboost[3]
    preboost[4] = E > 0 ? sum(rhEI)/E : 1
    alpha[4] = b*preboost[4]
    c = γ*I
    preboost[5] = 1-ellI/I
    alpha[5] = @indicator(I > ellI, c*preboost[5])
    preboost[6] = 1.0
    alpha[6] = ω*R
    decay = ψ*I + χ*I +
        a - alpha[1] - alpha[2] +
        b - alpha[3] - alpha[4] +
        c - alpha[5]
    decay, rhEI, rhIE
end

regular_part!(
    cols, guide, node, ll,
    S, E, I, R;
    kwargs...,
) = begin
    n = guide[node]
    t = n.tbeg
    tf = n.tend
    if t < tf
        alpha = similar(Vector{Prob},6)
        preboost = similar(Vector{Prob},6)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellE, ellI = ell(cols)
        while t < tf
            decay,rhEI,rhIE = event_rates!(
                alpha, preboost, t, guide, node,
                cols, ellE, ellI,
                S, E, I, R;
                kwargs...
            )
            k, s = rcateg(alpha)
            step = -log(rand())/s
            if t+step < tf
                ll -= decay*step + log(preboost[k])
                if k==1
                    S -= 1
                    E += 1
                    ll += log(1-ellE/E)
                elseif k==2
                    b,_,p = rcateg(rhIE,cols[Infec],true)
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
                    b,_,p = rcateg(rhEI,cols[Expos],true)
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
    end
    ll, S, E, I, R
end

end
