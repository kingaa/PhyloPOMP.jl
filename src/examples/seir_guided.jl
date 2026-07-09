"""
    GuidedSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using a filter guide and "soft" proposals.  That is, when
population-process events occur, the guide can steer those events
preferentially onto (or away from) particular branches, but the overall
event rate remains equal to that of the underlying population process.
"""
module GuidedSEIR

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Expos Infec
using .Demes: Expos, Infec, DemeSet

include("seir_trees.jl")
include("seir_funs.jl")

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
            alpha, pi;
            S, E, I, R,
            ellE, ellI,
            kwargs...,
            onE=ellE,
            onI=ellI,
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
                relhaz!(rh,t,guide,node)
                b, p = choose_branch(rh,n,cols,Infec,Expos)
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
                relhaz!(rh,t,guide,node)
                b, p = choose_branch(rh,n,cols,Expos,Infec)
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
