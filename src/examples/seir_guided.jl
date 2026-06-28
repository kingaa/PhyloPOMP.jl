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

event_rates!(
    alpha, ellI,
    S, E, I, R;
    β, σ, γ, ω, ψ, χ, pop,
    _...,
) = begin
    alpha[1] = β*S*I/pop
    alpha[2] = σ*E
    alpha[3] = @indicator(I > ellI, γ*(I-ellI))
    alpha[4] = ω*R
    ψ*I + χ*I + γ*ellI + @indicator(I ≤ ellI, γ*(I-ellI))
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
        alpha = similar(Vector{Prob},4)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellE, ellI = ell(cols)
        while t < tf
            decay = event_rates!(
                alpha, ellI,
                S, E, I, R;
                kwargs...,
            )
            k, s = rcateg(alpha)
            step = -log(rand())/s
            if t+step < tf
                ll -= decay*step
                if k==1
                    b,p = choose_branch(t,guide,node,I,cols,Infec,Expos)
                    ll -= log(p)
                    S -= 1
                    E += 1
                    if b == 0
                        ll += log(1-ellE/E)
                    else
                        ellE, ellI = swap!(cols,Infec,Expos,b)
                        ll += log(1-ellI/I)-log(E)
                    end
                elseif k==2
                    b,p = choose_branch(t,guide,node,E,cols,Expos,Infec)
                    ll -= log(p)
                    E -= 1
                    I += 1
                    if b == 0
                        ll += log(1-ellI/I)
                    else
                        ellE, ellI = swap!(cols,Expos,Infec,b)
                        ll -= log(I)
                    end
                elseif k==3
                    ll -= log(1-ellI/I)
                    I -= 1
                    R += 1
                elseif k==4
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
