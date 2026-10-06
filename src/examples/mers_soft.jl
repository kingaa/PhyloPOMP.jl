"""
    SoftMERS

Filter for the two-host MERS model with soft proposals.
Population-event rates equal the model's; the guide only picks among tracked branches.
"""
module SoftMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

const mers_tree = parse_newick(mers_newick, t0=0, demes=Demes)

include("mers_funs.jl")

## Soft: onC = I_c-offC, so pi[3]+pi[4] = 1 and transmission adds no decay.
## The guide enters only through choose_branch.

regular_part!(
    cols, state, guide, n, t, tf;
    kwargs...,
) = begin
    node = guide[n]
    (;S_c, I_c, S_h, I_h) = state
    alpha = similar(Vector{Prob},12)
    pi = similar(Vector{Prob},12)
    rh = relhaz_alloc(guide,n)
    step::Time = zero(Time)
    decay::Prob = zero(Prob)
    ll::Prob = zero(Prob)
    ellC, ellH = ell(cols)
    while t < tf
        offC = I_c*no_move_share(ellC,ellH,I_c,I_h)
        offH = I_h*no_move_share(ellH,ellC,I_h,I_c)
        decay = event_rates!(
            alpha, pi, (;S_c, I_c, S_h, I_h), cols;
            kwargs...,
            onC=I_c-offC, offC=offC,
            onH=I_h-offH, offH=offH,
        )
        k, s = rcateg(alpha)
        step = -log(rand())/s
        if k > 0 && t+step < tf
            ll -= decay*step + log(pi[k])
            if k==1                             # camel → camel
                S_c -= 1
                I_c += 1
                ## no branch point is observed here, so the birth must
                ## not have joined two tracked camel lineages
                ll += log(1 - ellC*(ellC-1)/I_c/(I_c-1))
            elseif k==2                         # human → human
                S_h -= 1
                I_h += 1
                ll += log(1 - ellH*(ellH-1)/I_h/(I_h-1))
            elseif k==3                         # camel → human, identity
                S_h -= 1
                I_h += 1
                ll += log(1 - ellH/I_h)
            elseif k==4                         # camel → human, tracked-branch swap
                relhaz!(rh,t,guide,n)
                b, p = choose_branch(rh,node,cols,Camel,Human)
                ll -= log(p)
                S_h -= 1
                I_h += 1
                ellC, ellH = swap!(cols,Camel,Human,b)
                ll += log(1 - ellC/I_c) - log(I_h)
            elseif k==5                         # human → camel, identity
                S_c -= 1
                I_c += 1
                ll += log(1 - ellC/I_c)
            elseif k==6                         # human → camel, tracked-branch swap
                relhaz!(rh,t,guide,n)
                b, p = choose_branch(rh,node,cols,Human,Camel)
                ll -= log(p)
                S_c -= 1
                I_c += 1
                ellC, ellH = swap!(cols,Human,Camel,b)
                ll += log(1 - ellH/I_h) - log(I_c)
            elseif k==7                         # camel removal
                I_c -= 1
            elseif k==8                         # human removal
                I_h -= 1
            elseif k==9                         # S_c birth
                S_c += 1
            elseif k==10                        # S_h birth
                S_h += 1
            elseif k==11                        # S_c death
                S_c -= 1
            elseif k==12                        # S_h death
                S_h -= 1
            end
            t += step
        else
            step = tf - t
            ll -= decay*step
            break
        end
    end
    @assert I_c ≥ ellC && I_h ≥ ellH
    ll, (;S_c, I_c, S_h, I_h)
end

end
