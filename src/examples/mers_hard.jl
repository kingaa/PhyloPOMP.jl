"""
    HardMERS

A module containing an implementation of the phylodynamic filter for the
MERS-CoV two-host (Camel/Human) model, using a filter guide and "hard"
proposals. That is, color changes can be proposed on branches at
*unnormalized* auxiliary intensities that may exceed (or fall short of) the
event rates in the underlying population process; any rate not consumed by
these auxiliary intensities is returned to the population process through
the compensating decay `alpha_population - sum_j beta_j`.

This is the MERS analogue of `HardSEIR` (seir_hard.jl). It shares its
singular (genealogical-event) part, its rate functions and its
`filter_pomp` with `SoftMERS` via `mers_funs.jl`, which in turn carries
those pieces verbatim from `GuidedMERS` (mers_guided.jl); the three
modules differ *only* in `regular_part!`.

It differs from `SoftMERS` precisely in that its per-branch auxiliary
intensities use the raw relative hazards directly (`onC =
sum_relhaz(...)`) rather than renormalizing them to preserve the naive
tracked-branch mass `ell/I`; the shortfall or excess is absorbed by the
decay term returned from `transmission!`. All three kernels target the
same likelihood.
"""
module HardMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

const mers_tree = parse_newick(mers_newick, t0=0, demes=Demes)

include("mers_funs.jl")

## Hard proposals: the aggregate tracked-branch intensity is the *raw*
## sum of relative hazards, `onC = sum_relhaz(rh,node,cols,Camel,Human)`,
## while the identity share keeps its naive value, `offC = I_c-ellC`.
## `pi[3]+pi[4]` therefore need not equal 1, and the leftover
## `rate_hc - alpha[3] - alpha[4]` (of either sign) is returned by
## `transmission!` as decay.  `relhaz!` must run once per loop iteration,
## before `event_rates!`, because the rates themselves depend on it.
##
## Removal accounting is unchanged from the former on-disk version of
## this file: `pi[7] = 1-ellC/I_c`, `alpha[7] = γ_c*I_c*pi[7]` gated on
## `I_c > ellC`, and the generic `ll -= decay*step + log(pi[k])` line
## charges `-log(1-ellC/I_c)` when a removal fires.  (See mers_soft.jl
## for the proof that SoftMERS's former, differently-bookkept, removal
## convention is identical in value to this one.)

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
        relhaz!(rh,t,guide,n)
        decay = event_rates!(
            alpha, pi, (;S_c, I_c, S_h, I_h), cols;
            kwargs...,
            onC=sum_relhaz(rh,node,cols,Camel,Human),
            offC=I_c-ellC,
            onH=sum_relhaz(rh,node,cols,Human,Camel),
            offH=I_h-ellH,
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
