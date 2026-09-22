"""
    SoftMERS

A module containing an implementation of the phylodynamic filter for the
MERS-CoV two-host (Camel/Human) model, using so-called "soft" proposals.
Soft preserves the naive aggregate identity/tracked-branch mass split
(`(I-ell)/I` vs `ell/I`) but, when a tracked-branch event is proposed,
distributes that fixed mass among the tracked branches according to the
guide's relative hazards rather than uniformly. The overall event rate
for each population process therefore remains exactly that of the
underlying population process, matching the target's collapsed
identity/tracked split; only the *within-group* branch choice is guided.

This is the MERS analogue of `SoftSEIR` (seir_soft.jl). It shares its
singular (genealogical-event) part, its rate functions and its
`filter_pomp` with `HardMERS` via `mers_funs.jl`, which in turn carries
those pieces verbatim from `GuidedMERS` (mers_guided.jl); the three
modules differ *only* in `regular_part!`.

`GuidedMERS` differs from `SoftMERS` in exactly one respect: its
cross-deme proposals call `choose_branch(t, guide, node, I, cols, i, j)`,
which normalizes the identity weight `I-ell` jointly against the *raw*
per-branch relative hazards. The identity/tracked split is therefore
itself guide-driven there, whereas here it is pinned at `(I-ell)/I` vs
`ell/I` and only the choice *within* the tracked group is guided. All
three kernels target the same likelihood.
"""
module SoftMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

const mers_tree = parse_newick(mers_newick, t0=0, demes=Demes)

include("mers_funs.jl")

## Soft proposals: the aggregate tracked-branch mass is `ell/I`, i.e.
## exactly what the naive (unguided) filter uses, so `onC = ellC` and
## `onH = ellH`; `offC`/`offH` then default to `I-ell`.  Since
## `pi[3]+pi[4] = pi[5]+pi[6] = 1`, the transmission processes contribute
## nothing to the decay, and the guide enters only through
## `choose_branch`, which apportions the tracked mass among the tracked
## branches in proportion to their relative hazards.
##
## Removal accounting: the former on-disk version of this file used
## `alpha[7] = γ_c*(I_c-ellC)` with `pi[7] = 1` and then charged
## `-log(1-ellC/I_c)` by hand at the removal event.  The shared
## `removal!` in mers_funs.jl instead sets `pi[7] = 1-ellC/I_c` and
## `alpha[7] = γ_c*I_c*pi[7]`, letting the generic `ll -= log(pi[k])`
## line do the charging.  The two conventions are identical in value:
##   * alpha:  γ_c*I_c*(1-ellC/I_c) = γ_c*(I_c-ellC), in both cases
##             gated on `I_c > ellC`;
##   * charge: -log(pi[7]) = -log(1-ellC/I_c), the same term;
##   * decay:  if I_c > ellC, `rate_c - alpha[7]` = γ_c*ellC, which is
##             the old `γ_c*ellC` term; if I_c ≤ ellC, alpha[7] = 0 and
##             `rate_c - alpha[7]` = γ_c*I_c, which is the old
##             `γ_c*ellC + γ_c*(I_c-ellC)`.
## The same argument holds slot-for-slot for human removal.

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
        decay = event_rates!(
            alpha, pi, (;S_c, I_c, S_h, I_h), cols;
            kwargs...,
            onC=ellC,
            onH=ellH,
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
