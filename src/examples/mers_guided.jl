"""
    GuidedMERS

A module containing an implementation of the phylodynamic filter for the
MERS-CoV two-host (Camel/Human) model, using "guided" proposals. Unlike
`SoftMERS`, which preserves the naive identity/tracked-branch mass split
and only guides the choice *within* the tracked-branch group, `GuidedMERS`
jointly normalizes the identity weight (the untracked source-host count)
together with every tracked branch's relative hazard, then draws from
a single categorical over all of them at once. The identity/tracked split
itself is therefore guide-driven here, not fixed to the naive ratio --
this is what makes `GuidedMERS` mathematically distinct from `SoftMERS`.

This is the MERS analogue of `GuidedSEIR` (seir_guided.jl); the joint
chooser uses `choose_branch` from `guide.jl`.
"""
module GuidedMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

const mers_tree = parse_newick(first(mers_trees), t0=0, demes=Demes)

include("mers_funs.jl")

event_rates!(
    alpha, ellC, ellH,
    S_c, I_c, S_h, I_h;
    β_cc, β_ch, β_hc, β_hh,
    γ_c, γ_h, χ_c, χ_h, B_c, B_h, N_c, N_h,
    _...,
) = begin
    alpha[1]  = β_cc*S_c*I_c/N_c                              # cc transmission
    alpha[2]  = β_hh*S_h*I_h/N_h                              # hh transmission
    alpha[3]  = β_hc*S_h*I_c/N_c                              # camel → human
    alpha[4]  = β_ch*S_c*I_h/N_h                              # human → camel
    alpha[5]  = @indicator(I_c > ellC, γ_c*(I_c-ellC))     # camel removal
    alpha[6]  = @indicator(I_h > ellH, γ_h*(I_h-ellH))     # human removal
    alpha[7]  = B_c                                            # S_c birth
    alpha[8]  = B_h                                            # S_h birth
    alpha[9]  = B_c*S_c/N_c                                    # S_c death
    alpha[10] = B_h*S_h/N_h                                    # S_h death
    χ_c*I_c + χ_h*I_h +
        γ_c*ellC + @indicator(I_c ≤ ellC, γ_c*(I_c-ellC)) +
        γ_h*ellH + @indicator(I_h ≤ ellH, γ_h*(I_h-ellH))
end

regular_part!(
    cols, guide, node, ll,
    S_c, I_c, S_h, I_h;
    kwargs...,
) = begin
    n = guide[node]
    t = n.tbeg
    tf = n.tend
    if t < tf
        alpha = similar(Vector{Prob}, 10)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellC, ellH = ell(cols)
        while t < tf
            decay = event_rates!(
                alpha, ellC, ellH,
                S_c, I_c, S_h, I_h;
                kwargs...,
            )
            k, s = rcateg(alpha)
            step = -log(rand())/s
            if t+step < tf
                ll -= decay*step
                if k==1                                     # cc
                    S_c -= 1
                    I_c += 1
                    ll += log(1 - ellC*(ellC-1)/I_c/(I_c-1))
                elseif k==2                                 # hh
                    S_h -= 1
                    I_h += 1
                    ll += log(1 - ellH*(ellH-1)/I_h/(I_h-1))
                elseif k==3                                 # camel → human
                    b, p = choose_branch(t, guide, node, I_c, cols, Camel, Human)
                    ll -= log(p)
                    S_h -= 1
                    I_h += 1
                    if b == 0
                        ll += log(1 - ellH/I_h)
                    else
                        ellC, ellH = swap!(cols, Camel, Human, b)
                        ll += log(1 - ellC/I_c) - log(I_h)
                    end
                elseif k==4                                 # human → camel
                    b, p = choose_branch(t, guide, node, I_h, cols, Human, Camel)
                    ll -= log(p)
                    S_c -= 1
                    I_c += 1
                    if b == 0
                        ll += log(1 - ellC/I_c)
                    else
                        ellC, ellH = swap!(cols, Human, Camel, b)
                        ll += log(1 - ellH/I_h) - log(I_c)
                    end
                elseif k==5                                 # camel removal
                    ll -= log(1 - ellC/I_c)
                    I_c -= 1
                elseif k==6                                 # human removal
                    ll -= log(1 - ellH/I_h)
                    I_h -= 1
                elseif k==7                                 # S_c birth
                    S_c += 1
                elseif k==8                                 # S_h birth
                    S_h += 1
                elseif k==9                                 # S_c death
                    S_c -= 1
                elseif k==10                                # S_h death
                    S_h -= 1
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
        @assert I_c ≥ ellC && I_h ≥ ellH
    end
    ll, S_c, I_c, S_h, I_h
end

end
