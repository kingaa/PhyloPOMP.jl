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

This is the MERS analogue of `SoftSEIR` (seir_soft.jl), and is
mathematically distinct from `GuidedMERS` (mers_guided.jl), which instead
jointly normalizes the identity weight and the branch hazards together
(so the identity/tracked split itself, not just the within-group choice,
is guide-driven).
"""
module SoftMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

const mers_tree = parse_newick(first(mers_trees), t0=0, demes=Demes)

include("mers_funs.jl")

event_rates!(
    alpha, pi, cols,
    S_c, I_c, S_h, I_h;
    β_cc, β_ch, β_hc, β_hh,
    γ_c, γ_h, χ_c, χ_h, B_c, B_h, N_c, N_h,
    _...,
) = begin
    ellC, ellH = ell(cols)
    @assert I_c ≥ ellC && I_h ≥ ellH
    alpha[1] = β_cc*S_c*I_c/N_c
    alpha[2] = β_hh*S_h*I_h/N_h
    alpha[3] = alpha[4] = β_hc*S_h*I_c/N_c
    alpha[5] = alpha[6] = β_ch*S_c*I_h/N_h
    alpha[7] = @indicator(I_c > ellC, γ_c*(I_c-ellC))
    alpha[8] = @indicator(I_h > ellH, γ_h*(I_h-ellH))
    alpha[9] = B_c
    alpha[10] = B_h
    alpha[11] = B_c*S_c/N_c
    alpha[12] = B_h*S_h/N_h

    pi[1:2] .= one(Prob)
    pi[3] = @indicator(I_c > 0, (I_c-ellC)/I_c)
    pi[4] = @indicator(I_c > 0, ellC/I_c)
    pi[5] = @indicator(I_h > 0, (I_h-ellH)/I_h)
    pi[6] = @indicator(I_h > 0, ellH/I_h)
    pi[7:12] .= one(Prob)

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
        alpha = similar(Vector{Prob}, 12)
        pi = similar(Vector{Prob}, 12)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellC, ellH = ell(cols)
        while t < tf
            decay = event_rates!(
                alpha, pi, cols,
                S_c, I_c, S_h, I_h;
                kwargs...,
            )
            k, s = rcateg(alpha .* pi)
            step = -log(rand())/s
            if t+step < tf
                ll -= decay*step + log(pi[k])
                if k==1
                    S_c -= 1
                    I_c += 1
                    ll += log(1 - ellC*(ellC-1)/I_c/(I_c-1))
                elseif k==2
                    S_h -= 1
                    I_h += 1
                    ll += log(1 - ellH*(ellH-1)/I_h/(I_h-1))
                elseif k==3                                 # camel → human, identity
                    S_h -= 1
                    I_h += 1
                    ll += log(1 - ellH/I_h)
                elseif k==4                                 # camel → human, tracked-branch swap
                    lins = [guide[node].linmap[b] for b ∈ cols[Camel]]
                    rh = relhaz(t, guide, node, Camel, Human, lins)
                    b, _, share = rcateg(rh, cols[Camel], true)
                    ll -= log(share)
                    S_h -= 1
                    I_h += 1
                    ellC, ellH = swap!(cols, Camel, Human, b)
                    ll += log(1 - ellC/I_c) - log(I_h)
                elseif k==5                                 # human → camel, identity
                    S_c -= 1
                    I_c += 1
                    ll += log(1 - ellC/I_c)
                elseif k==6                                 # human → camel, tracked-branch swap
                    lins = [guide[node].linmap[b] for b ∈ cols[Human]]
                    rh = relhaz(t, guide, node, Human, Camel, lins)
                    b, _, share = rcateg(rh, cols[Human], true)
                    ll -= log(share)
                    S_c -= 1
                    I_c += 1
                    ellC, ellH = swap!(cols, Human, Camel, b)
                    ll += log(1 - ellH/I_h) - log(I_c)
                elseif k==7                                 # camel removal
                    ll -= log(1 - ellC/I_c)
                    I_c -= 1
                elseif k==8                                 # human removal
                    ll -= log(1 - ellH/I_h)
                    I_h -= 1
                elseif k==9                                 # S_c birth
                    S_c += 1
                elseif k==10                                # S_h birth
                    S_h += 1
                elseif k==11                                # S_c death
                    S_c -= 1
                elseif k==12                                # S_h death
                    S_h -= 1
                else
                    @assert false "impossible event" # COV_EXCL_LINE
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
