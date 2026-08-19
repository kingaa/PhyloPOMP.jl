"""
    HardMERS

A module containing an implementation of the phylodynamic filter for the
MERS-CoV two-host (Camel/Human) model, using a filter guide and "hard"
proposals. That is, color changes can be proposed on branches at
*unnormalized* auxiliary intensities that may exceed (or fall short of) the
event rates in the underlying population process; any rate not consumed by
these auxiliary intensities is returned to the population process through
the compensating decay `alpha_population - sum_j beta_j`.

This is the MERS analogue of `HardSEIR` (seir_hard.jl), and differs from
`SoftMERS` precisely in that its per-branch auxiliary intensities use the
raw relative hazards directly rather than renormalizing them to preserve
the naive tracked-branch mass `ell/I`.
"""
module HardMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

const mers_tree = parse_newick(first(mers_trees), t0=0, demes=Demes)

include("mers_funs.jl")

event_rates!(
    alpha, pi, rh, n,
    cols, ellC, ellH,
    S_c, I_c, S_h, I_h;
    β_cc, β_ch, β_hc, β_hh,
    γ_c, γ_h, χ_c, χ_h, B_c, B_h, N_c, N_h,
    _...,
) = begin
    rh_c2h = relhaz(rh,n,cols,Camel,Human)
    rh_h2c = relhaz(rh,n,cols,Human,Camel)

    rate_cc = β_cc*S_c*I_c/N_c
    pi[1] = 1.0
    alpha[1] = rate_cc

    rate_hh = β_hh*S_h*I_h/N_h
    pi[2] = 1.0
    alpha[2] = rate_hh

    rate_hc = β_hc*S_h*I_c/N_c                                   # camel → human
    pi[3] = @indicator(I_c > 0, (I_c-ellC)/I_c)          # no-swap (identity)
    alpha[3] = rate_hc*pi[3]
    onC = sum(rh_c2h)
    pi[4] = @indicator(I_c > 0, onC/I_c)                 # swap (aggregate)
    alpha[4] = rate_hc*pi[4]

    rate_ch = β_ch*S_c*I_h/N_h                                   # human → camel
    pi[5] = @indicator(I_h > 0, (I_h-ellH)/I_h)          # no-swap (identity)
    alpha[5] = rate_ch*pi[5]
    onH = sum(rh_h2c)
    pi[6] = @indicator(I_h > 0, onH/I_h)                 # swap (aggregate)
    alpha[6] = rate_ch*pi[6]

    rate_c = γ_c*I_c                                          # camel removal
    pi[7] = @indicator(I_c > 0, 1-ellC/I_c)
    alpha[7] = @indicator(I_c > ellC, rate_c*pi[7])

    rate_h = γ_h*I_h                                          # human removal
    pi[8] = @indicator(I_h > 0, 1-ellH/I_h)
    alpha[8] = @indicator(I_h > ellH, rate_h*pi[8])

    pi[9]  = 1.0; alpha[9]  = B_c                        # S_c birth
    pi[10] = 1.0; alpha[10] = B_h                        # S_h birth
    pi[11] = 1.0; alpha[11] = B_c*S_c/N_c                        # S_c death
    pi[12] = 1.0; alpha[12] = B_h*S_h/N_h                        # S_h death

    decay = χ_c*I_c + χ_h*I_h +
        rate_hc - alpha[3] - alpha[4] +
        rate_ch - alpha[5] - alpha[6] +
        rate_c - alpha[7] +
        rate_h - alpha[8]
    decay, rh_c2h, rh_h2c
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
        rh = relhaz_alloc(guide,node)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellC, ellH = ell(cols)
        while t < tf
            relhaz!(rh,t,guide,node)
            decay, rh_c2h, rh_h2c = event_rates!(
                alpha, pi, rh, n,
                cols, ellC, ellH,
                S_c, I_c, S_h, I_h;
                kwargs...,
            )
            k, s = rcateg(alpha)
            step = -log(rand())/s
            if t+step < tf
                ll -= decay*step + log(pi[k])
                if k==1                                     # cc
                    S_c -= 1
                    I_c += 1
                    ll += log(1 - ellC*(ellC-1)/I_c/(I_c-1))
                elseif k==2                                 # hh
                    S_h -= 1
                    I_h += 1
                    ll += log(1 - ellH*(ellH-1)/I_h/(I_h-1))
                elseif k==3                                 # camel → human, no-swap
                    S_h -= 1
                    I_h += 1
                    ll += log(1 - ellH/I_h)
                elseif k==4                                 # camel → human, swap
                    b, _, p = rcateg(rh_c2h, cols[Camel], true)
                    ll -= log(p)
                    S_h -= 1
                    I_h += 1
                    ellC, ellH = swap!(cols, Camel, Human, b)
                    ll += log(1 - ellC/I_c) - log(I_h)
                elseif k==5                                 # human → camel, no-swap
                    S_c -= 1
                    I_c += 1
                    ll += log(1 - ellC/I_c)
                elseif k==6                                 # human → camel, swap
                    b, _, p = rcateg(rh_h2c, cols[Human], true)
                    ll -= log(p)
                    S_c -= 1
                    I_c += 1
                    ellC, ellH = swap!(cols, Human, Camel, b)
                    ll += log(1 - ellH/I_h) - log(I_c)
                elseif k==7                                 # camel removal
                    I_c -= 1
                elseif k==8                                 # human removal
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
