"""
    NaiveMTBD

Genealogy-conditioned particle filter for the two-type linear
birth-death-sampling (MTBD) model, with naive proposals: when an
unobserved event could either leave the tracked lineages alone or move
one of them to the other deme, the choice is made without looking at
the tip types downstream.

The filter degenerates on large trees.  See [`GuidedMTBD`](@ref).
"""
module NaiveMTBD

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mtbd_funs.jl")

singular_part!(
    cols, ll, geneal, node,
    I1, I2;
    lambda11, lambda12, lambda21, lambda22,
    psi1, psi2,
    _...,
) = begin
    ell1, ell2 = ell(cols)
    n = geneal[node]
    @assert I1 ≥ ell1 && I2 ≥ ell2
    if n.type == Root
        ## the root lineage is carried by an untracked host, chosen at random
        i, _, p = rcateg([I1-ell1, I2-ell2], DemeSet, true)
        if ismissing(i)
            ll = Prob(-Inf)
            ell1, ell2 = plant!(cols, Camel, n.lineage)
            I1 += 1
        else
            ll -= log(p)
            ell1, ell2 = plant!(cols, i, n.lineage)
        end
    elseif n.type == Sample
        if n.lineage ∉ cols[n.deme]
            ## the tip type contradicts this particle's colouring;
            ## recolour so the state stays valid, and kill the particle
            ll = Prob(-Inf)
            if n.deme == Camel
                ell1, ell2 = swap!(cols, Human, Camel, n.lineage)
                I1 += 1
            else
                ell1, ell2 = swap!(cols, Camel, Human, n.lineage)
                I2 += 1
            end
        end
        ell1, ell2 = chop!(cols, n.deme, n.lineage)
        ## destructive sampling (r = 1): total sampling rate in the deme
        if n.deme == Camel
            ll += log(psi1*I1)
            I1 -= 1
        else
            ll += log(psi2*I2)
            I2 -= 1
        end
    elseif n.type == Node
        children = map(c -> geneal[c].lineage, n.children)
        if n.lineage ∈ cols[Camel]
            a11 = lambda11*I1          # total 1 -> 1 birth rate
            a12 = lambda12*I1          # total 1 -> 2 birth rate
            k, _, p = rcateg([a11, 0.5*a12, 0.5*a12], true)
            if k == 0
                ll = Prob(-Inf)
                ell1, ell2 = fork!(cols, Camel, n.lineage, (Camel, Camel), children)
                I1 += 1
            else
                ll -= log(p)
                if k == 1
                    ell1, ell2 = fork!(cols, Camel, n.lineage, (Camel, Camel), children)
                    I1 += 1
                    ll += log(a11) - log(I1*(I1-1)/2)
                else
                    orient = k == 2 ? (Camel, Human) : (Human, Camel)
                    ell1, ell2 = fork!(cols, Camel, n.lineage, orient, children)
                    I2 += 1
                    ll += log(a12) - log(I1*I2)
                end
            end
        elseif n.lineage ∈ cols[Human]
            a22 = lambda22*I2
            a21 = lambda21*I2
            k, _, p = rcateg([a22, 0.5*a21, 0.5*a21], true)
            if k == 0
                ll = Prob(-Inf)
                ell1, ell2 = fork!(cols, Human, n.lineage, (Human, Human), children)
                I2 += 1
            else
                ll -= log(p)
                if k == 1
                    ell1, ell2 = fork!(cols, Human, n.lineage, (Human, Human), children)
                    I2 += 1
                    ll += log(a22) - log(I2*(I2-1)/2)
                else
                    orient = k == 2 ? (Camel, Human) : (Human, Camel)
                    ell1, ell2 = fork!(cols, Human, n.lineage, orient, children)
                    I1 += 1
                    ll += log(a21) - log(I1*I2)
                end
            end
        else
            @assert false "impossible node deme" # COV_EXCL_LINE
        end
    else
        @assert false "impossible node type" # COV_EXCL_LINE
    end
    @assert I1 ≥ ell1 && I2 ≥ ell2
    ll, I1, I2
end

## Rates of unobserved events and the naive split pi.
## A cross-type birth proposes "no move" with share 1/(ell+1),
## so it stays possible when ell = I.
## A migration carries its lineage, so "no move" needs an untracked migrant.
event_rates!(
    alpha, pi, cols,
    I1, I2;
    lambda11, lambda12, lambda21, lambda22,
    mu1, mu2, psi1, psi2, m12, m21,
    _...,
) = begin
    ell1, ell2 = ell(cols)
    @assert I1 ≥ ell1 && I2 ≥ ell2
    alpha[1] = lambda11*I1                          # 1 -> 1 birth
    alpha[2] = lambda22*I2                          # 2 -> 2 birth
    alpha[3] = alpha[4] = lambda12*I1               # 1 -> 2 birth
    alpha[5] = alpha[6] = lambda21*I2               # 2 -> 1 birth
    alpha[7] = @indicator(I1 > ell1, mu1*(I1-ell1)) # untracked death, 1
    alpha[8] = @indicator(I2 > ell2, mu2*(I2-ell2)) # untracked death, 2
    alpha[9] = alpha[10] = m12*I1                   # migration 1 -> 2
    alpha[11] = alpha[12] = m21*I2                  # migration 2 -> 1

    pi[1:2] .= one(Prob)
    pi[3] = 1/(ell1+1)                              # 1 -> 2, no move
    pi[4] = ell1/(ell1+1)                           # 1 -> 2, move
    pi[5] = 1/(ell2+1)                              # 2 -> 1, no move
    pi[6] = ell2/(ell2+1)                           # 2 -> 1, move
    pi[7:8] .= one(Prob)
    pi[9] = @indicator(I1 > ell1, 1/(ell1+1))       # migration 1 -> 2, no move
    pi[10] = one(Prob) - pi[9]                      # migration 1 -> 2, move
    pi[11] = @indicator(I2 > ell2, 1/(ell2+1))      # migration 2 -> 1, no move
    pi[12] = one(Prob) - pi[11]                     # migration 2 -> 1, move

    ## decay: tracked hosts must survive, and nothing is sampled off the tree
    psi1*I1 + psi2*I2 + mu1*ell1 + mu2*ell2
end

regular_part!(
    cols, ll,
    t, dt,
    I1, I2;
    args...,
) = begin
    tf = t + dt
    if t < tf
        alpha = similar(Vector{Prob}, 12)
        pi = similar(Vector{Prob}, 12)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ell1, ell2 = ell(cols)
        while t < tf
            decay = event_rates!(alpha, pi, cols, I1, I2; args...)
            k, s = rcateg(alpha .* pi)
            step = -log(rand())/s
            if t+step < tf
                ll -= decay*step + log(pi[k])
                if k == 1                               # 1 -> 1 birth
                    I1 += 1
                    ## no branch point is observed, so the birth must not
                    ## have joined two tracked lineages
                    ll += log(1-ell1*(ell1-1)/I1/(I1-1))
                elseif k == 2                           # 2 -> 2 birth
                    I2 += 1
                    ll += log(1-ell2*(ell2-1)/I2/(I2-1))
                elseif k == 3                           # 1 -> 2, no move
                    I2 += 1
                    ll += log(1-ell2/I2)
                elseif k == 4                           # 1 -> 2, move
                    ll += log(ell1)
                    I2 += 1
                    ell1, ell2 = swap!(cols, Camel, Human, rand(cols[Camel]))
                    ll += log((1-ell1/I1)/I2)
                elseif k == 5                           # 2 -> 1, no move
                    I1 += 1
                    ll += log(1-ell1/I1)
                elseif k == 6                           # 2 -> 1, move
                    ll += log(ell2)
                    I1 += 1
                    ell1, ell2 = swap!(cols, Human, Camel, rand(cols[Human]))
                    ll += log((1-ell2/I2)/I1)
                elseif k == 7                           # untracked death, 1
                    ll -= log(1-ell1/I1)
                    I1 -= 1
                elseif k == 8                           # untracked death, 2
                    ll -= log(1-ell2/I2)
                    I2 -= 1
                elseif k == 9                           # migration 1 -> 2, no move
                    I1 -= 1
                    I2 += 1
                    ll += log(1-ell2/I2)
                elseif k == 10                          # migration 1 -> 2, move
                    ll += log(ell1)
                    ell1, ell2 = swap!(cols, Camel, Human, rand(cols[Camel]))
                    I1 -= 1
                    I2 += 1
                    ll -= log(I2)
                elseif k == 11                          # migration 2 -> 1, no move
                    I2 -= 1
                    I1 += 1
                    ll += log(1-ell1/I1)
                elseif k == 12                          # migration 2 -> 1, move
                    ll += log(ell2)
                    ell1, ell2 = swap!(cols, Human, Camel, rand(cols[Human]))
                    I2 -= 1
                    I1 += 1
                    ll -= log(I1)
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
        @assert I1 ≥ ell1 && I2 ≥ ell2
    end
    ll, I1, I2
end

"""
    filter_pomp(gen; lambda11, lambda12, lambda21, lambda22,
                mu1, mu2, psi1, psi2, m12, m21, I1_0 = 1, I2_0 = 0)

Construct the naive MTBD filter for genealogy `gen`, which should be
parsed with `demes = NaiveMTBD.Demes`.  Rates default to the MERS
maximum-likelihood estimate ([`mers_mle`](@ref)); see
[`epi_params`](@ref) to build them from reproduction numbers.  The
process starts from `I1_0` type-1 and `I2_0` type-2 hosts.
"""
filter_pomp(
    gen::Genealogy;
    lambda11 = mers_mle.lambda11, lambda12 = mers_mle.lambda12,
    lambda21 = mers_mle.lambda21, lambda22 = mers_mle.lambda22,
    mu1 = mers_mle.mu1, mu2 = mers_mle.mu2,
    psi1 = mers_mle.psi1, psi2 = mers_mle.psi2,
    m12 = mers_mle.m12, m21 = mers_mle.m21,
    I1_0 = 1, I2_0 = 0,
) = begin
    check(gen)
    pomp(
        params = (
            lambda11 = Float64(lambda11), lambda12 = Float64(lambda12),
            lambda21 = Float64(lambda21), lambda22 = Float64(lambda22),
            mu1 = Float64(mu1), mu2 = Float64(mu2),
            psi1 = Float64(psi1), psi2 = Float64(psi2),
            m12 = Float64(m12), m21 = Float64(m21),
            I1_0 = Float64(I1_0), I2_0 = Float64(I2_0),
        ),
        t0 = timezero(gen),
        times = times(gen),
        init_state = (
            node = one(Name), ll = zero(Prob),
            cols = Coloring(Demes), I1 = zero(Int64), I2 = zero(Int64),
        ),
        rinit = function (; kwargs...)
            (; node = one(Name), ll = zero(Prob),
             cols = Coloring(Demes), mtbd_rinit(; kwargs...)...)
        end,
        rprocess = onestep(
            function (; node, ll, cols, geneal, I1, I2, t, dt, args...)
                cols = copy(cols)
                ll = zero(Prob)
                ll, I1, I2 = singular_part!(cols, ll, geneal, node, I1, I2; args...)
                if isfinite(ll)
                    ll, I1, I2 = regular_part!(cols, ll, t, dt, I1, I2; args...)
                end
                (; node = node + one(Name), ll, cols, I1, I2)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,),
    )
end

end
