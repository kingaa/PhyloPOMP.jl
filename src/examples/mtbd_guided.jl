"""
    GuidedMTBD

Genealogy-conditioned particle filter for the two-type linear
birth-death-sampling (MTBD) model, with guided proposals.  The guide
(see [`guide`](@ref)) carries information about tip types back up the
tree, and is used to steer three kinds of choice:

- the deme of the root lineage;
- the orientation of a cross-type branch point;
- whether an unobserved cross-type birth or migration moves one of the
  tracked lineages, and if so which one.

Event rates themselves are those of the model ("semisoft" proposals,
as in [`GuidedMERS`](@ref)), so only the colouring is guided.  Each
guided choice is corrected by its proposal probability, so the filter
estimates the same likelihood as [`NaiveMTBD`](@ref), with less
variance when tip types are informative.
"""
module GuidedMTBD

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mtbd_funs.jl")

"""
    knowledge!(v; deme, type, time)

The deme is known exactly where the genealogy supplies it: at sample
tips.  Branch points are left free, since in MTBD a branch point can be
within or across types.
"""
knowledge!(v; deme, type, time) = begin
    if ismissing(deme)
        false
    else
        demekron!(v, deme)
        true
    end
end

live_condition((; I1, I2), cols) = begin
    ell1, ell2 = ell(cols)
    I1 ≥ ell1 && I2 ≥ ell2
end

singular_root!(
    guidenode, cols, (; I1, I2);
    _...,
) = begin
    ell1, ell2 = ell(cols)
    i, _, p = rcateg(guidenode.present[:,1] .* [I1-ell1, I2-ell2], DemeSet, true)
    if ismissing(i)
        Prob(-Inf), (; I1, I2)
    else
        plant!(cols, i, guidenode.chillins[1])
        -log(p), (; I1, I2)
    end
end

singular_sample!(
    gennode, guidenode, cols, (; I1, I2);
    psi1, psi2, _...,
) = begin
    deme = gennode.deme
    if guidenode.parlin ∉ cols[deme]
        return Prob(-Inf), (; I1, I2)
    end
    chop!(cols, deme, guidenode.parlin)
    ## destructive sampling (r = 1): total sampling rate in the deme
    if deme == Camel
        ll = log(psi1*I1)
        I1 -= 1
    else
        ll = log(psi2*I2)
        I2 -= 1
    end
    ll, (; I1, I2)
end

singular_branch!(
    guidenode, cols, (; I1, I2);
    lambda11, lambda12, lambda21, lambda22, _...,
) = begin
    pr = guidenode.present          # pr[deme, child]
    par = guidenode.parlin
    kids = guidenode.chillins
    if par ∈ cols[Camel]
        a11 = lambda11*I1
        a12 = lambda12*I1
        k, _, p = rcateg([a11*pr[1,1]*pr[1,2],
                          a12*pr[1,1]*pr[2,2],
                          a12*pr[2,1]*pr[1,2]], true)
        if k == 1
            fork!(cols, Camel, par, (Camel, Camel), kids)
            I1 += 1
            ll = log(a11) - log(p) - log(I1*(I1-1)/2)
        elseif k == 2 || k == 3
            fork!(cols, Camel, par, k == 2 ? (Camel, Human) : (Human, Camel), kids)
            I2 += 1
            ll = log(a12) - log(p) - log(I1*I2)
        else
            ll = Prob(-Inf)
        end
    elseif par ∈ cols[Human]
        a22 = lambda22*I2
        a21 = lambda21*I2
        k, _, p = rcateg([a22*pr[2,1]*pr[2,2],
                          a21*pr[1,1]*pr[2,2],
                          a21*pr[2,1]*pr[1,2]], true)
        if k == 1
            fork!(cols, Human, par, (Human, Human), kids)
            I2 += 1
            ll = log(a22) - log(p) - log(I2*(I2-1)/2)
        elseif k == 2 || k == 3
            fork!(cols, Human, par, k == 2 ? (Camel, Human) : (Human, Camel), kids)
            I1 += 1
            ll = log(a21) - log(p) - log(I1*I2)
        else
            ll = Prob(-Inf)
        end
    else
        ll = Prob(-Inf)
    end
    ll, (; I1, I2)
end

singular_part!(
    cols, state, genealogy, guide, n;
    kwargs...,
) = begin
    if !live_condition(state, cols)
        return false, Prob(-Inf), state
    end
    guidenode = guide[n]
    if guidenode.type == Root
        ll, state = singular_root!(guidenode, cols, state; kwargs...)
    elseif guidenode.type == Sample
        ll, state = singular_sample!(genealogy[n], guidenode, cols, state; kwargs...)
    elseif guidenode.type == Node
        ll, state = singular_branch!(guidenode, cols, state; kwargs...)
    else
        @assert false "impossible node type" # COV_EXCL_LINE
    end
    isfinite(ll), ll, state
end

## Unobserved events.  Event types are drawn at the model's own rates;
## guidance enters only through `choose_move` (guide.jl), which decides whether a
## cross-type birth or migration moves a tracked lineage, and which.
event_rates!(
    alpha, (; I1, I2), cols;
    lambda11, lambda12, lambda21, lambda22,
    mu1, mu2, psi1, psi2, m12, m21,
    _...,
) = begin
    ell1, ell2 = ell(cols)
    alpha[1] = lambda11*I1                          # 1 -> 1 birth
    alpha[2] = lambda22*I2                          # 2 -> 2 birth
    alpha[3] = lambda12*I1                          # 1 -> 2 birth
    alpha[4] = lambda21*I2                          # 2 -> 1 birth
    alpha[5] = @indicator(I1 > ell1, mu1*(I1-ell1)) # untracked death, 1
    alpha[6] = @indicator(I2 > ell2, mu2*(I2-ell2)) # untracked death, 2
    alpha[7] = m12*I1                               # migration 1 -> 2
    alpha[8] = m21*I2                               # migration 2 -> 1
    ## decay: tracked hosts must survive, and nothing is sampled off the tree
    psi1*I1 + psi2*I2 + mu1*ell1 + mu2*ell2
end

## A birth from a type-i host that produces a new type-j host.  Target
## factors, from the exchangeable representation:
##   no move: the new host carries no tracked lineage,   1 - ell_j/I_j'
##   move b:  b passes to the new host and the parent is
##            left untracked,                  (1 - (ell_i-1)/I_i)/I_j'
cross_birth!(t, guide, node, cols, i, j, Ii, Ij) = begin
    li, lj = ell(cols,i), ell(cols,j)
    Ij += 1
    f0 = 1 - lj/Ij
    f1 = (1 - (li-1)/Ii)/Ij
    b, q = choose_move(t, guide, node, cols, i, j, f0, f1)
    if q == 0
        ll = Prob(-Inf)
    elseif b == 0
        ll = log(f0) - log(q)
    else
        swap!(cols, i, j, b)
        ll = log(f1) - log(q)
    end
    ll, Ii, Ij
end

## A migration from type i to type j.  The migrant carries any lineage it
## has, so "no move" needs an untracked migrant (impossible if ell_i = I_i):
##   no move: the arriving host carries no tracked lineage, 1 - ell_j/I_j'
##   move b:  the arriving host carries b,                  1/I_j'
migration!(t, guide, node, cols, i, j, Ii, Ij) = begin
    li, lj = ell(cols,i), ell(cols,j)
    f0 = Ii > li ? 1 - lj/(Ij+1) : zero(Prob)
    f1 = 1/(Ij+1)
    b, q = choose_move(t, guide, node, cols, i, j, f0, f1)
    Ii -= 1
    Ij += 1
    if q == 0
        ll = Prob(-Inf)
    elseif b == 0
        ll = log(f0) - log(q)
    else
        swap!(cols, i, j, b)
        ll = log(f1) - log(q)
    end
    ll, Ii, Ij
end

regular_part!(
    cols, state, guide, n, t, tf;
    kwargs...,
) = begin
    alpha = similar(Vector{Prob}, 8)
    ll::Prob = zero(Prob)
    (; I1, I2) = state
    while t < tf
        decay = event_rates!(alpha, (; I1, I2), cols; kwargs...)
        k, s = rcateg(alpha)
        step::Time = -log(rand())/s
        if k > 0 && t+step < tf
            te = t + step       # event time: the guide is evaluated here
            ell1, ell2 = ell(cols)
            if k == 1
                I1 += 1
                ll1 = log(1-ell1*(ell1-1)/I1/(I1-1))
            elseif k == 2
                I2 += 1
                ll1 = log(1-ell2*(ell2-1)/I2/(I2-1))
            elseif k == 3
                ll1, I1, I2 = cross_birth!(te, guide, n, cols, Camel, Human, I1, I2)
            elseif k == 4
                ll1, I2, I1 = cross_birth!(te, guide, n, cols, Human, Camel, I2, I1)
            elseif k == 5
                ll1 = -log(1-ell1/I1)
                I1 -= 1
            elseif k == 6
                ll1 = -log(1-ell2/I2)
                I2 -= 1
            elseif k == 7
                ll1, I1, I2 = migration!(te, guide, n, cols, Camel, Human, I1, I2)
            elseif k == 8
                ll1, I2, I1 = migration!(te, guide, n, cols, Human, Camel, I2, I1)
            else
                @assert false "impossible event" # COV_EXCL_LINE
            end
            ll += ll1 - decay*step
            t = te
        else
            ll -= decay*(tf - t)
            break
        end
    end
    ll, (; I1, I2)
end

"""
    filter_pomp(gen, m; lambda11, lambda12, lambda21, lambda22,
                mu1, mu2, psi1, psi2, m12, m21, I1_0 = 1, I2_0 = 0)

Construct the guided MTBD filter for genealogy `gen`, which should be
parsed with `demes = GuidedMTBD.Demes`, using the guiding finite-state
Markov process `m` (construct it with [`fsmarkov`](@ref)).  Rates
default to the MERS maximum-likelihood estimate ([`mers_mle`](@ref));
see [`epi_params`](@ref) to build them from reproduction numbers.
"""
filter_pomp(
    gen::Genealogy,
    m::FSMarkovProc;
    lambda11 = mers_mle.lambda11, lambda12 = mers_mle.lambda12,
    lambda21 = mers_mle.lambda21, lambda22 = mers_mle.lambda22,
    mu1 = mers_mle.mu1, mu2 = mers_mle.mu2,
    psi1 = mers_mle.psi1, psi2 = mers_mle.psi2,
    m12 = mers_mle.m12, m21 = mers_mle.m21,
    I1_0 = 1, I2_0 = 0,
) = begin
    check(gen)
    guidegen = guide(gen, m, knowledge!)
    pomp(
        params = (
            lambda11 = Float64(lambda11), lambda12 = Float64(lambda12),
            lambda21 = Float64(lambda21), lambda22 = Float64(lambda22),
            mu1 = Float64(mu1), mu2 = Float64(mu2),
            psi1 = Float64(psi1), psi2 = Float64(psi2),
            m12 = Float64(m12), m21 = Float64(m21),
            I1_0 = Float64(I1_0), I2_0 = Float64(I2_0),
        ),
        t0 = timezero(guidegen),
        times = times(guidegen),
        userdata = (genealogy = gen, guide = guidegen),
        init_state = (
            node = one(Name), ll = zero(Prob), live = true,
            cols = Coloring(Demes), state = (I1 = zero(Int64), I2 = zero(Int64)),
        ),
        rinit = function (; kwargs...)
            (node = one(Name), ll = zero(Prob), live = true,
             cols = Coloring(Demes), state = mtbd_rinit(; kwargs...))
        end,
        rprocess = onestep(
            function (; t, dt, genealogy, guide, node, ll, live, cols, state, kwargs...)
                if live
                    cols = copy(cols)
                    live, ll, state = singular_part!(cols, state, genealogy, guide, node; kwargs...)
                else
                    ll = Prob(-Inf)
                end
                if live && dt > 0
                    ll1, state = regular_part!(cols, state, guide, node, t, t+dt; kwargs...)
                    ll += ll1
                end
                (; node = node + 1, ll, live, cols, state)
            end,
        ),
        logdmeasure = function (; live, ll, _...)
            live ? ll : Prob(-Inf)
        end,
    )
end

end
