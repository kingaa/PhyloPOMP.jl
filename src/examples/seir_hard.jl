"""
    HardSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using a filter guide and "hard" proposals.  That is, color
changes can be proposed on branches at rates that exceed the event rates
in the underlying population process.
"""
module HardSEIR

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time, FSMarkovProc

@demes Demes Expos Infec
using .Demes: Expos, Infec, DemeSet

include("seir_trees.jl")
include("seir_funs.jl")

event_rates!(
    alpha, preboost, t, guide, node,
    cols, ellE, ellI;
    S, E, I, R,
    β, σ, γ, ω, ψ, χ, pop,
    _...,
) = begin
    rhEI = relhaz(t,guide,node,Expos,Infec,cols)
    rhIE = relhaz(t,guide,node,Infec,Expos,cols)
    a = β*S*I/pop
    preboost[1] = I > 0 ? 1-ellI/I : 0
    alpha[1] = a*preboost[1]
    preboost[2] = I > 0 ? sum(rhIE)/I : 1
    alpha[2] = a*preboost[2]
    b = σ*E
    preboost[3] = E > 0 ? 1-ellE/E : 0
    alpha[3] = b*preboost[3]
    preboost[4] = E > 0 ? sum(rhEI)/E : 1
    alpha[4] = b*preboost[4]
    c = γ*I
    preboost[5] = 1-ellI/I
    alpha[5] = @indicator(I > ellI, c*preboost[5])
    preboost[6] = 1.0
    alpha[6] = ω*R
    decay = ψ*I + χ*I +
        a - alpha[1] - alpha[2] +
        b - alpha[3] - alpha[4] +
        c - alpha[5]
    decay, rhEI, rhIE
end

regular_part!(
    cols, guide, node, ll;
    S, E, I, R,
    β, σ, γ, ω, ψ, χ, pop,
    _...,
) = begin
    n = guide[node]
    t = n.tbeg
    tf = n.tend
    if t < tf
        alpha = similar(Vector{Prob},6)
        preboost = similar(Vector{Prob},6)
        ellE, ellI = ell(cols)
        decay,rhEI,rhIE = event_rates!(
            alpha, preboost, t, guide, node,
            cols, ellE, ellI;
            S = S, E = E, I = I, R = R,
            β = β, σ = σ, γ = γ, ω = ω,
            ψ = ψ, χ = χ, pop = pop,
        )
        k, s = rcateg(alpha)
        step::Time = -log(rand())/s
        while (t+step < tf)
            ll -= decay*step + log(preboost[k])
            if k==1
                S -= 1
                E += 1
                ll += log(1-ellE/E)
            elseif k==2
                b,_,p = rcateg(rhIE,cols[Infec],true)
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
                b,_,p = rcateg(rhEI,cols[Expos],true)
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
            else
                @assert false "impossible error!" # COV_EXCL_LINE
            end
            t += step
            decay,rhEI,rhIE = event_rates!(
                alpha, preboost, t, guide, node,
                cols, ellE, ellI;
                S = S, E = E, I = I, R = R,
                β = β, σ = σ, γ = γ, ω = ω,
                ψ = ψ, χ = χ, pop = pop,
            )
            k, s = rcateg(alpha)
            step = -log(rand())/s
        end
        step = tf - t
        ll -= decay*step
    end
    @assert I ≥ ellI && E ≥ ellE
    ll, S, E, I, R
end

"""
    filter_pomp(g; β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
         pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08)

Constructs a pomp object based on the filter guide `g`.
"""
filter_pomp(
    gen::Genealogy,
    m::FSMarkovProc;
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
    pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    guidegen = guide(gen,m,knowledge!)
    pomp(
        params = (
            β = Float64(β), σ = Float64(σ), γ = Float64(γ),
            ω = Float64(ω), ψ = Float64(ψ), χ = Float64(χ),
            pop = Float64(pop),
            S0 = Float64(S0), E0 = Float64(E0),
            I0 = Float64(I0), R0 = Float64(R0),
        ),
        t0 = timezero(guidegen),
        times = times(guidegen),
        rinit = function (; S0, E0, I0, R0, pop, _...)
            m = pop/(S0+E0+I0+R0)
            (
                node = one(Name),
                ll = zero(Prob),
                cols = Coloring(Demes),
                S = round(Int64, m*Float64(S0)),
                E = round(Int64, m*Float64(E0)),
                I = round(Int64, m*Float64(I0)),
                R = round(Int64, m*Float64(R0)),
            )
        end,
        rprocess = onestep(
            function (
                ; node, ll, cols, guide,
                S, E, I, R,
                args...,
                )
                cols = copy(cols)
                ll = zero(Prob)
                ll, S, E, I, R = singular_part!(
                    cols, guide, node, ll;
                    S = S, E = E, I = I, R = R,
                    args...,
                )
                if isfinite(ll)
                    ll, S, E, I, R = regular_part!(
                        cols, guide, node, ll;
                        S = S, E = E, I = I, R = R,
                        args...,
                    )
                end
                (; node = node+1, ll = ll, cols = cols,
                 S = S, E = E, I = I, R = R)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (guide = guidegen,),
    )
end

end
