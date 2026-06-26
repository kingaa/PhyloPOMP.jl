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
    alpha, ellI;
    S, E, I, R,
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
    cols, guide, node, ll;
    S, E, I, R,
    β, σ, γ, ω, ψ, χ, pop,
    _...,
) = begin
    n = guide[node]
    t = n.tbeg
    tf = n.tend
    if t < tf
        alpha = similar(Vector{Prob},4)
        ellE, ellI = ell(cols)
        @assert I ≥ ellI && E ≥ ellE
        decay = event_rates!(
            alpha, ellI;
            S = S, E = E, I = I, R = R,
            β = β, σ = σ, γ = γ, ω = ω,
            ψ = ψ, χ = χ, pop = pop,
        )
        k, s = rcateg(alpha)
        step::Time = -log(rand())/s
        while (t+step < tf)
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
            @assert I ≥ ellI && E ≥ ellE
            t += step
            decay = event_rates!(
                alpha, ellI;
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
