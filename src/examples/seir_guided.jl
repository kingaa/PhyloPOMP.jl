module GuidedSEIR

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Time, Prob
import PartiallyObservedMarkovProcesses as POMP

@demes SEIR Expos Infec
using .SEIR: Expos, Infec, DemeSet

seir_singular!(
    cols, guide, node, ll;
    S, E, I, R, pop, β, ψ,
    _...,
) = begin
    (ellE,ellI) = ell(cols)
    @assert I ≥ ellI && E ≥ ellE
    n = guide[node]
    if n.type==Root
        if length(n.chillins) == 1
            if E-ellE+I-ellI > 0
                i, _, p = rcateg(n.present[:,1].*[E-ellE, I-ellI], DemeSet, true)
                ll -= log(p)
                (ellE,ellI) = plant!(cols,i,n.chillins[1])
            else
                ll += Float64(-Inf)
                (ellE,ellI) = plant!(cols,Infec,n.chillins[1])
                I += 1
            end
        else
            error("too many children ($(length(n.chillins)) > 1) at root $(n.name), t=$(n.time)")
        end
    elseif n.type==Sample
        if n.parlin ∉ cols[Infec]
            ll += Float64(-Inf)
            (ellE,ellI) = swap!(cols,Expos,Infec,n.parlin)
            E -= 1
            I += 1
        end
        ll += log(ψ*I)
        if length(n.chillins) == 0
            (ellE,ellI) = chop!(cols,Infec,n.parlin)
            ll += log(1-ellI/I)
        elseif length(n.chillins) == 1
            (ellE,ellI) = chop!(cols,Infec,n.parlin,Infec,n.chillins[1])
            ll -= log(I)
        else
            error("too many children ($(length(n.chillins)) > 1) at sample $(n.name), t=$(n.time)")
        end
    elseif n.type==Node
        if n.parlin ∉ cols[Infec]
            ll += Float64(-Inf)
            (ellE,ellI) = swap!(cols,Expos,Infec,n.parlin)
            E -= 1
            I += 1
        end
        if length(n.chillins) == 2
            ll += log(β*S*I/pop)
            k, _, p = rcateg([n.present[1,1]*n.present[2,2], n.present[1,2]*n.present[2,1]], true)
            ll -= log(p)
            if k==1
                (ellE,ellI) = fork!(cols,Infec,n.parlin,(Expos,Infec),n.chillins)
            else
                (ellE,ellI) = fork!(cols,Infec,n.parlin,(Infec,Expos),n.chillins)
            end
            S -= 1
            E += 1
            ll -= log(E*I)
        else
            error("too many children ($(length(n.chillins)) ≠ 2) at node $(n.name), t=$(n.time)")
        end
    else
        @assert false "impossible node type"
    end
    (; ll = ll, S = S, E = E, I = I, R = R)
end

event_rates!(
    alpha, ellI;
    S, E, I, R,
    β, σ, γ, ω, ψ, pop,
    _...,
) = begin
    alpha[1] = β*S*I/pop
    alpha[2] = σ*E
    alpha[3] = @indicator(I > ellI, γ*I)
    alpha[4] = ω*R
    ψ*I + @indicator(I ≤ ellI, γ*I)
end

seir_regular!(
    cols, guide, node, ll;
    S, E, I, R,
    β, σ, γ, ω, ψ, pop,
    _...,
) = begin
    n = guide[node]
    t = n.tbeg
    tf = n.tend
    if t < tf
        alpha = similar(Vector{Float64},4)
        (ellE,ellI) = ell(cols)
        @assert I ≥ ellI && E ≥ ellE
        decay = event_rates!(
            alpha, ellI;
            S = S, E = E, I = I, R = R,
            β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, pop = pop,
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
                    (ellE,ellI) = swap!(cols,Infec,Expos,b)
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
                    (ellE,ellI) = swap!(cols,Expos,Infec,b)
                    ll -= log(I)
                end
            elseif k==3
                I -= 1
                R += 1
            elseif k==4
                R -= 1
                S += 1
            else
                @assert false "impossible error!"
            end
            @assert I ≥ ellI && E ≥ ellE
            t += step
            decay = event_rates!(
                alpha, ellI;
                S = S, E = E, I = I, R = R,
                β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, pop = pop,
            )
            k, s = rcateg(alpha)
            step = -log(rand())/s
        end
        step = tf - t
        ll -= decay*step
    end
    @assert I ≥ ellI && E ≥ ellE
    (; ll = ll, S = S, E = E, I = I, R = R)
end

seir(
    guide::Guide;
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02,
    pop = 100,
    S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    pomp(
        params = (
            β = Float64(β), σ = Float64(σ), γ = Float64(γ),
            ω = Float64(ω), ψ = Float64(ψ),
            pop = Float64(pop),
            S0 = Float64(S0), E0 = Float64(E0),
            I0 = Float64(I0), R0 = Float64(R0),
        ),
        t0 = timezero(guide),
        times = times(guide),
        rinit = function (; S0, E0, I0, R0, pop, _...)
            m = pop/(S0+E0+I0+R0)
            (
                node = one(Name),
                ll = zero(Float64),
                cols = Coloring(SEIR),
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
                ll = zero(Float64)
                ll, S, E, I, R = seir_singular!(
                    cols, guide, node, ll;
                    S = S, E = E, I = I, R = R,
                    args...,
                )
                if isfinite(ll)
                    ll, S, E, I, R = seir_regular!(
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
        userdata = (guide = guide,),
    )
end

end
