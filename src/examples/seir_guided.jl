module GuidedSEIR

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Time, Prob
import PartiallyObservedMarkovProcesses as POMP

@demes SEIR Expos Infec
using .SEIR: Expos, Infec, T as DemeType

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
                i, _, p = rcateg(n.present[:,1].*[E-ellE, I-ellI], DemeType, true)
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

choice1(g, x, cols, rh) = begin
    s = sum(rh[g.linmap[i]] for i ∈ cols; init = zero(Prob))
    v = [x-ell(cols), s]
    v./sum(v)
end

choice2(g, cols, rh) = begin
    lins = map(i->g.linmap[i],collect(cols))
    k,_,p = rcateg(rh[lins],true)
    g.alllins[lins[k]], p
end

event_rates!(
    alpha, pi, cols, node, rh;
    S, E, I, R,
    β, σ, γ, ω, ψ, pop,
    _...,
) = begin
    (ellE,ellI) = ell(cols)
    @assert I ≥ ellI && E ≥ ellE
    alpha[2] = alpha[1] = β*S*I/pop
    alpha[4] = alpha[3] = σ*E
    alpha[5] = @indicator(I > ellI, γ*I)
    alpha[6] = ω*R
    cI = choice1(node,I,cols[Infec],rh[(Infec,Expos)])
    pi[1] = @indicator(I > 0, cI[1])
    pi[2] = @indicator(I > 0, cI[2])
    cE = choice1(node,E,cols[Expos],rh[(Expos,Infec)])
    pi[3] = @indicator(E > 0, cE[1])
    pi[4] = @indicator(E > 0, cE[2])
    pi[6] = pi[5] = 1.0
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
        alpha = similar(Vector{Float64}, 6)
        pi = similar(Vector{Float64}, 6)
        (ellE,ellI) = ell(cols)
        @assert I ≥ ellI && E ≥ ellE
        rh = relhaz(guide,t,node)
        decay = event_rates!(
            alpha, pi, cols, n, rh;
            S = S, E = E, I = I, R = R,
            β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, pop = pop,
        )
        k, s = rcateg(alpha .* pi)
        step::Time = -log(rand())/s
        while (t+step < tf)
            ll -= decay*step+log(pi[k])
            if k==1
                S -= 1
                E += 1
                ll += log(1-ellE/E)
            elseif k==2
                b,p = choice2(n,cols[Infec],rh[(Infec,Expos)])
                ll -= log(p)
                (ellE,ellI) = swap!(cols,Infec,Expos,b)
                S -= 1
                E += 1
                ll += log(1-ellI/I)-log(E)
            elseif k==3
                E -= 1
                I += 1
                ll += log(1-ellI/I)
            elseif k==4
                b,p = choice2(n,cols[Expos],rh[(Expos,Infec)])
                ll -= log(p)
                (ellE,ellI) = swap!(cols,Expos,Infec,b)
                E -= 1
                I += 1
                ll -= log(I)
            elseif k==5
                I -= 1
                R += 1
            elseif k==6
                R -= 1
                S += 1
            else
                @assert false "impossible error!"
            end
            @assert I ≥ ellI && E ≥ ellE
            t += step
            decay = event_rates!(
                alpha, pi, cols, n, rh;
                S = S, E = E, I = I, R = R,
                β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, pop = pop,
            )
            k, s = rcateg(alpha .* pi)
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
