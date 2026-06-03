module NaiveSEIR

using ..PhyloPOMP
import PartiallyObservedMarkovProcesses as POMP

@demes SEIR Expos Infec
using .SEIR: Expos, Infec, T as DemeType

seir_singular(
    geneal, node, ll, cols;
    S, E, I, R, pop, β, ψ,
    _...,
) = begin
    n = geneal[node]
    (ellE,ellI) = ell(cols)
    @assert I ≥ ellI && E ≥ ellE
    if n.type==PhyloPOMP.Root
        if length(n.children) == 1
            if E-ellE+I-ellI > 0
                i, _, p = rcateg([E-ellE, I-ellI], DemeType, true)
                ll -= log(p)
                (ellE,ellI) = plant!(cols,i,n.lineage)
            else
                ll += Float64(-Inf)
                (ellE,ellI) = plant!(cols,Infec,n.lineage)
                I += 1
            end
        else
            error("too many children ($(length(n.children)) > 1) at root $(n.name), t=$(n.time)")
        end
    elseif n.type==PhyloPOMP.Sample
        if n.lineage ∉ cols[Infec]
            ll += Float64(-Inf)
            (ellE,ellI) = swap!(cols,Expos,Infec,n.lineage)
            E -= 1
            I += 1
        end
        ll += log(ψ*I)
        if length(n.children) == 0
            (ellE,ellI) = chop!(cols,Infec,n.lineage)
            ll += log(1-ell(cols,Infec)/I)
        elseif length(n.children) == 1
            (ellE,ellI) = chop!(cols,Infec,n.lineage,Infec,geneal[n.children[1]].lineage)
            ll -= log(I)
        else
            error("too many children ($(length(n.children)) > 1) at sample $(n.name), t=$(n.time)")
        end
    elseif n.type==PhyloPOMP.Node
        if n.lineage ∉ cols[Infec]
            ll += Float64(-Inf)
            (ellE,ellI) = swap!(cols,Expos,Infec,n.lineage)
            E -= 1
            I += 1
        end
        if length(n.children) == 2
            ll += log(β*S*I/pop)
            k, _, p = rcateg([1, 1], true)
            ll -= log(p)
            if k==1
                (ellE,ellI) = fork!(
                    cols,
                    Infec,Expos,Infec,
                    n.lineage,
                    geneal[n.children[1]].lineage,
                    geneal[n.children[2]].lineage,
                )
            else
                (ellE,ellI) = fork!(
                    cols,
                    Infec,Expos,Infec,
                    n.lineage,
                    geneal[n.children[2]].lineage,
                    geneal[n.children[1]].lineage,
                )
            end
            S -= 1
            E += 1
            ll -= log(E*I)
        else
            error("too many children ($(length(n.children)) ≠ 2) at node $(n.name), t=$(n.time)")
        end
    else
        @assert false "impossible node type"
    end
    (; ll = ll, cols = cols, S = S, E = E, I = I, R = R)
end

event_rates!(
    alpha, pi, cols;
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
    pi[1] = @indicator(I > 0, 1-ellI/I)
    pi[2] = @indicator(I > 0, ellI/I)
    pi[3] = @indicator(E > 0, 1-ellE/E)
    pi[4] = @indicator(E > 0, ellE/E)
    pi[6] = pi[5] = 1.0
    ψ*I + @indicator(I ≤ ellI, γ*I)
end

seir_regular(
    ll, cols;
    t, dt,
    S, E, I, R,
    β, σ, γ, ω, ψ, pop,
    _...,
) = begin
    alpha = similar(Vector{Float64}, 6)
    pi = similar(Vector{Float64}, 6)
    (ellE,ellI) = ell(cols)
    @assert I ≥ ellI && E ≥ ellE
    tf = t+dt
    if t < tf
        penalty = event_rates!(
            alpha, pi, cols;
            S = S, E = E, I = I, R = R,
            β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, pop = pop,
        )
        k, s = rcateg(alpha .* pi)
        step::PhyloPOMP.Time = -log(rand())/s
        while (t+step < tf)
            ll -= penalty*step+log(pi[k])
            if k==1
                S -= 1
                E += 1
                ll += log(1-ellE/E)
            elseif k==2
                ll += log(ellI)
                b = rand(cols[Infec])
                (ellE,ellI) = swap!(cols,Infec,Expos,b)
                S -= 1
                E += 1
                ll += log(1-ellI/I)-log(E)
            elseif k==3
                E -= 1
                I += 1
                ll += log(1-ellI/I)
            elseif k==4
                ll += log(ellE)
                b = rand(cols[Expos])
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
            penalty = event_rates!(
                alpha, pi, cols;
                S = S, E = E, I = I, R = R,
                β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, pop = pop,
            )
            k, s = rcateg(alpha .* pi)
            step = -log(rand())/s
        end
        step = tf - t
        ll -= penalty*step
    end
    @assert I ≥ ellI && E ≥ ellE
    (; ll = ll, cols = cols, S = S, E = E, I = I, R = R)
end

seir(
    gen::Genealogy;
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
        t0 = timezero(gen),
        times = times(gen),
        rinit = function (; S0, E0, I0, R0, pop, _...)
            m = pop/(S0+E0+I0+R0)
            (
                node = one(PhyloPOMP.Name),
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
                ; node, ll, cols, geneal,
                S, E, I, R,
                args...,
                )
                cols = copy(cols)
                ll = zero(Float64)
                ll, cols, S, E, I, R = seir_singular(
                    geneal, node, ll, cols;
                    S = S, E = E, I = I, R = R,
                    args...,
                )
                ll, cols, S, E, I, R = seir_regular(
                    ll, cols;
                    S = S, E = E, I = I, R = R,
                    args...,
                )
                (; node = node+1, ll = ll, cols = cols,
                 S = S, E = E, I = I, R = R)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,),
    )
end

end
