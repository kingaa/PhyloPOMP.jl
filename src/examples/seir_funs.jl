## This file contains definitions that are used in the guided filters
## for the SEIR model.

"""
    knowledge!(v; deme, type, time)

This function should return true if the guide probabilities will be fixed
at this node and false otherwise. If the former, it should fill the vector
`v` with an appropriate probability vector.
"""
knowledge!(
    v; deme, type, time,
) = begin
    if type==Sample || type==Node
        demekron!(v,Infec)
        true
    else
        false
    end
end

check(
    gen::Genealogy,
) = begin
    for node ∈ eachindex(gen)
        n = gen[node]
        if n.type == Root
            @assert length(n.children)==1 "wrong number of children ($(length(n.children)) != 1) at root $(n.name), t=$(n.slate)"
        elseif n.type == Sample
            @assert length(n.children)<2 "too many children ($(length(n.children)) > 1) at sample $(n.name), t=$(n.time)"
        elseif n.type == Node
            @assert length(n.children)==2 "wrong number of children ($(length(n.children)) ≠ 2) at node $(n.name), t=$(n.time)"
        end
    end
    nothing
end

transmission!(
    alpha, pi;
    β, S, I, pop,
    onI, offI = I - onI,
    _...,
) = begin
    rate = β*S*I/pop
    pi[1] = @indicator(I > 0, offI/I)
    pi[2] = @indicator(I > 0, onI/I)
    alpha[1] = rate*pi[1]
    alpha[2] = rate*pi[2]
    rate-sum(alpha)
end

progression!(
    alpha, pi;
    σ, E,
    onE, offE = E - onE,
    _...,
) = begin
    rate = σ*E
    pi[1] = @indicator(E > 0, offE/E)
    pi[2] = @indicator(E > 0, onE/E)
    alpha[1] = rate*pi[1]
    alpha[2] = rate*pi[2]
    rate-sum(alpha)
end

recovery!(
    alpha, pi;
    γ, I,
    ellI,
    _...,
) = begin
    rate = γ*I
    pi[1] = @indicator(I > ellI, 1-ellI/I)
    alpha[1] = rate*pi[1]
    rate-sum(alpha)
end

waning!(
    alpha, pi;
    ω, R,
    _...,
) = begin
    rate = ω*R
    pi[1] = 1
    alpha[1] = rate
    zero(Float64)
end

sampling(
    ;ψ, χ, I,
    _...,
) = begin
    (ψ+χ)*I
end

event_rates!(alpha, pi; kwargs...) = begin
    decay = zero(Prob)
    decay += transmission!(
        @view(alpha[1:2]),@view(pi[1:2]);
        kwargs...,
    )
    decay += progression!(
        @view(alpha[3:4]),@view(pi[3:4]);
        kwargs...,
    )
    decay += recovery!(
        @view(alpha[5]),@view(pi[5]);
        kwargs...,
    )
    decay += waning!(
        @view(alpha[6]),@view(pi[6]);
        kwargs...,
    )
    decay += sampling(;kwargs...,)
    decay
end

singular_part!(
    cols, guide, node, ll, live,
    S, E, I, R;
    pop, β, ψ, χ,
    _...,
) = begin
    ellE, ellI = ell(cols)
    if I < ellI || E < ellE
        live = false
    end
    if live
        n = guide[node]
        if n.type==Root
            i, _, p = rcateg(n.present[:,1].*[E-ellE, I-ellI], DemeSet, true)
            ll -= log(p)
            if ismissing(i)
                live = false
            else
                ellE, ellI = plant!(cols,i,n.chillins[1])
            end
        elseif n.type==Sample
            if n.parlin ∉ cols[Infec]
                live = false
            elseif length(n.chillins) == 0
                k,_,p = rcateg([ψ, χ],true)
                ll -= log(p)
                ellE, ellI = chop!(cols,Infec,n.parlin)
                if k==0
                    live = false
                elseif k==1         # non-destructive sample
                    ll += log(ψ*(I-ellI));
                elseif k==2         # destructive sample
                    ll += log(χ*I)
                    I -= 1
                end
            elseif length(n.chillins) == 1
                ellE, ellI = chop!(cols,Infec,n.parlin,Infec,n.chillins[1])
                ll += log(ψ)
            end
        elseif n.type==Node
            if n.parlin ∉ cols[Infec]
                live = false
            else
                ll += log(β*S*I/pop)
                k, _, p = rcateg([n.present[1,1]*n.present[2,2], n.present[1,2]*n.present[2,1]], true)
                ll -= log(p)
                @assert k ≠ 0
                if k==1
                    ellE, ellI = fork!(cols,Infec,n.parlin,(Expos,Infec),n.chillins)
                else
                    ellE, ellI = fork!(cols,Infec,n.parlin,(Infec,Expos),n.chillins)
                end
                if S > 0
                    S -= 1
                end
                E += 1
                ll -= log(E*I)
            end
        else
            @assert false "impossible node type" # COV_EXCL_LINE
        end
    end
    if !live
        ll = Prob(-Inf)
    end
    ll, S, E, I, R, live
end

"""
    filter_pomp(g, m; β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02,
         χ = 0.0, pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08)

Constructs a pomp object based on the genealogy `g` and finite-state
Markov guiding process `m`.
"""
filter_pomp(
    gen::Genealogy,
    m::FSMarkovProc;
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
    pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    check(gen)
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
                live = true
            )
        end,
        rprocess = onestep(
            function (
                ; t, dt,
                node, ll, cols, guide,
                S, E, I, R, live,
                kwargs...,
                )
                tf = t+dt
                cols = copy(cols)
                ll = zero(Prob)
                ll, S, E, I, R, live = singular_part!(
                    cols, guide, node, ll, live,
                    S, E, I, R;
                    kwargs...,
                )
                if live && t < tf && isfinite(ll)
                    ll, S, E, I, R = regular_part!(
                        cols, guide, node, ll,
                        t, tf, S, E, I, R;
                        kwargs...,
                    )
                end
                (;node = node+1, ll, cols, S, E, I, R, live)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (guide = guidegen,),
    )
end
