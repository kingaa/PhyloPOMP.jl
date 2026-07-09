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

singular_part!(
    cols, guide, node, ll,
    S, E, I, R;
    pop, β, ψ, χ,
    _...,
) = begin
    ellE, ellI = ell(cols)
    n = guide[node]
    @assert I ≥ ellI && E ≥ ellE
    if n.type==Root
        if length(n.chillins) == 1
            if E-ellE+I-ellI > 0
                i, _, p = rcateg(n.present[:,1].*[E-ellE, I-ellI], DemeSet, true)
                ll -= log(p)
                ellE, ellI = plant!(cols,i,n.chillins[1])
            else
                ## even though this realization is incompatible with the data,
                ## it is necessary to correct the coloring to avoid downstream errors.
                ll += Prob(-Inf)
                ellE, ellI = plant!(cols,Infec,n.chillins[1])
                I += 1
            end
        else
            error("too many children ($(length(n.chillins)) > 1) at root $(n.name), t=$(n.time)")
        end
    elseif n.type==Sample
        if n.parlin ∉ cols[Infec]
            ## even though this realization is incompatible with the data,
            ## it is necessary to correct the coloring to avoid downstream errors.
            ll += Prob(-Inf)
            ellE, ellI = swap!(cols,Expos,Infec,n.parlin)
            E -= 1
            I += 1
        end
        if length(n.chillins) == 0
            ll += log((ψ+χ)*I)
            ellE, ellI = chop!(cols,Infec,n.parlin)
            k,_ = rcateg([ψ, χ])
            if k==1             # non-destructive sampling
                ll += log(1-ellI/I)
            elseif k==2         # destructive sampling
                I -= 1
            else
                @assert false "impossible choice" # COV_EXCL_LINE
            end
        elseif length(n.chillins) == 1
            ellE, ellI = chop!(cols,Infec,n.parlin,Infec,n.chillins[1])
            ll += log(ψ)
        else
            error("too many children ($(length(n.chillins)) > 1) at sample $(n.name), t=$(n.time)")
        end
    elseif n.type==Node
        if n.parlin ∉ cols[Infec]
            ## even though this realization is incompatible with the data,
            ## it is necessary to correct the coloring to avoid downstream errors.
            ll += Prob(-Inf)
            ellE, ellI = swap!(cols,Expos,Infec,n.parlin)
            E -= 1
            I += 1
        end
        if length(n.chillins) == 2
            ll += log(β*S*I/pop)
            k, _, p = rcateg([n.present[1,1]*n.present[2,2], n.present[1,2]*n.present[2,1]], true)
            ll -= log(p)
            if k==1
                ellE, ellI = fork!(cols,Infec,n.parlin,(Expos,Infec),n.chillins)
            else
                ellE, ellI = fork!(cols,Infec,n.parlin,(Infec,Expos),n.chillins)
            end
            S -= 1
            E += 1
            ll -= log(E*I)
        else
            error("wrong number of children ($(length(n.chillins)) ≠ 2) at node $(n.name), t=$(n.time)")
        end
    else
        @assert false "impossible node type" # COV_EXCL_LINE
    end
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
                kwargs...,
                )
                t = guide[node].tbeg
                tf = guide[node].tend
                cols = copy(cols)
                ll = zero(Prob)
                ll, S, E, I, R = singular_part!(
                    cols, guide, node, ll,
                    S, E, I, R;
                    kwargs...,
                )
                if t < tf && isfinite(ll)
                    ll, S, E, I, R = regular_part!(
                        cols, guide, node, ll,
                        t, tf, S, E, I, R;
                        kwargs...,
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
