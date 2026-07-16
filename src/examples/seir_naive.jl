"""
    NaiveSEIR

A module containing an implementation of the phylodynamic filter for
an SEIR model, using so-called naïve proposals. These proposals are
non-anticipatory. In particular, population events occurring at any
time result in color changes with a probability proportional to the
fraction of appropriate lineages represented in the genealogy at that
time.
"""
module NaiveSEIR

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time

@demes Demes Expos Infec
using .Demes: Expos, Infec, DemeSet

include("seir_trees.jl")

singular_part!(
    cols, geneal, node, ll,
    S, E, I, R;
    pop, β, ψ, χ,
    _...,
) = begin
    ellE, ellI = ell(cols)
    n = geneal[node]
    @assert I ≥ ellI && E ≥ ellE
    if n.type==Root
        if length(n.children) == 1
            if E-ellE+I-ellI > 0
                i, _, p = rcateg([E-ellE, I-ellI], DemeSet, true)
                ll -= log(p)
                ellE, ellI = plant!(cols,i,n.lineage)
            else
                ## even though this realization is incompatible with the data,
                ## it is necessary to correct the coloring to avoid downstream errors.
                ll += Prob(-Inf)
                ellE, ellI = plant!(cols,Infec,n.lineage)
                I += 1
            end
        else
            error("too many children ($(length(n.children)) > 1) at root $(n.name), t=$(n.time)")
        end
    elseif n.type==Sample
        if n.lineage ∉ cols[Infec]
            ## even though this realization is incompatible with the data,
            ## it is necessary to correct the coloring to avoid downstream errors.
            ll += Prob(-Inf)
            ellE, ellI = swap!(cols,Expos,Infec,n.lineage)
            E -= 1
            I += 1
        end
        if length(n.children) == 0
            k,_,p = rcateg([ψ, χ],true)
            ll -= log(p)
            ellE, ellI = chop!(cols,Infec,n.lineage)
            if k==1             # non-destructive sample
                ll += log(ψ*(I-ellI));
            elseif k==2         # destructive sample
                ll += log(χ*I)
                I -= 1
            end
        elseif length(n.children) == 1
            chillin = geneal[n.children[1]].lineage
            ellE, ellI = chop!(cols,Infec,n.lineage,Infec,chillin)
            ll += log(ψ)
        else
            error("too many children ($(length(n.children)) > 1) at sample $(n.name), t=$(n.time)")
        end
    elseif n.type==Node
        if n.lineage ∉ cols[Infec]
            ## even though this realization is incompatible with the data,
            ## it is necessary to correct the coloring to avoid downstream errors.
            ll += Prob(-Inf)
            ellE, ellI = swap!(cols,Expos,Infec,n.lineage)
            E -= 1
            I += 1
        end
        if length(n.children) == 2
            chillins = map(n.children) do i
                geneal[i].lineage
            end
            ll += log(β*S*I/pop)
            k, _, p = rcateg([1, 1], true)
            ll -= log(p)
            if k==1
                ellE, ellI = fork!(cols,Infec,n.lineage,(Expos,Infec),chillins)
            else
                ellE, ellI = fork!(cols,Infec,n.lineage,(Infec,Expos),chillins)
            end
            S -= 1
            E += 1
            ll -= log(E*I)
        else
            error("too many children ($(length(n.children)) ≠ 2) at node $(n.name), t=$(n.time)")
        end
    else
        @assert false "impossible node type" # COV_EXCL_LINE
    end
    ll, S, E, I, R
end

event_rates!(
    alpha, pi, cols,
    S, E, I, R;
    β, σ, γ, ω, ψ, χ, pop,
    _...,
) = begin
    ellE, ellI = ell(cols)
    @assert I ≥ ellI && E ≥ ellE
    alpha[2] = alpha[1] = β*S*I/pop
    alpha[4] = alpha[3] = σ*E
    alpha[5] = @indicator(I > ellI, γ*(I-ellI))
    alpha[6] = ω*R
    pi[1] = @indicator(I > 0, 1-ellI/I)
    pi[2] = @indicator(I > 0, ellI/I)
    pi[3] = @indicator(E > 0, 1-ellE/E)
    pi[4] = @indicator(E > 0, ellE/E)
    pi[6] = pi[5] = 1.0
    ψ*I + χ*I + γ*ellI + @indicator(I ≤ ellI, γ*(I-ellI))
end

regular_part!(
    cols, ll,
    t, dt,
    S, E, I, R;
    kwargs...,
) = begin
    tf = t+dt
    if t < tf
        alpha = similar(Vector{Prob}, 6)
        pi = similar(Vector{Prob}, 6)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellE, ellI = ell(cols)
        while t < tf
            decay = event_rates!(
                alpha, pi, cols,
                S, E, I, R;
                kwargs...,
            )
            k, s = rcateg(alpha .* pi)
            step = -log(rand())/s
            if t+step < tf
                ll -= decay*step+log(pi[k])
                if k==1
                    S -= 1
                    E += 1
                    ll += log(1-ellE/E)
                elseif k==2
                    ll += log(ellI)
                    b = rand(cols[Infec])
                    ellE, ellI = swap!(cols,Infec,Expos,b)
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
                    ellE, ellI = swap!(cols,Expos,Infec,b)
                    E -= 1
                    I += 1
                    ll -= log(I)
                elseif k==5
                    ll -= log(1-ellI/I)
                    I -= 1
                    R += 1
                elseif k==6
                    R -= 1
                    S += 1
                end
                t += step
            else
                step = tf - t
                ll -= decay*step
                break
            end
        end
        @assert I ≥ ellI && E ≥ ellE
    end
    ll, S, E, I, R
end

"""
    filter_pomp(gen; β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
         pop = 100, S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08)

Constructs a pomp object based on the genealogy `gen` and implementing
the naïve genealogical filter.
"""
filter_pomp(
    gen::Genealogy;
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
    pop = 100,
    S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    pomp(
        params = (
            β = Float64(β), σ = Float64(σ), γ = Float64(γ),
            ω = Float64(ω), ψ = Float64(ψ), χ = Float64(χ),
            pop = Float64(pop),
            S0 = Float64(S0), E0 = Float64(E0),
            I0 = Float64(I0), R0 = Float64(R0),
        ),
        t0 = timezero(gen),
        times = times(gen),
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
                ; node, ll, cols, geneal,
                t, dt,
                S, E, I, R,
                args...,
                )
                cols = copy(cols)
                ll = zero(Prob)
                ll, S, E, I, R = singular_part!(
                    cols, geneal, node, ll,
                    S, E, I, R;
                    args...,
                )
                if isfinite(ll)
                    ll, S, E, I, R = regular_part!(
                        cols, ll, t, dt,
                        S, E, I, R;
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
        userdata = (geneal = gen,),
    )
end

end
