module NaiveMERS

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Prob, Time

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

const mers_tree = parse_newick(mers_newick,t0=0,demes=Demes)

singular_part!(
    cols, ll, geneal, node,
    Sc, Ic, Sh, Ih;
    βcc, βch, βhc, βhh,
    χc, χh,
    Nc, Nh,
    _...,
) = begin
    ellc, ellh = ell(cols)
    n = geneal[node]
    @assert Ic ≥ ellc && Ih ≥ ellh
    if n.type == Root
        @assert length(n.children)==1 "wrong number of children ($(length(n.children)) != 1) at root $(n.name), t=$(n.slate)"
        i, _, p = rcateg([Ic-ellc, Ih-ellh], DemeSet, true)
        ll -= log(p)
        if ismissing(i)
            ll = Prob(-Inf)
            ellc, ellh = plant!(cols,Camel,n.lineage)
            Ic += 1
        else
            ellc, ellh = plant!(cols, i, n.lineage)
        end
    elseif n.type == Sample
        @assert length(n.children)==0 "too many children ($(length(n.children)) > 0) at sample $(n.name), t=$(n.slate)"
        if n.lineage ∉ cols[n.deme]
            ll = Prob(-Inf)
            if n.deme == Camel
                ellc, ellh = swap!(cols,Human,Camel,n.lineage)
                Ic += 1
            else
                ellc, ellh = swap!(cols,Camel,Human,n.lineage)
                Ih += 1
            end
        end
        ellc, ellh = chop!(cols, n.deme, n.lineage)
        if n.deme == Camel
            ll += log(χc * Ic)
            Ic -= 1
        elseif n.deme == Human
            ll += log(χh * Ih)
            Ih -= 1
        else
            @assert false "impossible sample deme" # COV_EXCL_LINE
        end
    elseif n.type == Node
        @assert length(n.children)==2 "wrong number of children ($(length(n.children)) ≠ 2) at node $(n.name), t=$(n.time)"
        children = map(n.children) do i
            geneal[i].lineage
        end
        if n.lineage ∈ cols[Camel]
            λcc = Nc > 0 ? βcc * Sc * Ic / Nc : 0.0
            λhc = Nc > 0 ? βhc * Sh * Ic / Nc : 0.0
            k,_,p = rcateg([λcc, 0.5*λhc, 0.5*λhc], true)
            ll -= log(p)
            if k == 0
                ll = Prob(-Inf)
                ellc, ellh = fork!(cols, Camel, n.lineage, (Camel, Camel), children)
                Ic += 1
            elseif k == 1       # camel-camel
                ellc, ellh = fork!(cols, Camel, n.lineage, (Camel, Camel), children)
                if Sc > 0
                    Sc -= 1
                end
                Ic += 1
                ll += log(λcc) - log(Ic * (Ic - 1) / 2)
            else                # camel-human
                if k == 2
                    ellc, ellh = fork!(cols, Camel, n.lineage, (Camel, Human), children)
                elseif k == 3
                    ellc, ellh = fork!(cols, Camel, n.lineage, (Human, Camel), children)
                else
                    @assert "impossible rcateg output" # COV_EXCL_LINE
                end
                if Sh > 0
                    Sh -= 1
                end
                Ih += 1
                ll += log(λhc) - log(Ic * Ih)
            end
        elseif n.lineage ∈ cols[Human]
            λhh = Nh > 0 ? βhh * Sh * Ih / Nh : 0.0
            λch = Nh > 0 ? βch * Sc * Ih / Nh : 0.0
            k,_,p = rcateg([λhh, 0.5*λch, 0.5*λch], true)
            ll -= log(p)
            if k == 0
                ll = Prob(-Inf)
                ellc, ellh = fork!(cols, Human, n.lineage, (Human, Human), children)
                Ih += 1
            elseif k == 1       # human-human
                ellc, ellh = fork!(cols, Human, n.lineage, (Human, Human), children)
                if Sh > 0
                    Sh -= 1
                end
                Ih += 1
                ll += log(λhh) - log(Ih * (Ih - 1) / 2)
            else                # human-camel
                if k == 2
                    ellc, ellh = fork!(cols, Human, n.lineage, (Camel, Human), children)
                elseif k == 3
                    ellc, ellh = fork!(cols, Human, n.lineage, (Human, Camel), children)
                else
                    @assert "impossible rcateg output" # COV_EXCL_LINE
                end
                if Sc > 0
                    Sc -= 1
                end
                Ic += 1
                ll += log(λch) - log(Ic * Ih)
            end
        else
            @assert false "impossible node deme" # COV_EXCL_LINE
        end
    else
        @assert false "impossible node type" # COV_EXCL_LINE
    end
    @assert Ic ≥ ellc && Ih ≥ ellh
    ll, Sc, Ic, Sh, Ih
end

event_rates!(
    alpha, pi, cols,
    Sc, Ic, Sh, Ih;
    βcc, βch, βhc, βhh,
    γc, γh, χc, χh, Bc, Bh, Nc, Nh,
    _...,
) = begin
    ellc, ellh = ell(cols)
    @assert Ic ≥ ellc && Ih ≥ ellh
    alpha[1] = βcc*Sc*Ic/Nc
    alpha[2] = βhh*Sh*Ih/Nh
    alpha[3] = alpha[4] = βhc*Sh*Ic/Nc
    alpha[5] = alpha[6] = βch*Sc*Ih/Nh
    alpha[7] = @indicator(Ic > ellc, γc*(Ic-ellc))
    alpha[8] = @indicator(Ih > ellh, γh*(Ih-ellh))
    alpha[9] = Bc
    alpha[10] = Bh
    alpha[11] = Bc*Sc/Nc
    alpha[12] = Bh*Sh/Nh

    pi[1:2] .= one(Prob)
    pi[3] = @indicator(Ic > 0, 1-ellc/Ic)
    pi[4] = @indicator(Ic > 0, ellc/Ic)
    pi[5] = @indicator(Ih > 0, 1-ellh/Ih)
    pi[6] = @indicator(Ih > 0, ellh/Ih)
    pi[7:12] .= one(Prob)

    χc * Ic + χh * Ih +
        γc*ellc + @indicator(Ic ≤ ellc, γc*(Ic-ellc)) +
        γh*ellh + @indicator(Ih ≤ ellh, γh*(Ih-ellh))
end

regular_part!(
    cols, ll,
    t, dt,
    Sc, Ic, Sh, Ih;
    args...,
) = begin
    tf = t + dt
    if t < tf
        alpha = similar(Vector{Prob}, 12)
        pi = similar(Vector{Prob}, 12)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellc, ellh = ell(cols)
        while t+step < tf
            decay = event_rates!(
                alpha, pi, cols,
                Sc, Ic, Sh, Ih;
                args...,
            )
            k, s = rcateg(alpha .* pi)
            step = -log(rand())/s
            if t+step < tf
                ll -= decay*step+log(pi[k])
                if k == 1
                    Sc -= 1
                    Ic += 1
                    ll += log(1-(ellc*(ellc-1)/Ic/(Ic-1)))
                elseif k == 2
                    Sh -= 1
                    Ih += 1
                    ll += log(1-(ellh*(ellh-1)/Ih/(Ih-1)))
                elseif k == 3
                    Sh -= 1
                    Ih += 1
                    ll += log(1 - ellh / Ih)
                elseif k == 4
                    ll += log(ellc)
                    Sh -= 1
                    Ih += 1
                    b = rand(cols[Camel])
                    ellc, ellh = swap!(cols, Camel, Human, b)
                    ll += log((1 - ellc / Ic) / Ih)
                elseif k == 5
                    Sc -= 1
                    Ic += 1
                    ll += log(1 - ellc / Ic)
                elseif k == 6
                    ll += log(ellh)
                    Sc -= 1
                    Ic += 1
                    b = rand(cols[Human])
                    ellc, ellh = swap!(cols, Human, Camel, b)
                    ll += log((1 - ellh / Ih) / Ic)
                elseif k == 7
                    ll -= log(1-ellc/Ic)
                    Ic -= 1
                elseif k == 8
                    ll -= log(1-ellh/Ih)
                    Ih -= 1
                elseif k == 9
                    Sc += 1
                elseif k == 10
                    Sh += 1
                elseif k == 11
                    Sc -= 1
                elseif k == 12
                    Sh -= 1
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
        @assert Ic >= ellc && Ih >= ellh
    end
    ll, Sc, Ic, Sh, Ih
end

"""
    filter_pomp(; ...)

Construct a Julia POMP object for the phylopomp MERS genealogy-conditioned
filter. Parameter names and event order follow R phylopomp's MERS model.
"""
filter_pomp(
    ;βcc = 4.0, βch = 0.0, βhc = 1.0, βhh = 4.0,
    γc = 1.0, γh = 1.0,
    χc = 1.0, χh = 0.0,
    Bc = 0.1, Bh = 0.03,
    Sc0 = 1.0, Sh0 = 1.0,
    Ic0 = 0.01, Ih0 = 0.0,
    Nc = 10000, Nh = 10000,
) = begin
    gen = mers_tree
    pomp(
        params = map(
            Float64,
            (;βcc, βch, βhc, βhh,
             γc, γh, χc, χh, Bc, Bh,
             Sc0, Sh0, Ic0, Ih0, Nc, Nh,
             )
        ),
        t0 = timezero(gen),
        times = times(gen),
        rinit = function (; Sc0, Sh0, Ic0, Ih0, Nc, Nh, _...)
            fc = Nc / (Sc0 + Ic0)
            fh = Nh / (Sh0 + Ih0)
            (
                node = one(Name),
                ll = zero(Float64),
                cols = Coloring(Demes),
                Sc = round(Int64, fc * Sc0),
                Ic = round(Int64, fc * Ic0),
                Sh = round(Int64, fh * Sh0),
                Ih = round(Int64, fh * Ih0),
            )
        end,
        rprocess = onestep(
            function (; node, ll, cols, geneal,
                      Sc, Ic, Sh, Ih,
                      t, dt, args...,
                      )
                cols = copy(cols)
                ll = zero(Float64)
                ll, Sc, Ic, Sh, Ih = singular_part!(
                    cols, ll, geneal, node,
                    Sc, Ic, Sh, Ih;
                    args...,
                )
                if isfinite(ll)
                    ll, Sc, Ic, Sh, Ih = regular_part!(
                        cols, ll, t, dt,
                        Sc, Ic, Sh, Ih;
                        args...,
                    )
                end
                (; node = node + one(Name), ll, cols, Sc, Ic, Sh, Ih)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,)
    )
end

end
