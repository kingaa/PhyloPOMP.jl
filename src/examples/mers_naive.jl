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
    Beta_cc, Beta_ch, Beta_hc, Beta_hh,
    chi_c, chi_h,
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
            ll += log(chi_c * Ic)
            Ic -= 1
        elseif n.deme == Human
            ll += log(chi_h * Ih)
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
            lambda_cc = Nc > 0 ? Beta_cc * Sc * Ic / Nc : 0.0
            lambda_hc = Nc > 0 ? Beta_hc * Sh * Ic / Nc : 0.0
            k,_,p = rcateg([lambda_cc, 0.5*lambda_hc, 0.5*lambda_hc], true)
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
                ll += log(lambda_cc) - log(Ic * (Ic - 1) / 2)
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
                ll += log(lambda_hc) - log(Ic * Ih)
            end
        elseif n.lineage ∈ cols[Human]
            lambda_hh = Nh > 0 ? Beta_hh * Sh * Ih / Nh : 0.0
            lambda_ch = Nh > 0 ? Beta_ch * Sc * Ih / Nh : 0.0
            k,_,p = rcateg([lambda_hh, 0.5*lambda_ch, 0.5*lambda_ch], true)
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
                ll += log(lambda_hh) - log(Ih * (Ih - 1) / 2)
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
                ll += log(lambda_ch) - log(Ic * Ih)
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
    Beta_cc, Beta_ch, Beta_hc, Beta_hh,
    gamma_c, gamma_h, chi_c, chi_h, Bc, Bh, Nc, Nh,
    _...,
) = begin
    ellc, ellh = ell(cols)
    @assert Ic ≥ ellc && Ih ≥ ellh
    alpha[1] = Beta_cc*Sc*Ic/Nc
    alpha[2] = Beta_hh*Sh*Ih/Nh
    alpha[3] = alpha[4] = Beta_hc*Sh*Ic/Nc
    alpha[5] = alpha[6] = Beta_ch*Sc*Ih/Nh
    alpha[7] = @indicator(Ic > ellc, gamma_c*(Ic-ellc))
    alpha[8] = @indicator(Ih > ellh, gamma_h*(Ih-ellh))
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

    chi_c * Ic + chi_h * Ih +
        gamma_c*ellc + @indicator(Ic ≤ ellc, gamma_c*(Ic-ellc)) +
        gamma_h*ellh + @indicator(Ih ≤ ellh, gamma_h*(Ih-ellh))
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
    ;Beta_cc = 4.0, Beta_ch = 0.0, Beta_hc = 1.0, Beta_hh = 4.0,
    gamma_c = 1.0, gamma_h = 1.0,
    chi_c = 1.0, chi_h = 0.0,
    Bc = 0.1, Bh = 0.03,
    Sc0 = 1.0, Sh0 = 1.0,
    Ic0 = 0.01, Ih0 = 0.0,
    Nc = 10000, Nh = 10000,
) = begin
    gen = mers_tree
    pomp(
        params = (
            Beta_cc = Float64(Beta_cc), Beta_ch = Float64(Beta_ch),
            Beta_hc = Float64(Beta_hc), Beta_hh = Float64(Beta_hh),
            gamma_c = Float64(gamma_c), gamma_h = Float64(gamma_h),
            chi_c = Float64(chi_c), chi_h = Float64(chi_h),
            Bc = Float64(Bc), Bh = Float64(Bh),
            Sc0 = Float64(Sc0), Sh0 = Float64(Sh0),
            Ic0 = Float64(Ic0), Ih0 = Float64(Ih0),
            Nc = Float64(Nc), Nh = Float64(Nh),
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
