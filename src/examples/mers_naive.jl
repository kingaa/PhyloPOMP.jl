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
    if Ic < ellc || Ih < ellh
        ll += Prob(-Inf)
    elseif n.type == Root
        if length(n.children) == 1
            freec = Ic - ellc
            freeh = Ih - ellh
            if freec + freeh > 0
                i, _, p = rcateg([freec, freeh], DemeSet, true)
                ll -= log(p)
                ellc, ellh = plant!(cols, i, n.lineage)
            else
                ll += Prob(-Inf)
                ellc, ellh = plant!(cols,Camel,n.lineage)
                Ic += 1
            end
        else
            error("too many children ($(length(n.children)) > 1) at root $(n.name)")
        end
    elseif n.type == Sample
        deme = n.deme
        if ismissing(deme) # FIXME: trap for this problem at an earlier stage
            error("MERS samples must have deme metadata Camel or Human")
        end
        if n.lineage ∉ cols[deme]
            ll += Prob(-Inf)
            if deme == Camel
                ellE, ellI = swap!(cols,Human,Camel,n.lineage)
                Ic += 1
                Ih -= 1
            elseif deme == Human
                ellE, ellI = swap!(cols,Camel,Human,n.lineage)
                Ic -= 1
                Ih += 1
            else
                @assert false "impossible sample deme" # COV_EXCL_LINE
            end
        end
        if length(n.children) == 0
            ellc, ellh = chop!(cols, deme, n.lineage)
            if deme == Camel
                ll += log(chi_c * Ic)
                Ic -= 1
            elseif deme == Human
                ll += log(chi_h * Ih)
                Ih -= 1
            else
                @assert false "impossible sample deme" # COV_EXCL_LINE
            end
        else
            error("MERS sampling is destructive but sample $(n.name) has $(length(n.children)) children")
        end
    elseif n.type == Node
        if length(n.children) != 2
            error("too many children ($(length(n.children)) != 2) at node $(n.name)")
        end
        children = map(n.children) do i
            geneal[i].lineage
        end
        if n.lineage ∈ cols[Camel]
            lambda_cc = Nc > 0 ? Beta_cc * Sc * Ic / Nc : 0.0
            lambda_hc = Nc > 0 ? Beta_hc * Sh * Ic / Nc / 2.0 : 0.0
            k,lambda = rcateg([lambda_cc, lambda_hc, lambda_hc], false)
            if k == 1
                @assert Sc > 0
                ellc, ellh = fork!(cols, Camel, n.lineage, (Camel, Camel), children)
                Sc -= 1
                Ic += 1
                ll += log(lambda) - log(Ic * (Ic - 1) / 2)
            else
                @assert Sh > 0
                if k == 2
                    ellc, ellh = fork!(cols, Camel, n.lineage, (Camel, Human), children)
                elseif k == 3
                    ellc, ellh = fork!(cols, Camel, n.lineage, (Human, Camel), children)
                else
                    @assert "impossible rcateg output" # COV_EXCL_LINE
                end
                Sh -= 1
                Ih += 1
                ll += log(lambda) - log(Ic * Ih)
            end
        elseif n.lineage ∈ cols[Human]
            lambda_hh = Nh > 0 ? Beta_hh * Sh * Ih / Nh : 0.0
            lambda_ch = Nh > 0 ? Beta_ch * Sc * Ih / Nh / 2.0 : 0.0
            k,lambda,p = rcateg([lambda_hh, lambda_ch, lambda_ch], true)
            if k == 1
                @assert Sh > 0
                ellc, ellh = fork!(cols, Human, n.lineage, (Human, Human), children)
                Sh -= 1
                Ih += 1
                ll += log(lambda) - log(Ih * (Ih - 1) / 2)
            else
                @assert Sc > 0
                if k == 2
                    ellc, ellh = fork!(cols, Human, n.lineage, (Camel, Human), children)
                elseif k == 3
                    ellc, ellh = fork!(cols, Human, n.lineage, (Human, Camel), children)
                else
                    @assert "impossible rcateg output" # COV_EXCL_LINE
                end
                Sc -= 1
                Ic += 1
                ll += log(lambda) - log(Ic * Ih)
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
    alpha[10] = alpha[9] = Bc
    alpha[12] = alpha[11] = Bh

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
        ellc, ellh = ell(cols)
        decay = event_rates!(
            alpha, pi, cols,
            Sc, Ic, Sh, Ih;
            args...,
        )
        k, s = rcateg(alpha .* pi)
        step::Time = -log(rand())/s
        @assert Ic >= ellc && Ih >= ellh
        while t+step < tf
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
            decay = event_rates!(
                alpha, pi, cols,
                Sc, Ic, Sh, Ih;
                args...,
            )
            k, s = rcateg(alpha .* pi)
            step = -log(rand())/s
        end
        step = tf - t
        ll -= decay*step
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
    ;Beta_cc = 4.0, Beta_ch = 0.0, Beta_hc = 0.0, Beta_hh = 4.0,
    gamma_c = 1.0, gamma_h = 1.0,
    chi_c = 1.0, chi_h = 0.0,
    Bc = 0.0, Bh = 0.0,
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
