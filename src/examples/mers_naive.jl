module NaiveMERS

export mers

using ..PhyloPOMP
using ..PhyloPOMP: Root, Node, Sample, Name, Time

@demes Demes Camel Human
using .Demes: Camel, Human, DemeSet

include("mers_tree.jl")

## finite_log(x) = x > 0 ? log(x) : Float64(-Inf)
## safe_rate(x) = isfinite(x) && x > 0 ? Float64(x) : 0.0

function mers_singular!(
    cols, ll, geneal, node;
    Sc, Ic, Sh, Ih,
    Beta_cc, Beta_ch, Beta_hc, Beta_hh,
    chi_c, chi_h,
    Nc, Nh,
    _...,
    )
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
            k,lambda,p = rcateg([lambda_cc, lambda_hc, lambda_hc], false)
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
                ellc, ellh = fork!(cols, Human, n.lineage, (Human, Humnan), children)
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
    ll, Sc, Ic, Sh, Ih
end

function mers_event_rates!(
    rate, logpi, ellc, ellh;
    Sc, Ic, Sh, Ih,
    Beta_cc, Beta_ch, Beta_hc, Beta_hh,
    gamma_c, gamma_h, chi_c, chi_h, Bc, Bh, Nc, Nh,
    _...,
    )
    fill!(rate, 0.0)
    fill!(logpi, Float64(-Inf))
    penalty = 0.0

    alpha = Nc > 0 ? Beta_cc * Sc * Ic / Nc : 0.0
    disc = Ic > 0 ? ellc * (ellc - 1) / Ic / (Ic + 1) : 1.0
    penalty += alpha * disc
    rate[1] = safe_rate(alpha * (1 - disc)); logpi[1] = 0.0

    alpha = Nh > 0 ? Beta_hh * Sh * Ih / Nh : 0.0
    disc = Ih > 0 ? ellh * (ellh - 1) / Ih / (Ih + 1) : 1.0
    penalty += alpha * disc
    rate[2] = safe_rate(alpha * (1 - disc)); logpi[2] = 0.0

    alpha = Nh > 0 ? Beta_ch * Sc * Ih / Nh : 0.0
    pi = Ih > 0 ? (Ih - 0.5 * ellh) / Ih : 1.0
    rate[3] = safe_rate(alpha * pi); logpi[3] = log(pi)
    pi = Ih > 0 ? 1 - pi : 0.0
    rate[4] = safe_rate(alpha * pi); logpi[4] = (pi > 0 && ellh > 0) ? log(pi) - log(ellh) : Float64(-Inf)

    alpha = Nc > 0 ? Beta_hc * Sh * Ic / Nc : 0.0
    pi = Ic > 0 ? (Ic - 0.5 * ellc) / Ic : 1.0
    rate[5] = safe_rate(alpha * pi); logpi[5] = log(pi)
    pi = Ic > 0 ? 1 - pi : 0.0
    rate[6] = safe_rate(alpha * pi); logpi[6] = (pi > 0 && ellc > 0) ? log(pi) - log(ellc) : Float64(-Inf)

    alpha = gamma_c * Ic
    if Ic > ellc
        rate[7] = safe_rate(alpha)
    else
        penalty += alpha
    end
    logpi[7] = 0.0

    alpha = gamma_h * Ih
    if Ih > ellh
        rate[8] = safe_rate(alpha)
    else
        penalty += alpha
    end
    logpi[8] = 0.0

    penalty += chi_c * Ic + chi_h * Ih
    rate[9] = safe_rate(Bc); logpi[9] = 0.0
    rate[10] = safe_rate(Bh); logpi[10] = 0.0
    rate[11] = safe_rate(Bc); logpi[11] = 0.0
    rate[12] = Nh > 0 ? safe_rate(Bh / Nh * Sh) : 0.0; logpi[12] = 0.0

    penalty, sum(rate)
end

function mers_regular!(
    cols, ll;
    t, dt,
    Sc, Ic, Sh, Ih,
    args...,
    )
    tf = t + dt
    rate = zeros(Float64, 12)
    logpi = zeros(Float64, 12)
    ellc, ellh = ell(cols)

    while t < tf
        penalty, event_rate = mers_event_rates!(
            rate, logpi, ellc, ellh;
            Sc, Ic, Sh, Ih,
            args...,
        )

        if event_rate <= 0
            ll -= penalty * (tf - t)
            break
        end

        step = -log(rand()) / event_rate
        if t + step >= tf
            ll -= penalty * (tf - t)
            break
        end

        k, _ = rcateg(rate)
        ll -= penalty * step + logpi[k]

        if k == 1
            Sc -= 1; Ic += 1
        elseif k == 2
            Sh -= 1; Ih += 1
        elseif k == 3
            Sc -= 1; Ic += 1
            ll += log(1 - ellc / Ic)
        elseif k == 4
            Sc -= 1; Ic += 1
            b = rand(cols[Human])
            ellc, ellh = swap!(cols, Human, Camel, b)
            ll += log(1 - ellh / Ih) - log(Ic)
        elseif k == 5
            Sh -= 1; Ih += 1
            ll += log(1 - ellh / Ih)
        elseif k == 6
            Sh -= 1; Ih += 1
            b = rand(cols[Camel])
            ellc, ellh = swap!(cols, Camel, Human, b)
            ll += log(1 - ellc / Ic) - log(Ih)
        elseif k == 7
            Ic -= 1
        elseif k == 8
            Ih -= 1
        elseif k == 9
            Sc += 1
        elseif k == 10
            Sh += 1
        elseif k == 11
            Sc > 0 && (Sc -= 1)
        elseif k == 12
            Sh > 0 && (Sh -= 1)
        else
            @assert false "impossible event"
        end

        @assert Ic >= ellc && Ih >= ellh
        t += step
    end
    ll, Sc, Ic, Sh, Ih
end

"""
    mers(gen; ...)

Construct a Julia POMP object for the phylopomp MERS genealogy-conditioned
filter. Parameter names and event order follow R phylopomp's MERS model.
"""
function mers(
    Beta_cc = 4.0, Beta_ch = 0.0, Beta_hc = 0.0, Beta_hh = 4.0,
    gamma_c = 1.0, gamma_h = 1.0,
    chi_c = 1.0, chi_h = 0.0,
    Bc = 0.0, Bh = 0.0,
    Sc0 = 1.0, Sh0 = 1.0,
    Ic0 = 0.01, Ih0 = 0.0,
    Nc = 10000, Nh = 10000,
    )
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
                ll, Sc, Ic, Sh, Ih = mers_singular!(
                    cols, ll, geneal, node;
                    Sc, Ic, Sh, Ih,
                    args...,
                )
                if isfinite(ll)
                    ll, Sc, Ic, Sh, Ih = mers_regular!(
                        cols, ll;
                        t, dt,
                        Sc, Ic, Sh, Ih,
                        args...,
                    )
                end
                (; node = node + one(Name), ll, cols, Sc, Ic, Sh, Ih)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = parse_newick(mers_tree),)
    )
end
