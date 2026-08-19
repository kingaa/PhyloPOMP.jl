## This file contains definitions that are used in the guided filters
## for the MERS-CoV two-host (Camel/Human) model.
##
## It is the MERS analogue of seir_funs.jl and is `include`d by both
## GuidedMERS (mers_guided.jl) and HardMERS (mers_hard.jl).
##
## The two structural differences from seir_funs.jl are:
##
##   1. A guide node (`GuideNode`) carries no `.deme` field, but a MERS
##      Sample's host species (Camel or Human) is *essential* data (it
##      selects χ_c vs χ_h and which color to chop). We therefore pass
##      the underlying `Genealogy` alongside the `Guide` and read the
##      sampled deme from `geneal[node].deme`. `guide[node]` and
##      `geneal[node]` share the same node index.
##
##   2. An internal Node (coalescence) is NOT deme-fixed in MERS: the
##      parent can transmit within-host (cc, hh) or across-host (hc, ch).
##      Hence `knowledge!` fixes the guide only at Samples --- exactly the
##      behavior of PhyloPOMP's default `known_deme!`.

"""
    knowledge!(v; deme, type, time)

Return `true` if the guide probabilities are fixed at this node (and, if so,
fill `v` with the appropriate probability vector); return `false` otherwise.

For MERS the deme is fixed exactly when the genealogy supplies host
metadata, i.e. at the sample tips. This is identical to the default
`known_deme!`.
"""
knowledge!(
    v; deme, type, time,
) = begin
    if ismissing(deme)
        false
    else
        demekron!(v, deme)
        true
    end
end

singular_part!(
    cols, guide, geneal, node, ll,
    S_c, I_c, S_h, I_h;
    β_cc, β_ch, β_hc, β_hh,
    χ_c, χ_h, N_c, N_h,
    _...,
) = begin
    ellC, ellH = ell(cols)
    n = guide[node]
    @assert I_c ≥ ellC && I_h ≥ ellH
    if n.type==Root
        if length(n.chillins)==1
            freeC = I_c - ellC
            freeH = I_h - ellH
            if freeC + freeH > 0
                i, _, p = rcateg(n.present[:,1].*[freeC, freeH], DemeSet, true)
                ll -= log(p)
                if ismissing(i)
                    ll = Prob(-Inf)
                    ellC, ellH = plant!(cols, Camel, n.chillins[1])
                    I_c += 1
                else
                    ellC, ellH = plant!(cols, i, n.chillins[1])
                end
            else
                ll += Prob(-Inf)
                ellC, ellH = plant!(cols, Camel, n.chillins[1])
                I_c += 1
            end
        else
            error("too many children ($(length(n.chillins)) > 1) at root, node $node")
        end
    elseif n.type==Sample
        deme = geneal[node].deme
        if ismissing(deme)
            error("MERS samples must have deme metadata Camel or Human (node $node)")
        end
        if n.parlin ∉ cols[deme]
            ll += Prob(-Inf)
            if deme==Camel
                ellC, ellH = swap!(cols, Human, Camel, n.parlin)
                I_c += 1
                I_h -= 1
            elseif deme==Human
                ellC, ellH = swap!(cols, Camel, Human, n.parlin)
                I_c -= 1
                I_h += 1
            else
                @assert false "impossible sample deme" # COV_EXCL_LINE
            end
        end
        if length(n.chillins)==0
            ellC, ellH = chop!(cols, deme, n.parlin)
            if deme==Camel
                ll += log(χ_c * I_c)
                I_c -= 1
            elseif deme==Human
                ll += log(χ_h * I_h)
                I_h -= 1
            else
                @assert false "impossible sample deme" # COV_EXCL_LINE
            end
        else
            error("MERS sampling is destructive but sample at node $node has $(length(n.chillins)) children")
        end
    elseif n.type==Node
        if length(n.chillins) ≠ 2
            error("wrong number of children ($(length(n.chillins)) ≠ 2) at node $node")
        end
        if n.parlin ∈ cols[Camel]
            rate_cc = N_c > 0 ? β_cc*S_c*I_c/N_c : 0.0
            rate_hc = N_c > 0 ? β_hc*S_h*I_c/N_c : 0.0
            g1 = n.present[1,1] * n.present[1,2]   # (Camel,Camel)
            g2 = n.present[1,1] * n.present[2,2]   # (Camel,Human)
            g3 = n.present[2,1] * n.present[1,2]   # (Human,Camel)
            k, _, p = rcateg([rate_cc*g1, rate_hc*g2, rate_hc*g3], true)
            ll -= log(p)
            if k==0
                ll = Prob(-Inf)
                ellC, ellH = fork!(cols, Camel, n.parlin, (Camel, Camel), n.chillins)
                I_c += 1
            elseif k==1
                ellC, ellH = fork!(cols, Camel, n.parlin, (Camel, Camel), n.chillins)
                if S_c > 0
                    S_c -= 1
                end
                I_c += 1
                ll += log(rate_cc) - log(I_c*(I_c-1)/2)
            else
                if k==2
                    ellC, ellH = fork!(cols, Camel, n.parlin, (Camel, Human), n.chillins)
                elseif k==3
                    ellC, ellH = fork!(cols, Camel, n.parlin, (Human, Camel), n.chillins)
                else
                    @assert false "impossible rcateg output" # COV_EXCL_LINE
                end
                if S_h > 0
                    S_h -= 1
                end
                I_h += 1
                ll += log(rate_hc) - log(I_c*I_h)
            end
        elseif n.parlin ∈ cols[Human]
            rate_hh = N_h > 0 ? β_hh*S_h*I_h/N_h : 0.0
            rate_ch = N_h > 0 ? β_ch*S_c*I_h/N_h : 0.0
            g1 = n.present[2,1] * n.present[2,2]   # (Human,Human)
            g2 = n.present[1,1] * n.present[2,2]   # (Camel,Human)
            g3 = n.present[2,1] * n.present[1,2]   # (Human,Camel)
            k, _, p = rcateg([rate_hh*g1, rate_ch*g2, rate_ch*g3], true)
            ll -= log(p)
            if k==0
                ll = Prob(-Inf)
                ellC, ellH = fork!(cols, Human, n.parlin, (Human, Human), n.chillins)
                I_h += 1
            elseif k==1
                ellC, ellH = fork!(cols, Human, n.parlin, (Human, Human), n.chillins)
                if S_h > 0
                    S_h -= 1
                end
                I_h += 1
                ll += log(rate_hh) - log(I_h*(I_h-1)/2)
            else
                if k==2
                    ellC, ellH = fork!(cols, Human, n.parlin, (Camel, Human), n.chillins)
                elseif k==3
                    ellC, ellH = fork!(cols, Human, n.parlin, (Human, Camel), n.chillins)
                else
                    @assert false "impossible rcateg output" # COV_EXCL_LINE
                end
                if S_c > 0
                    S_c -= 1
                end
                I_c += 1
                ll += log(rate_ch) - log(I_c*I_h)
            end
        else
            @assert false "impossible node deme" # COV_EXCL_LINE
        end
    else
        @assert false "impossible node type" # COV_EXCL_LINE
    end
    @assert I_c ≥ ellC && I_h ≥ ellH
    ll, S_c, I_c, S_h, I_h
end

"""
    filter_pomp(gen, m; β_cc = 4.0, ...)

Constructs a pomp object for the MERS genealogy-conditioned filter, based on
the filter guide built from genealogy `gen` and the guiding finite-state
Markov process `m` (construct `m` with `fsmarkov`).

Pass the built-in `mers_tree` as `gen` to reproduce the default data set; a call
`filter_pomp(mers_tree, fsmarkov(...))` constructs the guided filter.
"""
filter_pomp(
    gen::Genealogy,
    m::FSMarkovProc;
    β_cc = 4.0, β_ch = 0.0, β_hc = 0.0, β_hh = 4.0,
    γ_c = 1.0, γ_h = 1.0,
    χ_c = 1.0, χ_h = 0.0,
    B_c = 0.0, B_h = 0.0,
    S_c0 = 1.0, S_h0 = 1.0,
    I_c0 = 0.01, I_h0 = 0.0,
    N_c = 10000, N_h = 10000,
) = begin
    guidegen = guide(gen, m, knowledge!)
    pomp(
        params = (
            β_cc = Float64(β_cc), β_ch = Float64(β_ch),
            β_hc = Float64(β_hc), β_hh = Float64(β_hh),
            γ_c = Float64(γ_c), γ_h = Float64(γ_h),
            χ_c = Float64(χ_c), χ_h = Float64(χ_h),
            B_c = Float64(B_c), B_h = Float64(B_h),
            S_c0 = Float64(S_c0), S_h0 = Float64(S_h0),
            I_c0 = Float64(I_c0), I_h0 = Float64(I_h0),
            N_c = Float64(N_c), N_h = Float64(N_h),
        ),
        t0 = timezero(guidegen),
        times = times(guidegen),
        rinit = function (; S_c0, S_h0, I_c0, I_h0, N_c, N_h, _...)
            m_c = N_c / (S_c0 + I_c0)
            m_h = N_h / (S_h0 + I_h0)
            (
                node = one(Name),
                ll = zero(Prob),
                cols = Coloring(Demes),
                S_c = round(Int64, m_c*Float64(S_c0)),
                I_c = round(Int64, m_c*Float64(I_c0)),
                S_h = round(Int64, m_h*Float64(S_h0)),
                I_h = round(Int64, m_h*Float64(I_h0)),
            )
        end,
        rprocess = onestep(
            function (
                ; node, ll, cols, guide, geneal,
                S_c, I_c, S_h, I_h,
                kwargs...,
                )
                cols = copy(cols)
                ll = zero(Prob)
                ll, S_c, I_c, S_h, I_h = singular_part!(
                    cols, guide, geneal, node, ll,
                    S_c, I_c, S_h, I_h;
                    kwargs...,
                )
                if isfinite(ll)
                    ll, S_c, I_c, S_h, I_h = regular_part!(
                        cols, guide, node, ll,
                        S_c, I_c, S_h, I_h;
                        kwargs...,
                    )
                end
                (; node = node+1, ll = ll, cols = cols,
                 S_c = S_c, I_c = I_c, S_h = S_h, I_h = I_h)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (guide = guidegen, geneal = gen),
    )
end
