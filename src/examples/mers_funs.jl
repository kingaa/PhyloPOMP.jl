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

## ---------------------------------------------------------------------
## Support-safe proposal-law helpers, shared by SoftMERS, GuidedMERS, and
## HardMERS. Each mixes a base (possibly degenerate) law with a uniform
## floor of mass `epsilon` so that every *target-feasible* outcome keeps
## strictly positive proposal mass -- this is the epsilon-floored mixture
##     q_epsilon(j) = (1-epsilon) q_base(j) + epsilon/(ell+1)
## described in `mers_filter_suite.tex`. Feasibility (a structurally zero
## population rate, e.g. an empty susceptible pool) is never overridden by
## the floor: an infeasible outcome keeps exactly zero proposal mass.
## ---------------------------------------------------------------------

"""
    floored_shares(w, epsilon)

Normalize the nonnegative weights `w` to shares summing to 1 (falling back
to a uniform law when every entry is zero or non-finite), then mix with a
uniform floor of mass `epsilon` so every entry is strictly positive.
`epsilon == 0` disables the floor entirely: the result is exactly the base
law (including its own all-zero-weight fallback to uniform, which is a
separate fallback and not governed by `epsilon`).
"""
floored_shares(w::AbstractVector{<:Real}, epsilon::Real) = begin
    @assert 0 <= epsilon <= 1
    n = length(w)
    @assert n > 0
    cleaned = map(x -> (isfinite(x) && x > 0) ? Float64(x) : 0.0, w)
    s = sum(cleaned)
    base = s > 0 ? cleaned ./ s : fill(1/n, n)
    (1-epsilon) .* base .+ epsilon/n
end

"""
    floored_shares_feasible(w, feasible, epsilon)

As `floored_shares`, but mixes only over the entries marked `true` in
`feasible`; entries marked `false` are structurally infeasible (their
target rate is zero) and receive exactly zero share, so the floor never
introduces proposal mass for a target-incompatible outcome. `epsilon == 0`
disables the floor: every feasible entry gets exactly its base share (with
that base's own zero-weight-among-feasible fallback to a uniform law over
the feasible set), and every infeasible entry is exactly zero, same as
with a positive `epsilon`.
"""
floored_shares_feasible(
    w::AbstractVector{<:Real}, feasible::AbstractVector{Bool}, epsilon::Real,
) = begin
    @assert 0 <= epsilon <= 1
    @assert length(w) == length(feasible)
    n = count(feasible)
    @assert n > 0 "no feasible outcomes to mix over"
    cleaned = map(x -> (isfinite(x) && x > 0) ? Float64(x) : 0.0, w)
    s = sum(cleaned[i] for i ∈ eachindex(cleaned) if feasible[i])
    map(eachindex(cleaned)) do i
        if !feasible[i]
            0.0
        else
            base = s > 0 ? cleaned[i]/s : 1/n
            (1-epsilon)*base + epsilon/n
        end
    end
end

"""
    cross_proposal(I, ell, epsilon)

Aggregate identity/tracked-branch mixture probabilities for a cross-deme
transmission with `I` source hosts, `ell` of them tracked, mixed with a
uniform floor of mass `epsilon` over the `ell+1` reduced outcomes. Returns
`(q0, qsum)`: `q0` is the identity-outcome probability and `qsum` is the
*total* probability of the tracked-branch group. `ell == 0` returns the
degenerate law `(1, 0)` (no branches to track, so "identity" is certain).
This is the MERS analogue of `NaiveMERS.cross_proposal` (same formula),
shared here so `SoftMERS` and `HardMERS` reuse it for their identity-vs-
tracked-group split. `epsilon == 0` disables the floor: `(q0, qsum)`
reduces to the plain naive split `((I-ell)/I, ell/I)`.
"""
cross_proposal(I, ell, epsilon) = begin
    @assert I >= ell >= 0
    @assert 0 <= epsilon <= 1
    ell == 0 && return (one(Prob), zero(Prob))
    q0 = (1-epsilon)*(I-ell)/I + epsilon/(ell+1)
    qb = (1-epsilon)/I + epsilon/(ell+1)
    q0, ell*qb
end

"""
    floored_branch_law(rh, ell, I, epsilon)

Per-branch mixture probabilities `q_epsilon(b)`, `b = 1,...,ell`, for the
tracked-branch outcomes of a cross-deme transmission with `ell` tracked
source lineages (relative hazards `rh`) out of `I` source hosts. Preserves
the naive aggregate tracked-branch mass `ell/I` (the branches' relative
hazards only redistribute *within* that mass), then applies the
epsilon-floor mixture over all `ell+1` reduced outcomes; `sum(floored_branch_law(rh,ell,I,epsilon)) == cross_proposal(I,ell,epsilon)[2]`.
Falls back to a uniform conditional law when every hazard is zero or
non-finite. This is `SoftMERS`'s branch law: guide-weighted *within* the
naive split, never guide-weighted *across* the identity/tracked split.
`epsilon == 0` disables the floor: each branch gets exactly `(ell/I)*share`
with no mixture term (the zero-hazard fallback to a uniform conditional
law is separate and still applies).
"""
floored_branch_law(rh::AbstractVector{<:Real}, ell::Integer, I::Real, epsilon::Real) = begin
    @assert 0 <= epsilon <= 1
    ell == 0 && return Float64[]
    @assert I >= ell > 0
    cleaned = map(r -> (isfinite(r) && r > 0) ? Float64(r) : 0.0, rh)
    s = sum(cleaned)
    share = s > 0 ? cleaned ./ s : fill(1/ell, ell)
    (1-epsilon) .* (ell/I) .* share .+ epsilon/(ell+1)
end

"""
    hard_branch_law(rh, ell, I, epsilon)

Per-branch auxiliary mixture intensities `kappa_epsilon(b)`, `b =
1,...,ell`, using the *unnormalized* relative hazards `rh` directly as
intensity multipliers, so `sum(hard_branch_law(...))` may differ from the
naive tracked-branch mass `ell/I` -- the excess or deficit relative to the
population rate must be returned to the population process via the
compensating decay (`alpha_population - sum(beta_j)`) where this is used.
Every entry is at least `epsilon/(ell+1) > 0`, even when the corresponding
hazard is zero or non-finite, so this never assigns zero mass to a
target-feasible outcome. This is `HardMERS`'s branch law, distinct from
`floored_branch_law` (`SoftMERS`'s) precisely because it is not
renormalized to preserve `ell/I`. `epsilon == 0` disables the floor: a
branch with zero (or non-finite) relative hazard gets exactly zero
auxiliary intensity, rather than the strictly-positive floor value --
i.e. that outcome is genuinely allowed to fail.
"""
hard_branch_law(rh::AbstractVector{<:Real}, ell::Integer, I::Real, epsilon::Real) = begin
    @assert 0 <= epsilon <= 1
    ell == 0 && return Float64[]
    @assert I >= ell > 0
    cleaned = map(r -> (isfinite(r) && r > 0) ? Float64(r) : 0.0, rh)
    (1-epsilon) .* cleaned ./ I .+ epsilon/(ell+1)
end

singular_part!(
    cols, guide, geneal, node, ll,
    S_c, I_c, S_h, I_h;
    β_cc, β_ch, β_hc, β_hh,
    χ_c, χ_h, N_c, N_h,
    proposal_floor = 0.05,
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
                ## Support-safe: mix the guide's root-deme belief with a
                ## uniform floor over whichever deme(s) actually have free
                ## (untracked) hosts, so a vanishing guide weight for a
                ## feasible deme never zeroes its proposal mass.
                feasible = [freeC > 0, freeH > 0]
                shares = floored_shares_feasible(n.present[:,1], feasible, proposal_floor)
                i, _, p = rcateg(shares, DemeSet, true)
                ll -= log(p)
                ellC, ellH = plant!(cols, i, n.chillins[1])
            else
                ## incompatible with the data, but the coloring must still
                ## be corrected to avoid downstream errors.
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
            ## incompatible with the data; correct the coloring.
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
            ## MERS sampling is destructive (rate χ).
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
            if rate_cc <= 0 && rate_hc <= 0
                ## No target-feasible fork outcome (S_c = S_h = 0): incompatible
                ## with the data, but the coloring must still be corrected to
                ## avoid downstream errors.
                ll += Prob(-Inf)
                ellC, ellH = fork!(cols, Camel, n.parlin, (Camel, Camel), n.chillins)
                I_c += 1
            else
                ## Weight the labelling by (population rate) × (guide present-prob).
                ## Present rows are indexed [Camel=1, Human=2]; columns are the
                ## two children. These are proposal weights; the `-log(p)` term below
                ## removes the proposal so the target weight remains unchanged.
                ## Support-safe: the epsilon floor mixes only over the outcomes
                ## whose *population rate* is nonzero (`feasible`), so a
                ## vanishing guide weight never zeroes a target-feasible
                ## outcome, while a structurally impossible one (rate 0) is
                ## never assigned proposal mass by the floor.
                g1 = n.present[1,1] * n.present[1,2]   # (Camel,Camel)
                g2 = n.present[1,1] * n.present[2,2]   # (Camel,Human)
                g3 = n.present[2,1] * n.present[1,2]   # (Human,Camel)
                feasible = [rate_cc > 0, rate_hc > 0, rate_hc > 0]
                shares = floored_shares_feasible(
                    [rate_cc*g1, rate_hc*g2, rate_hc*g3], feasible, proposal_floor,
                )
                k, _, p = rcateg(shares, true)
                ll -= log(p)
                if k==1
                    @assert S_c > 0
                    ellC, ellH = fork!(cols, Camel, n.parlin, (Camel, Camel), n.chillins)
                    S_c -= 1
                    I_c += 1
                    ll += log(rate_cc) - log(I_c*(I_c-1)/2)
                else
                    @assert S_h > 0
                    if k==2
                        ellC, ellH = fork!(cols, Camel, n.parlin, (Camel, Human), n.chillins)
                    elseif k==3
                        ellC, ellH = fork!(cols, Camel, n.parlin, (Human, Camel), n.chillins)
                    else
                        @assert false "impossible rcateg output" # COV_EXCL_LINE
                    end
                    S_h -= 1
                    I_h += 1
                    ll += log(rate_hc) - log(I_c*I_h)
                end
            end
        elseif n.parlin ∈ cols[Human]
            rate_hh = N_h > 0 ? β_hh*S_h*I_h/N_h : 0.0
            rate_ch = N_h > 0 ? β_ch*S_c*I_h/N_h : 0.0
            if rate_hh <= 0 && rate_ch <= 0
                ll += Prob(-Inf)
                ellC, ellH = fork!(cols, Human, n.parlin, (Human, Human), n.chillins)
                I_h += 1
            else
                g1 = n.present[2,1] * n.present[2,2]   # (Human,Human)
                g2 = n.present[1,1] * n.present[2,2]   # (Camel,Human)
                g3 = n.present[2,1] * n.present[1,2]   # (Human,Camel)
                feasible = [rate_hh > 0, rate_ch > 0, rate_ch > 0]
                shares = floored_shares_feasible(
                    [rate_hh*g1, rate_ch*g2, rate_ch*g3], feasible, proposal_floor,
                )
                k, _, p = rcateg(shares, true)
                ll -= log(p)
                if k==1
                    @assert S_h > 0
                    ellC, ellH = fork!(cols, Human, n.parlin, (Human, Human), n.chillins)
                    S_h -= 1
                    I_h += 1
                    ll += log(rate_hh) - log(I_h*(I_h-1)/2)
                else
                    @assert S_c > 0
                    if k==2
                        ellC, ellH = fork!(cols, Human, n.parlin, (Camel, Human), n.chillins)
                    elseif k==3
                        ellC, ellH = fork!(cols, Human, n.parlin, (Human, Camel), n.chillins)
                    else
                        @assert false "impossible rcateg output" # COV_EXCL_LINE
                    end
                    S_c -= 1
                    I_c += 1
                    ll += log(rate_ch) - log(I_c*I_h)
                end
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
    filter_pomp(gen, m; β_cc = 4.0, ..., proposal_floor = 0.05)

Constructs a pomp object for the MERS genealogy-conditioned filter, based on
the filter guide built from genealogy `gen` and the guiding finite-state
Markov process `m` (construct `m` with `fsmarkov`). `proposal_floor` is the
uniform-mixture mass `epsilon` used throughout `singular_part!` and each
kernel's `regular_part!` to guarantee support-safe proposals at
population/lineage boundaries; it must satisfy `0 <= proposal_floor <= 1`.
Setting `proposal_floor = 0` disables the floor entirely -- every helper's
`epsilon == 0` case reduces exactly to the un-floored base law, so a
target-feasible outcome with vanishing guide/hazard weight can once again
receive zero proposal mass and the particle genuinely dies, rather than
being kept alive by manufactured mass. This is the "no floor, let it fail"
mode; the default `0.05` keeps the support-safe behavior described in
`mers_filter_suite.tex` §10.1.

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
    proposal_floor = 0.05,
) = begin
    0 <= proposal_floor <= 1 || throw(ArgumentError(
        "proposal_floor must satisfy 0 <= proposal_floor <= 1, got $proposal_floor",
    ))
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
            proposal_floor = Float64(proposal_floor),
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
