# Exact log-likelihood of a genealogy under a linear table: every rate is proportional to the count of its `from` deme
# and the table has no compartments outside its demes (LBDP, BDEI, BDSS, MTBD). These are multi-type birth-death
# processes; the likelihood solves their backward equations (Stadler & Bonhoeffer 2013; Kühnert et al. 2016) along
# the genealogy with an adaptive Dormand–Prince Runge–Kutta method. `lbdp_exact` is the one-deme closed form.

export mtbd_exact, mtbd_rates

"""
    mtbd_rates(model, θ) -> (b, m, d, ψ, χ)

Per-host rates of a linear table at parameters `θ`, by deme in the order of `model.demes`:
`b[i, j]`, a host in deme `i` infects a new host in deme `j` (the infector stays in `i`); `m[i, j]`, a host moves
from `i` to `j`; `d[i]`, removal without sampling; `ψ[i]`, sampling that keeps the host; `χ[i]`, sampling that
removes it.

Throws `ArgumentError` when the table has compartments other than its demes, when a rate is not proportional to
the count of its event's `from` deme, or for an event that is none of these five.
"""
function mtbd_rates(model::MGPModel, θ::NamedTuple)
    K = length(model.demes)
    length(model.compartments) == K && Set(model.compartments) == Set(model.demes) ||
        throw(ArgumentError("mtbd_rates: table `$(model.name)` has compartments other than its demes"))
    comps = Tuple(model.compartments)
    pos = [findfirst(==(s), model.compartments) for s in model.demes]
    state(a, c, other) = NamedTuple{comps}(ntuple(k -> k == pos[a] ? c : other, length(comps)))
    b = zeros(K, K); m = zeros(K, K); d = zeros(K); ψ = zeros(K); χ = zeros(K)
    for ev in model.events
        a = ev.from
        h = Float64(ev.hazard(state(a, 1, 0), θ))
        isapprox(Float64(ev.hazard(state(a, 3, 1), θ)), 3h; rtol = 1e-12) ||
            throw(ArgumentError("mtbd_rates: the rate of `$(ev.name)` is not proportional to the count of its deme"))
        if ev.type == BIRTH && sum(ev.r) == 2 && length(ev.into) == 1
            b[a, only(ev.into)] += h
        elseif ev.type == MIGRATION && length(ev.into) == 1
            m[a, only(ev.into)] += h
        elseif ev.type == DEATH
            d[a] += h
        elseif ev.type == SAMPLE && sum(ev.r) == 0
            χ[a] += h
        elseif ev.type == SAMPLE && sum(ev.r) == 1 && ev.r[a] == 1
            ψ[a] += h
        elseif h != 0
            throw(ArgumentError("mtbd_rates: event `$(ev.name)` is not a birth, move, death or sample of one host"))
        end
    end
    (b = b, m = m, d = d, ψ = ψ, χ = χ)
end

"""
    mtbd_exact(g, model; θ, founder, rtol = 1e-10) -> Float64

Exact log-likelihood of genealogy `g` under the linear table `model` (see `mtbd_rates`) at parameters `θ`.
One host at `timezero(g)`, in deme `founder`: a deme name (`:I`) or a vector of probabilities over `model.demes`,
which gives the mixture of the likelihoods. The genealogy ends at `g.time`. Not conditioned on survival or on the
number of samples. The conventions are those of `mgp_filter_pomp` with `x0` that one host, so the filter's
likelihood estimate is unbiased for the exponential of this value, and of `lbdp_exact` with `n0 = 1`.

Backward from `g.time`, with `p_i` the probability that a host in deme `i` has no sampled descendant and `D_i` the
density of the genealogy below a lineage in deme `i`:
`p_i' = -δ_i p_i + d_i + Σ_j (b_ij p_i + m_ij) p_j` and `D_i' = -δ_i D_i + Σ_j b_ij (D_i p_j + p_i D_j) + m_ij D_j`,
`δ_i = Σ_j (b_ij + m_ij) + d_i + ψ_i + χ_i`, `p_i(g.time) = 1`. A tip starts `D_i = χ_i + ψ_i p_i`; a sampled ancestor
multiplies by `ψ_i`; a branch point combines its children as `D_i = Σ_j b_ij (D1_i D2_j + D1_j D2_i)`. A sample whose
deme the genealogy records counts only in that deme. Each lineage's `D` is rescaled after every step and the
scale factors are added to the log-likelihood.

`rtol` is the relative error allowed in each Runge–Kutta step. Throws `ArgumentError` unless `g` has exactly one
root with one child, branch points with two children and samples with at most one.
"""
function mtbd_exact(g::Genealogy, model::MGPModel; θ::NamedTuple, founder, rtol::Real = 1e-10)
    R = mtbd_rates(model, θ)
    K = length(model.demes)
    π0 = _founder_probs(founder, model)
    δ = [sum(R.b[i, :]) + sum(R.m[i, :]) + R.d[i] + R.ψ[i] + R.χ[i] for i in 1:K]
    cap = g.nsample + 2
    Y = zeros(K, cap + 1)              # column 1: p; columns 2..n1: the active lineages' D
    Y[:, 1] .= 1.0
    ws = _DP5Work(K, cap + 1)
    owner = zeros(Int, cap + 1)        # column -> genealogy node
    col = Dict{Int,Int}()              # genealogy node -> column
    n1 = 1
    logL = 0.0
    nroot = 0
    h = 0.01 / max(maximum(δ), 1e-3)
    tcur = g.time
    v = zeros(K)
    for n in reverse(eachindex(g))
        nd = g[n]
        Δτ = tcur - nd.slate
        Δτ >= 0 || throw(ArgumentError("mtbd_exact: genealogy nodes are not in time order"))
        if Δτ > 0
            h, ls = _dp5!(Y, n1, Δτ, h, ws, R, δ, rtol)
            logL += ls
        end
        tcur = nd.slate
        nch = length(nd.children)
        if nd.type == Root
            nch == 1 || throw(ArgumentError("mtbd_exact: root $n has $nch children"))
            c = col[only(nd.children)]
            logL += log(sum(π0[i] * Y[i, c] for i in 1:K))
            n1 = _drop!(Y, owner, col, c, n1)
            nroot += 1
            continue
        elseif nd.type == Sample && nch == 0
            for i in 1:K
                v[i] = R.χ[i] + R.ψ[i] * Y[i, 1]
            end
        elseif nd.type == Sample && nch == 1
            c = col[only(nd.children)]
            for i in 1:K
                v[i] = R.ψ[i] * Y[i, c]
            end
            n1 = _drop!(Y, owner, col, c, n1)
        elseif nd.type == Node && nch == 2
            c1, c2 = col[nd.children[1]], col[nd.children[2]]
            for i in 1:K
                s = 0.0
                for j in 1:K
                    s += R.b[i, j] * (Y[i, c1] * Y[j, c2] + Y[j, c1] * Y[i, c2])
                end
                v[i] = s
            end
            n1 = _drop!(Y, owner, col, max(c1, c2), n1)     # the larger column first: _drop! moves the last one
            n1 = _drop!(Y, owner, col, min(c1, c2), n1)
        else
            throw(ArgumentError("mtbd_exact: $(nd.type) node $n has $nch children"))
        end
        if !ismissing(nd.deme)
            length(instances(typeof(nd.deme))) == K ||
                throw(ArgumentError("mtbd_exact: the genealogy's demes are not the $K demes of `$(model.name)`"))
            for i in 1:K
                i == Int(nd.deme) || (v[i] = 0.0)
            end
        end
        s = maximum(v)
        s > 0 || return -Inf
        n1 += 1
        Y[:, n1] .= v ./ s
        logL += log(s)
        owner[n1] = n
        col[n] = n1
    end
    nroot == 1 && n1 == 1 || throw(ArgumentError("mtbd_exact: the genealogy must have one root ($nroot found)"))
    logL
end

_founder_probs(founder::Symbol, model::MGPModel) = begin
    k = findfirst(==(founder), model.demes)
    k === nothing && throw(ArgumentError("mtbd_exact: `$founder` is not a deme of `$(model.name)`"))
    [i == k ? 1.0 : 0.0 for i in eachindex(model.demes)]
end
_founder_probs(founder::AbstractVector{<:Real}, model::MGPModel) = begin
    length(founder) == length(model.demes) && all(>=(0), founder) && isapprox(sum(founder), 1; atol = 1e-12) ||
        throw(ArgumentError("mtbd_exact: founder must be probabilities over the $(length(model.demes)) demes"))
    collect(Float64, founder)
end

## Remove column `c` (its lineage has joined its parent): the last active column moves into it.
function _drop!(Y, owner, col, c, n1)
    delete!(col, owner[c])
    if c != n1
        Y[:, c] .= @view Y[:, n1]
        owner[c] = owner[n1]
        col[owner[c]] = c
    end
    owner[n1] = 0
    n1 - 1
end

## Right-hand side of the backward equations for columns 1:n1 of Y (column 1 is p).
function _mtbd_rhs!(dY, Y, n1, R, δ)
    K = size(Y, 1)
    b, m, d = R.b, R.m, R.d
    @inbounds for i in 1:K
        s = -δ[i] * Y[i, 1] + d[i]
        for j in 1:K
            s += (b[i, j] * Y[i, 1] + m[i, j]) * Y[j, 1]
        end
        dY[i, 1] = s
    end
    @inbounds for c in 2:n1, i in 1:K
        s = -δ[i] * Y[i, c]
        for j in 1:K
            s += b[i, j] * (Y[i, c] * Y[j, 1] + Y[i, 1] * Y[j, c]) + m[i, j] * Y[j, c]
        end
        dY[i, c] = s
    end
    nothing
end

struct _DP5Work
    k::NTuple{7,Matrix{Float64}}
    tmp::Matrix{Float64}
    ynew::Matrix{Float64}
end
_DP5Work(K, n) = _DP5Work(ntuple(_ -> zeros(K, n), 7), zeros(K, n), zeros(K, n))

## Dormand–Prince 5(4) coefficients (Hairer, Nørsett & Wanner, Solving ODEs I, Table 5.2).
const _DP_A = ((1 / 5,),
               (3 / 40, 9 / 40),
               (44 / 45, -56 / 15, 32 / 9),
               (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
               (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
               (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84))
const _DP_E = (71 / 57600, 0.0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40)

## Integrate columns 1:n1 of Y over a backward-time span Δτ with adaptive steps, starting from step h.
## After every accepted step each lineage column is divided by its largest entry; returns (next h, Σ log scales).
function _dp5!(Y, n1, Δτ, h, ws, R, δ, rtol)
    k, tmp, ynew = ws.k, ws.tmp, ws.ynew
    K = size(Y, 1)
    atol = rtol / 100
    τ = 0.0
    logs = 0.0
    _mtbd_rhs!(k[1], Y, n1, R, δ)
    while τ < Δτ
        last = h >= Δτ - τ
        hh = last ? Δτ - τ : h
        for s in 1:6
            a = _DP_A[s]
            @inbounds for c in 1:n1, i in 1:K
                acc = Y[i, c]
                for r in 1:s
                    acc += hh * a[r] * k[r][i, c]
                end
                (s < 6 ? tmp : ynew)[i, c] = acc
            end
            _mtbd_rhs!(k[s+1], s < 6 ? tmp : ynew, n1, R, δ)
        end
        err = 0.0
        @inbounds for c in 1:n1, i in 1:K
            e = 0.0
            for r in 1:7
                e += _DP_E[r] * k[r][i, c]
            end
            sc = atol + rtol * max(abs(Y[i, c]), abs(ynew[i, c]))
            err = max(err, abs(hh * e) / sc)
        end
        if err <= 1
            τ = last ? Δτ : τ + hh
            @inbounds for c in 1:n1, i in 1:K
                Y[i, c] = ynew[i, c]
            end
            for c in 2:n1
                s = maximum(@view Y[:, c])
                if s > 0
                    @views Y[:, c] ./= s
                    logs += log(s)
                end
            end
            _mtbd_rhs!(k[1], Y, n1, R, δ)        # rescaling changed Y, so recompute instead of reusing k[7]
            grow = err == 0 ? 5.0 : min(5.0, max(0.2, 0.9 * err^(-0.2)))
            last || (h = hh * grow)
        else
            h = hh * max(0.2, 0.9 * err^(-0.2))
        end
    end
    h, logs
end
