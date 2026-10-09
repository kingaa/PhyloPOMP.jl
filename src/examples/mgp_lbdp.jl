# Linear birth-death-sampling process, as R phylopomp's LBDP. One deme (n), no susceptibles.
# psi: sampling that keeps the host; chi: sampling that removes it.
# R's lbdp_exact gives the exact log-likelihood of a genealogy under this model.

@mgp LBDP begin
    compartments = (n,)
    demes = (n,)
    params = (λ, μ, ψ, χ)

    @event birth           rate=λ*n pop=(n=+1) move=fork(n => n, n)
    @event death           rate=μ*n pop=(n=-1) move=chop(n)
    @event sampling        rate=ψ*n pop=()     move=sample(n)        kind=singular
    @event sampling_remove rate=χ*n pop=(n=-1) move=sample_remove(n) kind=singular
end

export lbdp_exact

"""
    lbdp_exact(g; λ, μ, ψ, χ = 0.0, n0 = 1) -> Float64

Exact log-likelihood of genealogy `g` under `LBDP`, ported from R phylopomp's `lbdp_exact`.
`n0` hosts at `timezero(g)`; the genealogy ends at `g.time`. With `nroot` roots it includes
`log(n0!/(n0 - nroot)!)`, a factor the filters do not have, so compare filters with `n0 = 1`.
Throws `ArgumentError` when the number of tip samples is not the number of branch points plus roots.
"""
function lbdp_exact(g::Genealogy; λ::Real, μ::Real, ψ::Real, χ::Real = 0.0, n0::Integer = 1)
    tf = g.time
    t0 = timezero(g)
    tbr = Float64[]; ttp = Float64[]; nrt = 0; nin = 0
    for i in eachindex(g)
        nd = g[i]
        if nd.type == Root
            nrt += 1
        elseif nd.type == Sample
            isempty(nd.children) ? push!(ttp, nd.slate) : (nin += 1)
        else
            push!(tbr, nd.slate)
        end
    end
    length(ttp) == length(tbr) + nrt ||
        throw(ArgumentError("lbdp_exact: $(length(ttp)) tip samples, $(length(tbr)) branch points, $nrt roots"))
    a = λ - μ + ψ + χ
    b = λ - μ - ψ - χ
    d = sqrt(b * b + 4 * λ * (ψ + χ))
    ## R's G = (C + aS)/(C + bS), C = d·cosh(w), S = sinh(w), and H = 1/(cosh(w) + (b/d)·sinh(w))², w = d(t - tf)/2 ≤ 0,
    ## with cosh and sinh divided out (they overflow once d(tf - t)/2 passes about 710, e.g. R0 near 1 with small p):
    ## G = ((d - a) + (d + a)e)/((d - b) + (d + b)e) and log H = 2w + log 4 - 2·log(((d - b) + (d + b)e)/d), e = exp(2w).
    ## d ± a and d ± b come from d² - a² = 4μ(ψ + χ) and d² - b² = 4λ(ψ + χ), without cancellation.
    apa, ama = _plus_minus(d, a, 4 * μ * (ψ + χ))
    apb, amb = _plus_minus(d, b, 4 * λ * (ψ + χ))
    G(t) = (e = exp(d * (t - tf)); (ama + apa * e) / (amb + apb * e))
    logH(t) = (w = d * (t - tf) / 2; 2w + log(4) - 2 * log((amb + apb * exp(2w)) / d))
    sum(log(n0 - k) for k in 0:nrt-1; init = 0.0) +
        (n0 - nrt) * log(G(t0)) + nrt * logH(t0) +
        (nin > 0 ? nin * log(ψ) : 0.0) +
        sum(log(2 * λ) + logH(t) for t in tbr; init = 0.0) +
        sum(log(ψ * G(t) + χ) - logH(t) for t in ttp; init = 0.0)
end

## (d + x, d - x) for d = sqrt(x² + q), q ≥ 0, each formed without subtracting nearly equal numbers.
_plus_minus(d, x, q) = x > 0 ? (d + x, q / (d + x)) : (q / (d - x), d - x)
