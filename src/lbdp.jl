using SpecialFunctions: logfactorial, logabsbinomial

"""
    lbdp_exact(g; λ, μ, ψ, χ, n0 = missing)

Computes the exact likelihood of genealogy `g` under the linear
birth-death process with the given parameters.  Here `λ` is the birth
rate; `μ`, the death rate; `ψ`, the non-destructive sampling rate;
`χ`, the destructive sampling rate, and `n0` is the initial population
size. If `n0 == missing`, the likelihood is maximized over this
parameter.  The function returns the log-likelihood and the
maximum-likelihood estimate of `n0`.
"""
lbdp_exact(
    g::Genealogy;
    λ, μ, ψ, χ, n0 = missing,
) = begin
    t = [timezero(g), times(g)...]
    t0 = t[begin]
    tf = t[end]
    len = length(t)-1
    tbr = getindex(t,nodes(g))
    ttp = getindex(t,tips(g))
    nrt = length(roots(g))
    nin = length(samples(g))-length(ttp)
    @assert length(ttp) == length(tbr)+nrt "internal inconsistency in genealogy"
    a = λ-μ+ψ+χ
    b = λ-μ-ψ-χ
    d = sqrt(b*b+4*λ*(ψ+χ)) ## guaranteed to be real
    G = function(t)
        omega = d*(t-tf)/2
        C = d*cosh(omega)
        S = sinh(omega)
        (C+a*S)/(C+b*S)
    end
    H = function(t)
        omega = d*(t-tf)/2
        g = cosh(omega)+b/d*sinh(omega)
        1/g/g
    end
    nhat= Int64(floor(nrt/(1-G(t0))))
    if ismissing(n0)
        n0 = nhat
    else
        n0 = Int64(n0)
    end
    @assert n0 ≥ 0 "`n0` must be a nonnegative integer."
    ll = first(logabsbinomial(n0,nrt))+
        logfactorial(nrt)+
        (n0-nrt)*log(G(t0))+
        nrt*log(H(t0))+
        (nin > 0 ? nin*sum(log(ψ)) : 0)+
        sum(log.(2λ.*H.(tbr)))+
        sum(log.((ψ.*G.(ttp).+χ)./H.(ttp)))
    ll, nhat
end
