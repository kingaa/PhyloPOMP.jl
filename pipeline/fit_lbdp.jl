# M3 check on LBDP: maximum likelihood by mif on the generated filter (`mgp_filter_pomp`) against the exact MLE,
# found by Nelder–Mead on `lbdp_exact`. Coordinates as the deep-learning nets use them: p_sample known and fixed,
# estimate R0 and the infectious period (λ = R0/IP, μ = (1 − p)/IP, χ = p/IP, ψ = 0).
#
#   julia -t 16 --project=pipeline pipeline/fit_lbdp.jl <trees.tsv> [ntrees] [quick]
#
# The genealogy ends at its last sample (parse_newick's default); both likelihoods use the same genealogy.
# Writes <dir of trees.tsv>/fit_lbdp.csv.

using PhyloPOMP, Optim, Printf, CSV, DataFrames
import PartiallyObservedMarkovProcesses as POMP
using Random: seed!

θ_of(R0, IP, p) = (λ = R0 / IP, μ = (1 - p) / IP, ψ = 0.0, χ = p / IP)
R0_of(θ) = θ.λ / (θ.μ + θ.χ)
IP_of(θ) = 1 / (θ.μ + θ.χ)

## μ and χ move by the same factor, so μ/χ = (1 − p)/p, hence p_sample, stays fixed; mif's geometric mean
## keeps the ratio too.
const SD = 0.02
perturb(scale, lag; λ, μ, χ, _...) = begin
    f = exp(scale * SD * randn())
    (λ = λ * exp(scale * SD * randn()), μ = μ * f, χ = χ * f)
end

exact_mle(g, p) = begin
    ## NaN as +Inf: a NaN minimum would otherwise win the argmin below
    nll(x) = (v = lbdp_exact(g; θ_of(exp(x[1]), exp(x[2]), p)...); isnan(v) ? Inf : -v)
    fits = [optimize(nll, x0, NelderMead(), Optim.Options(g_tol = 1e-12, iterations = 20_000))
            for x0 in ([log(2.0), log(5.0)], [log(1.2), log(2.0)], [log(4.0), log(8.0)])]
    best = fits[argmin(Optim.minimum.(fits))]
    x = Optim.minimizer(best)
    (R0 = exp(x[1]), IP = exp(x[2]), ll = -Optim.minimum(best))
end

main(args) = begin
    tsv = args[1]
    ntrees = length(args) ≥ 2 ? parse(Int, args[2]) : 4
    quick = "quick" in args
    Nmif, Np, Np_eval, nreps, nstart = quick ? (5, 200, 300, 2, 2) : (60, 1000, 2000, 10, 8)
    rows = [split(l, '\t') for l in Iterators.drop(eachline(tsv), 1)][1:ntrees]
    seed!(20261008)
    out = DataFrame()
    for r in rows
        g, p = parse_newick(String(r[7])), parse(Float64, r[4])
        truth = (R0 = parse(Float64, r[2]), IP = parse(Float64, r[3]))
        ex = exact_mle(g, p)
        G = mgp_filter_pomp(g, PhyloPOMP.LBDP; θ = θ_of(ex.R0, ex.IP, p), x0 = (n = 1,))
        ## the filter at the exact MLE: unbiased for exp(ex.ll)
        atmle = pfilter_loglik(G; Np = Np_eval, nreps)
        ## LBDP has no population cap: the filter simulates every hidden host, about exp(r·span) of them, so a start
        ## with fast growth r = (R0 − 1)/IP never finishes (tree 1: r = 2 does not end in minutes, r = 0.41 takes
        ## 0.03 s at Np = 50). Keep starts with r·span ≤ log(tips) + 3, i.e. at most about 20 hidden hosts per tip.
        span = g.time - timezero(g)
        rmax = (log(nsample(g)) + 3) / span
        box = runif_design((R0 = 1.0, IP = 1.0), (R0 = 5.0, IP = 10.0), 100nstart)
        starts = first(filter(s -> (s.R0 - 1) / s.IP ≤ rmax, box), nstart)
        nrow(starts) == nstart || error("only $(nrow(starts)) starts with r ≤ $rmax")
        design = DataFrame([θ_of(s.R0, s.IP, p) for s in eachrow(starts)])
        t = @elapsed fit = POMP.profile(G, design; profiled = (), Nmif, Np, Np_eval, nreps,
                                        perturbations = perturb, cooling = geometric_cooling(0.5))
        b = fit[argmax(fit.loglik), :]
        θb = (λ = b.λ, μ = b.μ, ψ = b.ψ, χ = b.χ)
        row = (tree = parse(Int, r[1]), ntips = nsample(g), p_sample = p,
               R0_true = truth.R0, IP_true = truth.IP,
               R0_exact = ex.R0, IP_exact = ex.IP, ll_exact = ex.ll,
               ll_filter_at_exact = atmle.loglik, se_at_exact = atmle.se,
               R0_mif = R0_of(θb), IP_mif = IP_of(θb), p_mif = θb.χ / (θb.μ + θb.χ),
               ll_mif = b.loglik, se_mif = b.se, ll_exact_at_mif = lbdp_exact(g; θb...), seconds = t)
        push!(out, row)
        @printf("tree %d (%d tips): exact R0 %.3f IP %.3f ll %.3f | filter at exact %.3f ± %.3f | mif R0 %.3f IP %.3f, exact ll there %.3f (%.0f s)\n",
                row.tree, row.ntips, ex.R0, ex.IP, ex.ll, atmle.loglik, atmle.se, row.R0_mif, row.IP_mif,
                row.ll_exact_at_mif, t)
        flush(stdout)
        CSV.write(joinpath(dirname(tsv), "fit_lbdp.csv"), out)    # after every tree, so a stopped run keeps its rows
    end
    CSV.write(joinpath(dirname(tsv), "fit_lbdp.csv"), out)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
