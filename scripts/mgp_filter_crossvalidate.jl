## Generic filter (mgp_filter_pomp) against R phylopomp's filters on the same genealogies.
##
##   Rscript scripts/mgp_filter_crossvalidate.R <model> <ntrees> <Np> <nrep> out.tsv
##   julia --project=test scripts/mgp_filter_crossvalidate.jl out.tsv [Np] [nrep]
##
## For each tree: checks that Julia parses the same node times as R, then prints R's and Julia's
## logmeanexp ± SE, their difference in combined SEs, and for LBDP the exact log-likelihood from R
## and from the Julia port of lbdp_exact. Not part of test/runtests.jl.
using PhyloPOMP
using Random: seed!
using Printf

file = ARGS[1]
Np = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2000
nrep = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20
seed!(parse(Int, get(ENV, "SEED", "20261007")))

## Parameters as in the R script; x0 as R's rinit (nearbyint of pop·f/Σf).
setup = Dict(
    "sir"  => (PhyloPOMP.SIR,  (β = 3.0, γ = 1.0, ψ = 0.3, N = 100.0), (S = 99, I = 1, R = 0)),
    "si2r" => (PhyloPOMP.SI2R, (β = 4.0, κ = 3.0, γ = 1.0, ω = 0.5, ψ = 0.3, η_L = 0.5, η_H = 1.0,
                                N = 100.0, r = 1.0), (S = 99, I_L = 1, I_H = 0, R = 0)),
    "lbdp" => (PhyloPOMP.LBDP, (λ = 1.5, μ = 0.5, ψ = 0.3, χ = 0.2), (n = 1,)),
    "bdei" => (PhyloPOMP.BDEI, (σ = 1.0, λ = 2.0, μ = 0.5, χ = 0.4), (E = 0, I = 1)),
    "bdss" => (PhyloPOMP.BDSS, (λ_nn = 1.0, λ_ns = 0.3, λ_sn = 1.5, λ_ss = 2.5, μ = 0.5, χ = 0.4),
               (N = 1, S = 0)),
)

@printf("%-5s %4s %5s  %18s  %18s  %7s  %10s  %10s\n", "model", "tree", "tips", "R", "Julia", "z", "exact R", "exact jl")
for line in eachline(file)
    f = split(line, '\t')
    model, k, s, T = f[1], f[2], f[3], parse(Float64, f[4])
    rtimes = parse.(Float64, split(f[5], ','))
    rest, rse, rex = parse(Float64, f[6]), parse(Float64, f[7]), f[8]
    M, θ, x0 = setup[model]
    g = parse_newick(String(s); t0 = 0.0, time = T)
    jtimes = sort([[g[i].slate for i in eachindex(g)]; g.time])
    (length(jtimes) == length(rtimes) && maximum(abs.(jtimes .- rtimes)) < 1e-9) ||
        error("tree $k: Julia node times differ from R's")
    g[1].type == PhyloPOMP.Root && g[1].slate == timezero(g) || error("tree $k: root is not at t0")
    p = mgp_filter_pomp(g, M; θ = θ, x0 = x0)
    ll = [logLik(pfilter(p, Np = Np)) for _ in 1:nrep]
    jest, jse = logmeanexp(ll, se = true)
    z = (jest - rest) / sqrt(jse^2 + rse^2)
    exr = rex == "NA" ? "" : @sprintf("%.6f", parse(Float64, rex))
    exj = model == "lbdp" ? @sprintf("%.6f", lbdp_exact(g; θ...)) : ""
    @printf("%-5s %4s %5d  %9.4f ± %6.4f  %9.4f ± %6.4f  %7.2f  %10s  %10s\n",
            model, k, nsample(g), rest, rse, jest, jse, z, exr, exj)
end
