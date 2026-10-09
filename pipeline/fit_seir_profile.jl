# M3 on SEIR: maximum likelihood and a profile confidence interval for β from one simulated
# genealogy, with the filter generated from the `@mgp` SEIR table (`mgp_filter_pomp`).
#
#   julia -t 16 --project=pipeline pipeline/fit_seir_profile.jl          # full run
#   julia -t 16 --project=pipeline pipeline/fit_seir_profile.jl quick    # smoke test, about 1 min
#
# Estimated: β and ψ. Held at their true values: σ, γ, ω, χ, N, and the initial state.
# Writes output/m3_seir/{search,profile,mcap}.csv.

using PhyloPOMP
import PartiallyObservedMarkovProcesses as POMP
using CSV, DataFrames
using Random: MersenneTwister, seed!

quick = "quick" in ARGS
outdir = joinpath(@__DIR__, "..", "output", "m3_seir")
mkpath(outdir)

## The tree: the diagnostics tutorial's (33 tips).
θ_true = (β=4.0, σ=1.0, γ=1.0, ω=1.0, ψ=0.10, χ=0.0, N=100.0)
x0 = (S=99, E=0, I=1, R=0)
rng = MersenneTwister(20260813)
g = let g = nothing
    for _ in 1:100
        g = simulate(PhyloPOMP.SEIR, θ_true; x0=x0, graft=[0,1], tmax=20.0, rng=rng)
        nsample(g) ≥ 20 && break
    end
    g
end
G = mgp_filter_pomp(g, PhyloPOMP.SEIR; θ=θ_true, x0=x0)
println("tree: $(nsample(g)) samples; threads: $(Threads.nthreads())")

## Run sizes. Np_eval and nreps set the standard error of each evaluated point.
Nmif, Np, Np_eval, nreps = quick ? (5, 200, 300, 2) : (50, 1000, 2000, 10)
nstart, ngrid, nprof = quick ? (2, 6, 1) : (12, 12, 3)   # mcap needs 0.75·ngrid ≥ 4
cool = geometric_cooling(0.5)
seed!(20261008)

## 1. Global search over (β, ψ): mif from random starts, then replicate pfilter at each end point.
starts = runif_design((β=1.0, ψ=0.02), (β=10.0, ψ=0.5), nstart)
t1 = @elapsed search = POMP.profile(
    G, starts;
    profiled=(), Nmif, Np, Np_eval, nreps,
    perturbations=@perturbn(β ~ LogNormal(0.02), ψ ~ LogNormal(0.02)),
    cooling=cool,
)
sort!(search, :loglik, rev=true)
CSV.write(joinpath(outdir, "search.csv"), search)
println("search ($(round(t1/60, digits=1)) min):")
show(stdout, search[:, [:β, :ψ, :loglik, :se, :ess]]; allrows=true); println(); flush(stdout)

## 2. Profile over β: ψ re-fitted by mif at each grid value, from nprof starts.
pd = profile_design(
    β=range(2.0, 8.0, length=ngrid);
    lower=(ψ=0.03,), upper=(ψ=0.3,), nprof=nprof,
)
t2 = @elapsed prof = POMP.profile(
    G, pd;
    Nmif, Np, Np_eval, nreps,
    perturbations=@perturbn(ψ ~ LogNormal(0.02)),
    cooling=cool,
)
CSV.write(joinpath(outdir, "profile.csv"), prof)
println("profile ($(round(t2/60, digits=1)) min)"); flush(stdout)

## 3. mcap on the best point per grid value.
best = combine(groupby(filter(:loglik => isfinite, prof), :β), sdf -> sdf[argmax(sdf.loglik), :])
show(stdout, best[:, [:β, :ψ, :loglik, :se]]; allrows=true); println()
m = mcap(best.loglik, best.β)
CSV.write(joinpath(outdir, "mcap.csv"),
          DataFrame(mle=m.mle, lo=m.ci[1], hi=m.ci[2], se=m.se, se_stat=m.se_stat, se_mc=m.se_mc,
                    truth=θ_true.β, covers=m.ci[1] ≤ θ_true.β ≤ m.ci[2]))
println(m)
println("true β = $(θ_true.β): ",
        any(isnan, m.ci) ? "interval undetermined" :
        m.ci[1] ≤ θ_true.β ≤ m.ci[2] ? "inside the 95% interval" : "outside the 95% interval")
