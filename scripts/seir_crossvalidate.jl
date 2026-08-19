# Milestone 3: distribution-level cross-validation of the Julia forward
# simulator (src/simulate.jl) against R phylopomp's runSEIR.
#
# Run scripts/seir_crossvalidate.R first to produce the R-side Newick
# trees, then this script:
#   1. generates N independent SEIR draws on the Julia side directly via
#      `simulate`,
#   2. parses the N R-side Newick trees written by the R script (through
#      the SAME `parse_newick`/`nsample`/`samples` code used for Julia's
#      own trees, so both sides are measured by identical code -- the only
#      thing R contributes is the raw simulation, not any of the analysis),
#   3. compares `nsample` (primary: the statistic most sensitive to the
#      event-rate/sampling-process structure being validated) and
#      first-sample time conditional on nsample>=1 (secondary: probes the
#      waiting-time/hazard-normalization machinery, which a size statistic
#      does not) via a hand-rolled two-sample Kolmogorov-Smirnov test,
#   4. reports the extinction fraction on each side SEPARATELY (not folded
#      into the KS test): with ~30% of mass at nsample=0 on both sides, the
#      ECDFs agree exactly there, which deflates D and makes a KS pass
#      weaker evidence than the p-value alone would suggest.
#
# This script deliberately lives outside test/runtests.jl: it needs a
# working R + phylopomp installation, takes minutes to run, and a
# stochastic KS test has no place in CI (see check_milestone2.md's
# precedent for documenting such placement decisions explicitly).
#
# Usage:
#   Rscript scripts/seir_crossvalidate.R 1000 seir_crossvalidate_r_trees.txt
#   julia --project=. scripts/seir_crossvalidate.jl 1000 seir_crossvalidate_r_trees.txt

using PhyloPOMP
using Random: MersenneTwister

const N = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 1000
const RFILE = length(ARGS) ≥ 2 ? ARGS[2] : "seir_crossvalidate_r_trees.txt"

## Parameters MATCHED EXACTLY to scripts/seir_crossvalidate.R -- see that
## file's header comment for why E0=0 (single founding lineage, avoiding
## the untested multi-root path) is a deliberate choice, not an oversight.
θ = (β=4.0, σ=1.0, γ=1.0, ω=1.0, ψ=0.30, χ=0.0, N=100.0)
x0 = (S=99, E=0, I=1, R=0)
graft = [0,1]
tmax = 20.0

## --- Julia draws ---------------------------------------------------------
rng = MersenneTwister(20260813)
j_nsample = Vector{Int}(undef, N)
j_first = Vector{Union{Missing,Float64}}(undef, N)
for i ∈ 1:N
    g = simulate(PhyloPOMP.SEIR, θ; x0=x0, graft=graft, tmax=tmax, rng=rng)
    j_nsample[i] = nsample(g)
    j_first[i] = nsample(g) > 0 ? minimum(g[k].slate for k ∈ samples(g)) : missing
end

## --- R draws (parsed through the Julia analysis code) --------------------
isfile(RFILE) || error(
    "$RFILE not found -- run `Rscript scripts/seir_crossvalidate.R $N $RFILE` first"
)
lines = readlines(RFILE)
length(lines) == N || error(
    "expected $N lines from $RFILE, got $(length(lines)) -- N must match "*
    "between the R and Julia invocations"
)
r_nsample = Vector{Int}(undef, N)
r_first = Vector{Union{Missing,Float64}}(undef, N)
for i ∈ 1:N
    s = strip(lines[i])
    if isempty(s)
        r_nsample[i] = 0
        r_first[i] = missing
    else
        g = parse_newick(String(s), t0=0.0)
        r_nsample[i] = nsample(g)
        r_first[i] = minimum(g[k].slate for k ∈ samples(g))
    end
end

## --- extinction fractions (reported separately, not via KS -- see header) -
j_ext, r_ext = count(==(0),j_nsample)/N, count(==(0),r_nsample)/N
## normal-approximation CI for a proportion from N draws.
propci(p,n) = (se = sqrt(p*(1-p)/n); (max(0.0,p-1.96se), min(1.0,p+1.96se)))
j_ext_ci, r_ext_ci = propci(j_ext,N), propci(r_ext,N)
## two-proportion z-test for a difference (pooled SE).
pdiff = j_ext - r_ext
ppool = (count(==(0),j_nsample)+count(==(0),r_nsample))/(2N)
se_pool = sqrt(ppool*(1-ppool)*(2/N))
z_ext = pdiff/se_pool

## --- hand-rolled two-sample KS test (asymptotic p-value) ------------------
## No new dependency (HypothesisTests) added for a one-off validation
## script -- this is the standard textbook formula (Marsaglia/Kolmogorov),
## accurate at n,m = O(1000).
function ks_two_sample(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n, m = length(x), length(y)
    xs, ys = sort(x), sort(y)
    pts = sort(unique(vcat(xs,ys)))
    D = 0.0
    for v ∈ pts
        Fx = searchsortedlast(xs,v)/n
        Fy = searchsortedlast(ys,v)/m
        D = max(D, abs(Fx-Fy))
    end
    λ = D*sqrt(n*m/(n+m))
    p = clamp(2*sum((-1)^(k-1)*exp(-2*k^2*λ^2) for k ∈ 1:100), 0.0, 1.0)
    D, p
end

D_n, p_n = ks_two_sample(Float64.(j_nsample), Float64.(r_nsample))

j_first_c = collect(skipmissing(j_first))
r_first_c = collect(skipmissing(r_first))
D_f, p_f = ks_two_sample(j_first_c, r_first_c)

## --- summary stats (no Statistics dependency; a handful of one-liners) ---
mymean(x) = sum(x)/length(x)
function myquantile(x,q)
    xs = sort(x)
    xs[clamp(round(Int, q*(length(xs)-1))+1, 1, length(xs))]
end
qs(x) = (myquantile(x,0.25), myquantile(x,0.5), myquantile(x,0.75))

## --- report ----------------------------------------------------------------
report = """
# SEIR forward-simulator cross-validation: Julia vs R phylopomp

N = $N draws per side. Julia seed: `MersenneTwister(20260813)` (single
continuing stream). R seed: `set.seed(20260813)` (see
`scripts/seir_crossvalidate.R`). Parameters matched exactly on both sides
(β=4, σ=1, γ=1, ψ=0.30, χ=0, ω=1, N=100, S0=0.99, E0=0, I0=0.01, R0=0,
t0=0, tmax=20); see that script's header comment for why E0=0 (single
founding lineage) was chosen deliberately.

## Extinction fraction (reported separately from KS -- see script header)

| side  | fraction | 95% CI |
|-------|----------|--------|
| Julia | $(round(j_ext,digits=3)) | ($(round(j_ext_ci[1],digits=3)), $(round(j_ext_ci[2],digits=3))) |
| R     | $(round(r_ext,digits=3)) | ($(round(r_ext_ci[1],digits=3)), $(round(r_ext_ci[2],digits=3))) |

Two-proportion z-statistic (difference): z = $(round(z_ext,digits=3))
($(abs(z_ext) < 1.96 ? "not significant at α=0.05" : "SIGNIFICANT at α=0.05 -- investigate"))

## nsample distribution (primary statistic)

|       | n   | mean | Q25/median/Q75 |
|-------|-----|------|-----------------|
| Julia | $N | $(round(mymean(j_nsample),digits=1)) | $(qs(j_nsample)) |
| R     | $N | $(round(mymean(r_nsample),digits=1)) | $(qs(r_nsample)) |

Two-sample KS: D = $(round(D_n,digits=4)), p = $(round(p_n,digits=4))
($(p_n > 0.05 ? "no evidence of a distributional difference at α=0.05" : "SIGNIFICANT difference at α=0.05 -- investigate"))

## First-sample-time distribution (secondary statistic, conditional on nsample≥1)

|       | n (nonextinct) | mean | Q25/median/Q75 |
|-------|-----|------|-----------------|
| Julia | $(length(j_first_c)) | $(round(mymean(j_first_c),digits=3)) | $(qs(j_first_c)) |
| R     | $(length(r_first_c)) | $(round(mymean(r_first_c),digits=3)) | $(qs(r_first_c)) |

Two-sample KS: D = $(round(D_f,digits=4)), p = $(round(p_f,digits=4))
($(p_f > 0.05 ? "no evidence of a distributional difference at α=0.05" : "SIGNIFICANT difference at α=0.05 -- investigate"))
"""

println(report)
write("seir_crossvalidate_results.md", report)
println("Report written to seir_crossvalidate_results.md")
