# Cross-validation of the Julia forward simulator (src/simulate.jl) against
# R phylopomp's runSEIR / runMERS, on statistics that depend on internal-node
# TIMES as well as counts (the committed seir_crossvalidate.jl checked counts
# and first-sample time only; see check_milestone4.md for why that is not
# enough).
#
# Run scripts/mgp_crossvalidate.R first (see its header for the commands).
# Parameters here are matched EXACTLY to that script. Nothing is written to
# disk; the report is printed.
#
# Statistics per tree (non-extinct unless noted): nsample, number of internal
# nodes, number of inline (degree-1) samples, first/last/mean sample time,
# first/mean internal-node time, total branch length; for MERS also the
# camel/human sample counts. Two-sample KS (asymptotic p) on each.
#
# Baseline (2026-10-03, commit 120ca6a): SEIR N=1000 all p ≥ 0.08;
# MERS N=2000 all p ≥ 0.13. Fixed seeds on both sides, so a code change that
# does not alter RNG consumption must reproduce the Julia column exactly.
using PhyloPOMP
using PhyloPOMP: Root, Node, Sample
using PhyloPOMP.SoftMERS.Demes: Camel, Human
using Random: MersenneTwister
using Printf

length(ARGS) ≥ 3 || error("usage: mgp_crossvalidate.jl seir|mers N rfile [rfile_unobscured]")
model = ARGS[1]; N = parse(Int, ARGS[2]); rfile = ARGS[3]
rfile2 = length(ARGS) ≥ 4 ? ARGS[4] : nothing

if model == "seir"
    θ = (β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.30, χ = 0.0, N = 100.0)
    x0 = (S = 99, E = 0, I = 1, R = 0); graft = [0, 1]; tmax = 20.0
    M = PhyloPOMP.SEIR; D = PhyloPOMP.Unstructured; smap = nothing
elseif model == "mers"
    θ = (β_cc = 3.0, β_ch = 0.5, β_hc = 0.5, β_hh = 3.0, γ_c = 1.0, γ_h = 1.0,
         χ_c = 0.3, χ_h = 0.3, B_c = 0.5, B_h = 0.5, N_c = 20.0, N_h = 20.0)
    x0 = (S_c = 19, I_c = 1, S_h = 20, I_h = 0); graft = [1, 0]; tmax = 10.0
    M = PhyloPOMP.MERS; D = PhyloPOMP.SoftMERS.Demes; smap = [Camel, Human]
else
    error("unknown model: $model")
end

function stats(g)
    s = samples(g); n = nodes(g)
    st = [g[i].slate for i in s]; nt = [g[i].slate for i in n]
    tbl = sum(g[i].slate - g[g[i].parent].slate for i in eachindex(g) if !isnothing(g[i].parent); init = 0.0)
    (nsample = length(s),
     ninternal = length(n),
     ninline = count(length(g[i].children) == 1 for i in s),
     first_sample = isempty(st) ? missing : minimum(st),
     last_sample = isempty(st) ? missing : maximum(st),
     mean_sample = isempty(st) ? missing : sum(st) / length(st),
     first_internal = isempty(nt) ? missing : minimum(nt),
     mean_internal = isempty(nt) ? missing : sum(nt) / length(nt),
     tbl = isempty(s) ? missing : tbl,
     ncamel = count(g[i].deme === Camel for i in s),
     nhuman = count(g[i].deme === Human for i in s))
end

rng = MersenneTwister(20261003)
J = [stats(simulate(M, θ; x0, graft, tmax, rng, demeset = D, samplemap = smap)) for _ in 1:N]

lines = readlines(rfile)
length(lines) == N || error("$rfile has $(length(lines)) lines, expected $N")
R = map(lines) do s
    s = strip(s)
    isempty(s) && return stats(Genealogy{D}(0.0, tmax))
    stats(parse_newick(String(s); demes = D, t0 = 0.0, time = tmax))
end
Rsp = isnothing(rfile2) ? nothing : map(readlines(rfile2)) do s
    s = strip(s)
    isempty(s) && return (ncamel = 0, nhuman = 0)
    g = parse_newick(String(s); demes = D, t0 = 0.0, time = tmax)
    (ncamel = count(g[i].deme === Camel for i in samples(g)),
     nhuman = count(g[i].deme === Human for i in samples(g)))
end

function ks_two_sample(x, y)
    n, m = length(x), length(y)
    xs, ys = sort(x), sort(y)
    Dst = 0.0
    for v in sort(unique(vcat(xs, ys)))
        Dst = max(Dst, abs(searchsortedlast(xs, v) / n - searchsortedlast(ys, v) / m))
    end
    λ = Dst * sqrt(n * m / (n + m))
    p = clamp(2 * sum((-1)^(k - 1) * exp(-2 * k^2 * λ^2) for k in 1:100), 0.0, 1.0)
    Dst, p
end
mean(x) = sum(x) / length(x)

println("model=$model  N=$N per side  (Julia seed 20261003, R seed 20261003)")
jext = count(j -> j.nsample == 0, J) / N; rext = count(r -> r.nsample == 0, R) / N
pp = (jext + rext) / 2; z = (jext - rext) / sqrt(pp * (1 - pp) * 2 / N)
@printf("extinction fraction: Julia %.3f  R %.3f  (two-proportion z = %.2f)\n", jext, rext, z)
@printf("%-22s %10s %10s %8s %8s\n", "statistic", "Julia mean", "R mean", "KS D", "KS p")
for name in (:nsample, :ninternal, :ninline, :first_sample, :last_sample, :mean_sample,
             :first_internal, :mean_internal, :tbl)
    jx = Float64[getfield(j, name) for j in J if !ismissing(getfield(j, name))]
    rx = Float64[getfield(r, name) for r in R if !ismissing(getfield(r, name))]
    (isempty(jx) || isempty(rx)) && continue
    if all(==(jx[1]), jx) && all(==(jx[1]), rx)
        @printf("%-22s %10.3f %10.3f  (constant on both sides)\n", name, mean(jx), mean(rx)); continue
    end
    Dst, p = ks_two_sample(jx, rx)
    @printf("%-22s %10.3f %10.3f %8.4f %8.4f  (n=%d/%d)%s\n", name, mean(jx), mean(rx), Dst, p,
            length(jx), length(rx), p < 0.05 ? "  <-- significant" : "")
end
if !isnothing(Rsp)
    for (name, a, b) in (("n camel samples", Float64[j.ncamel for j in J], Float64[r.ncamel for r in Rsp]),
                         ("n human samples", Float64[j.nhuman for j in J], Float64[r.nhuman for r in Rsp]))
        Dst, p = ks_two_sample(a, b)
        @printf("%-22s %10.3f %10.3f %8.4f %8.4f%s\n", name, mean(a), mean(b), Dst, p,
                p < 0.05 ? "  <-- significant" : "")
    end
end
