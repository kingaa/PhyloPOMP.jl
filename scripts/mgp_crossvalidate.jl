# Compare src/simulate.jl with R phylopomp's runSEIR, runMERS, runSIR, runSI2R
# and runMTBD2.
# Run scripts/mgp_crossvalidate.R first. The parameters match that script.
# Prints a two-sample KS test for each statistic. Writes nothing.
#
# Per tree: nsample, internal nodes, inline (degree-1) samples, first/last/mean
# sample time, first/mean internal-node time, total branch length. With a
# fourth file (MERS, SI2R, MTBD) it also compares sample counts per deme.
# Final compartment counts come from <rfile>.states when that file exists.
#
# The si2r case runs PhyloPOMP.SI2R with r = 1, so every sample removes the
# host, as in R's runSI2R.
#
# Baseline, 2026-10-03. SEIR N=1000, SEIR with χ=0.1 ("seirchi") N=4000,
# MERS N=2000. Of 37 statistics, two had p < 0.05: SEIR final S (0.029) and
# seirchi first-sample time (0.034). An independent R seed and six Julia
# seeds did not repeat either one (first-sample time, pooled z = -1.48 on
# 18,512 vs 6,196 trees). Seeds are fixed, so a change that draws the same
# random numbers must reprint the Julia column exactly. One p < 0.05 per
# 30 tests is about what chance gives. A real bug fails the same statistic
# on more than one seed. SEED=n in the environment of both scripts gives a
# fresh pair.
#
# Baseline, 2026-10-05, N=2000. SIR: none of 11 statistics with p < 0.05.
# MTBD: none of 13. SI2R (r = 1): mean internal-node time (0.029) and total
# branch length (0.017) of 14. With SEED=7 and N=6000 those two gave 0.86 and
# 0.62, and none of the 14 was below 0.05.
using PhyloPOMP
using PhyloPOMP: Root, Node, Sample
using PhyloPOMP.SoftMERS.Demes: Camel, Human
using Random: MersenneTwister
using Printf

length(ARGS) ≥ 3 || error("usage: mgp_crossvalidate.jl seir|seirchi|mers|sir|si2r|mtbd N rfile [rfile_unobscured]")
model = ARGS[1]; N = parse(Int, ARGS[2]); rfile = ARGS[3]
rfile2 = length(ARGS) ≥ 4 ? ARGS[4] : nothing

if model == "seir" || model == "seirchi"
    ψ, χ = model == "seir" ? (0.30, 0.0) : (0.20, 0.10)   # seirchi: destructive `culling` on
    θ = (β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = ψ, χ = χ, N = 100.0)
    x0 = (S = 99, E = 0, I = 1, R = 0); graft = [0, 1]; tmax = 20.0
    M = PhyloPOMP.SEIR; D = PhyloPOMP.Unstructured; smap = nothing
elseif model == "mers"
    θ = (β_cc = 3.0, β_ch = 0.5, β_hc = 0.5, β_hh = 3.0, γ_c = 1.0, γ_h = 1.0,
         χ_c = 0.3, χ_h = 0.3, B_c = 0.5, B_h = 0.5, N_c = 20.0, N_h = 20.0)
    x0 = (S_c = 19, I_c = 1, S_h = 20, I_h = 0); graft = [1, 0]; tmax = 10.0
    M = PhyloPOMP.MERS; D = PhyloPOMP.SoftMERS.Demes; smap = [Camel, Human]
elseif model == "sir"
    θ = (β = 4.0, γ = 1.0, ψ = 0.3, N = 100.0)
    x0 = (S = 99, I = 1, R = 0); graft = [1]; tmax = 20.0
    M = PhyloPOMP.SIR; D = PhyloPOMP.Unstructured; smap = nothing
elseif model == "si2r"
    θ = (β = 4.0, κ = 3.0, γ = 1.0, ω = 0.5, ψ = 0.3, η_L = 0.5, η_H = 1.0, N = 100.0, r = 1.0)
    x0 = (S = 99, I_L = 1, I_H = 0, R = 0); graft = [1, 0]; tmax = 20.0
    M = PhyloPOMP.SI2R; D = PhyloPOMP.SI2RDemes; smap = [D.I_L, D.I_H]
elseif model == "mtbd"
    θ = (lambda11 = 1.2, lambda12 = 0.3, lambda21 = 0.2, lambda22 = 0.9, m12 = 0.2, m21 = 0.1,
         mu1 = 0.5, mu2 = 0.5, psi1 = 0.3, psi2 = 0.3, r1 = 0.7, r2 = 1.0)
    x0 = (I1 = 1, I2 = 0); graft = [1, 0]; tmax = 6.0
    M = PhyloPOMP.MTBD; D = PhyloPOMP.MTBDDemes; smap = [D.I1, D.I2]
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
     tbl = isempty(s) ? missing : tbl)
end
## samples per deme, in `smap` order; R labels demes 1, 2, ... in the same order
demecount(g) = [count(g[i].deme === d for i in samples(g)) for d in smap]

## SEED in the environment overrides the default seed, as in the R script
seed = parse(Int, get(ENV, "SEED", "20261003"))
rng = MersenneTwister(seed)
## `simulate_trajectory` draws exactly what `simulate` draws, so the trees are
## the ones `simulate` would give; the final state is a free extra.
traj = [simulate_trajectory(M, θ; x0, graft, tmax, rng, demeset = D, samplemap = smap) for _ in 1:N]
J = [stats(tr.genealogy) for tr in traj]
Jdeme = isnothing(smap) ? nothing : [demecount(tr.genealogy) for tr in traj]
Jfinal = [tr.states[end] for tr in traj]

lines = readlines(rfile)
length(lines) == N || error("$rfile has $(length(lines)) lines, expected $N")
R = map(lines) do s
    s = strip(s)
    isempty(s) && return stats(Genealogy{D}(0.0, tmax))
    stats(parse_newick(String(s); demes = D, t0 = 0.0, time = tmax))
end
Rsp = isnothing(rfile2) ? nothing : map(readlines(rfile2)) do s
    s = strip(s)
    isempty(s) && return zeros(Int, length(smap))
    demecount(parse_newick(String(s); demes = D, t0 = 0.0, time = tmax))
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

println("model=$model  N=$N per side  (Julia seed $seed)")
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
## final population state (needs <rfile>.states from mgp_crossvalidate.R);
## R's yaml names SEIR compartments S E I R and MERS ones Sc Ic Sh Ih
statefile = rfile * ".states"
if isfile(statefile)
    Rfinal = map(readlines(statefile)) do s
        Dict(String(k) => parse(Float64, v) for (k, v) in (split(kv, "=") for kv in split(strip(s))))
    end
    length(Rfinal) == N || error("$statefile has $(length(Rfinal)) lines, expected $N")
    rname(c) = replace(String(c), "_" => "")
    println("final state at tmax (Julia via simulate_trajectory, R via yaml):")
    for c in keys(x0)
        a = Float64[getfield(x, c) for x in Jfinal]
        b = Float64[r[rname(c)] for r in Rfinal]
        if all(==(a[1]), a) && all(==(a[1]), b)
            @printf("%-22s %10.3f %10.3f  (constant on both sides)\n", "final $c", mean(a), mean(b)); continue
        end
        Dst, p = ks_two_sample(a, b)
        @printf("%-22s %10.3f %10.3f %8.4f %8.4f%s\n", "final $c", mean(a), mean(b), Dst, p,
                p < 0.05 ? "  <-- significant" : "")
    end
end
if !isnothing(Rsp)
    for (k, d) in enumerate(smap)
        name = "n $(Symbol(d)) samples"
        a = Float64[j[k] for j in Jdeme]; b = Float64[r[k] for r in Rsp]
        Dst, p = ks_two_sample(a, b)
        @printf("%-22s %10.3f %10.3f %8.4f %8.4f%s\n", name, mean(a), mean(b), Dst, p,
                p < 0.05 ? "  <-- significant" : "")
    end
end
