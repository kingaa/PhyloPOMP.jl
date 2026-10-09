# BD maximum likelihood with and without conditioning on "the epidemic reached t samples", which is how
# generate_trees.jl stop=tips makes its data (a run that dies out first is rerun). Conditioned likelihood:
# lbdp_exact / P(reach t | R0, p), P from the closed form below (checked against a power series, a direct Catalan sum and
# two simulations, 2026-10-09).
#
#   julia -t 4 --project=pipeline scripts/bd_conditioning.jl <run_dir> [ntest = 500]
#
# Uses the last ntest trees of <run_dir>/trees.tsv (compare.jl's test split when N = 10,000) and writes
# <run_dir>/conditioning.tsv. Result on output/val10k (500 trees): R0 MAPE 15.60 → 15.27, bias +6.76% → +4.36%,
# 95% profile coverage 93.8% → 95.4%; IP MAPE 13.27 → 13.04, bias +5.14% → +3.35%, coverage 93.8% → 93.2%.
#
using PhyloPOMP, Optim, Printf, Statistics, Random
include(joinpath(@__DIR__, "..", "pipeline", "fit_lbdp.jl"))    # θ_of, exact_mle

function log_reach(t, R0, p)
    q = R0 / (R0 + 1); c = 4q * (1 - q); A = 1 - c * (1 - p); B = c * p / A
    coef = (1 - sqrt(A)) / (2q)
    S = coef
    if t > 1
        coef = sqrt(A) * B / (4q); S += coef
        for k in 2:t-1
            coef *= B * (k - 1.5) / k; S += coef
        end
    end
    P = 1 - S
    P > 1e-8 && return log(P)
    ## small: the tail sum Σ_{k≥t} plus survival (subcritical here, B well below 1)
    tail = 0.0; k = t
    coef *= B * (k - 1.5) / k
    while coef > 1e-18 * max(tail, 1e-300) && k < 10^7
        tail += coef; k += 1; coef *= B * (k - 1.5) / k
    end
    log(max(0.0, 1 - 1 / R0) + tail)
end

## Monte Carlo check of log_reach on the embedded jump chain
mc_reach(t, R0, p, n, rng) = begin
    q = R0 / (R0 + 1); hits = 0
    for _ in 1:n
        alive, s = 1, 0
        while alive > 0 && s < t && alive < 20_000
            if rand(rng) < q
                alive += 1
            else
                alive -= 1; rand(rng) < p && (s += 1)
            end
        end
        hits += (s >= t || alive >= 20_000)
    end
    hits / n
end

const SEARCH = (R0 = (0.01, 100.0), IP = (0.01, 1000.0)); const CUT = 3.841458820694124
profile(ll, which, v) = begin
    lo, hi = log.(SEARCH[which == :R0 ? :IP : :R0])
    f(x) = -(which == :R0 ? ll(v, exp(x)) : ll(exp(x), v))
    -Optim.minimum(optimize(f, lo, hi, Brent()))
end
endpoint(ll, which, est, ℓ̂, bound) = begin
    gap(v) = 2 * (ℓ̂ - profile(ll, which, v)) - CUT
    gap(bound) < 0 && return bound
    a, b = log(est), log(bound)
    for _ in 1:40
        m = (a + b) / 2; gap(exp(m)) < 0 ? (a = m) : (b = m)
    end
    exp((a + b) / 2)
end
fit(ll) = begin
    nll(x) = -ll(exp(x[1]), exp(x[2]))
    fits = [optimize(nll, x0, NelderMead(), Optim.Options(g_tol = 1e-12, iterations = 20_000))
            for x0 in ([log(2.0), log(5.0)], [log(1.2), log(2.0)], [log(4.0), log(8.0)])]
    best = fits[argmin(Optim.minimum.(fits))]; x = Optim.minimizer(best); ℓ̂ = -Optim.minimum(best)
    R0, IP = exp(x[1]), exp(x[2])
    (R0 = R0, IP = IP, ll = ℓ̂,
     R0_lo = endpoint(ll, :R0, R0, ℓ̂, SEARCH.R0[1]), R0_hi = endpoint(ll, :R0, R0, ℓ̂, SEARCH.R0[2]),
     IP_lo = endpoint(ll, :IP, IP, ℓ̂, SEARCH.IP[1]), IP_hi = endpoint(ll, :IP, IP, ℓ̂, SEARCH.IP[2]))
end

function main(args)
    rng = Xoshiro(1)
    println("P(reach t): closed form vs simulation (200,000 runs)")
    for (t, R0, p) in ((50, 1.05, 0.5), (50, 1.5, 0.05), (120, 0.9, 0.9), (50, 3.0, 0.3), (199, 1.2, 1.0))
        @printf("  t %3d R0 %.2f p %.2f: %.5f vs %.5f   (1 − 1/R0 = %.4f)\n", t, R0, p, exp(log_reach(t, R0, p)),
                mc_reach(t, R0, p, 200_000, rng), max(0, 1 - 1 / R0))
    end
    flush(stdout)
    dir = args[1]
    ntest = length(args) >= 2 ? parse(Int, args[2]) : 500
    lines = collect(eachline(joinpath(dir, "trees.tsv")))
    hdr = split(lines[1], '\t'); col(n) = findfirst(==(n), hdr)
    rows = [split(l, '\t') for l in lines[end-ntest+1:end]]   # compare.jl's test split when N = 10,000: the last 500 trees
    n = length(rows)
    res = Vector{Any}(undef, n)
    Threads.@threads :dynamic for k in 1:n
        r = rows[k]
        g = parse_newick(String(r[col("newick")])); p = parse(Float64, r[col("p_sample")])
        t = nsample(g)
        llu(R0, IP) = lbdp_exact(g; θ_of(R0, IP, p)...)
        llc(R0, IP) = llu(R0, IP) - log_reach(t, R0, p)
        res[k] = (R0 = parse(Float64, r[col("R_0")]), IP = parse(Float64, r[col("Infectious_Period")]), p = p, t = t,
                  u = fit(llu), c = fit(llc))
    end
    clip(x, lo, hi) = clamp(x, lo, hi)
    for (name, key) in (("unconditioned", :u), ("conditioned", :c))
        for (par, box) in ((:R0, (1.0, 5.0)), (:IP, (1.0, 10.0)))
            tr = [getfield(r, par) for r in res]; est = [getfield(getfield(r, key), par) for r in res]
            re = (est .- tr) ./ tr; rec = (clip.(est, box...) .- tr) ./ tr
            lo = [getfield(getfield(r, key), Symbol(par, :_lo)) for r in res]
            hi = [getfield(getfield(r, key), Symbol(par, :_hi)) for r in res]
            cov = mean(lo .<= tr .<= hi)
            @printf("%-13s %-2s  MAPE %5.2f  bias %+6.2f  | clipped MAPE %5.2f bias %+6.2f | coverage95 %.1f%%\n",
                    name, par, 100mean(abs.(re)), 100mean(re), 100mean(abs.(rec)), 100mean(rec), 100cov)
        end
    end
    ## where the two differ: by true R0
    for (lo, hi) in ((1.0, 1.5), (1.5, 2.5), (2.5, 5.0))
        sel = [lo <= r.R0 < hi for r in res]
        reu = [(r.u.R0 - r.R0) / r.R0 for r in res[sel]]; rec = [(r.c.R0 - r.R0) / r.R0 for r in res[sel]]
        @printf("true R0 in [%.1f, %.1f): %3d trees  R0 bias unconditioned %+6.2f%%  conditioned %+6.2f%%  MAPE %5.2f → %5.2f\n",
                lo, hi, sum(sel), 100mean(reu), 100mean(rec), 100mean(abs.(reu)), 100mean(abs.(rec)))
    end
    open(joinpath(dir, "conditioning.tsv"), "w") do io
        println(io, "R0\tIP\tp\tt\tR0_u\tIP_u\tll_u\tR0_c\tIP_c\tll_c")
        for r in res
            println(io, join((r.R0, r.IP, r.p, r.t, r.u.R0, r.u.IP, r.u.ll, r.c.R0, r.c.IP, r.c.ll), '\t'))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
