# M7: deep-learning estimates against maximum likelihood on the test trees of a run, with Voznica et al. (2022)'s
# metrics (Methods Extended pp. 34–35).
#
#   julia -t 16 --project=pipeline pipeline/compare.jl <run_dir>
#
# Test trees: the same split as train_dl.jl (the last min(10000, N÷20) rows of cblv.h5).
# Methods:
#   - our networks: <net>_<arch>_test.csv written by train_dl.jl (infectious period already × rescale);
#   - phylodeep's pretrained BD models (reference): BD_SMALL/LARGE_CNN.h5 through `load_keras`, and BD_SMALL_FFNN.h5
#     with its scaler (pipeline/test/data/phylodeep_bd_small_ffnn_scaler.csv), if PHYLODEEP_MODELS (or the venv
#     default) has them;
#   - MLE, p_sample fixed. BD: exact LBDP likelihood (`lbdp_exact`), maximized over (R0, IP) by `exact_mle` of
#     fit_lbdp.jl. BDEI, BDSS: exact likelihood from the multi-type birth-death equations (`mtbd_exact`; BDSS founder
#     a superspreader with probability f, as generated), maximized by `exact_mle_mtbd` of mle_models.jl.
#     The genealogy ends at its last sample, the stopping time of a stop=tips simulation. The likelihood is not
#     conditioned on the run reaching that many samples, although generate_trees.jl restarts runs that die out first.
#     "MLE (clipped)": the estimate clipped to the prior box (R0 ∈ [1, 5], IP ∈ [1, 10], incubation period ∈ [0.2, 50],
#     X ∈ [3, 10], f ∈ [0.05, 0.2]), as Voznica clipped BEAST2 estimates for comparison with the networks.
#     BD only: 95% profile confidence intervals: values of one parameter, the other maximized out, where
#     2(ℓ̂ − ℓ_prof) ≤ 3.84 (χ²₁, 0.95); searched in R0 ∈ [0.01, 100], IP ∈ [0.01, 1000]; a side that never
#     crosses the cutoff is recorded at the search bound and flagged open.
# "True" values are the parameters each tree was simulated with (Voznica's "true (or target) parameter values").
# Metrics per method and parameter: MRE = mean |est − true|/true (MAPE = 100·MRE), MRB = mean (est − true)/true, RMSE,
# coverage of the 95% interval (BD MLE only), and Voznica's paired two-sided z-test of each method's relative errors
# against the clipped MLE's. Per tree also the exact log-likelihood at the estimate minus that at the true values,
# ℓ(est) − ℓ(true), as Voznica's likelihood comparison (Methods Extended p. 35, Supplementary Tab. 3).
# Writes <run_dir>/comparison.csv (one row per tree and method) and <run_dir>/metrics.csv.

using PhyloPOMP, HDF5, CSV, DataFrames, Statistics, Optim, Printf, TOML, Lux, Random
using Distributions: Normal, cdf
include(joinpath(@__DIR__, "fit_lbdp.jl"))       # θ_of, exact_mle
include(joinpath(@__DIR__, "mle_models.jl"))     # MLE_MODELS, founder_probs, exact_mle_mtbd
include(joinpath(@__DIR__, "dl_models.jl"))      # cnn_cblv, ffnn_ss, load_keras

const BOX = (R0 = (1.0, 5.0), IP = (1.0, 10.0))
const SEARCH = (R0 = (0.01, 100.0), IP = (0.01, 1000.0))
const CUT = 3.841458820694124                    # χ²₁ 0.95 quantile

loglik(g, p, R0, IP) = lbdp_exact(g; θ_of(R0, IP, p)...)

## profile log-likelihood of `which` (:R0 or :IP) at value v, the other parameter maximized (Brent on its log)
profile(g, p, which, v) = begin
    other = which == :R0 ? :IP : :R0
    lo, hi = log.(SEARCH[other])
    f(x) = -(which == :R0 ? loglik(g, p, v, exp(x)) : loglik(g, p, exp(x), v))
    -Optim.minimum(optimize(f, lo, hi, Brent()))
end

## one side of the 95% interval by bisection on the log scale between the MLE and the search bound
endpoint(g, p, which, est, ℓ̂, bound) = begin
    gap(v) = 2 * (ℓ̂ - profile(g, p, which, v)) - CUT
    gap(bound) < 0 && return bound, true          # never crosses: open side
    a, b = log(est), log(bound)
    for _ in 1:40
        m = (a + b) / 2
        gap(exp(m)) < 0 ? (a = m) : (b = m)
    end
    exp((a + b) / 2), false
end

mle_row(g, p) = begin
    ex = exact_mle(g, p)
    ci = Dict{Symbol,Any}()
    for (which, est) in ((:R0, ex.R0), (:IP, ex.IP))
        lo, lo_open = endpoint(g, p, which, est, ex.ll, SEARCH[which][1])
        hi, hi_open = endpoint(g, p, which, est, ex.ll, SEARCH[which][2])
        ci[which] = (lo, hi, lo_open | hi_open)
    end
    (; ex.R0, ex.IP, ex.ll, R0_lo = ci[:R0][1], R0_hi = ci[:R0][2], R0_open = ci[:R0][3],
       IP_lo = ci[:IP][1], IP_hi = ci[:IP][2], IP_open = ci[:IP][3])
end

## paired two-sided z-test of relative errors (Voznica): z = mean(d)/(sd(d)/√n), d = |RE_a| − |RE_b|
ztest(a, b) = begin
    d = a .- b
    z = mean(d) / (std(d) / sqrt(length(d)))
    z, 2 * (1 - cdf(Normal(), abs(z)))
end

## short names for the targets (BD keeps the names of the first version: R0, IP)
const SHORT = Dict("R_0" => "R0", "Infectious_Period" => "IP", "Incubation_Period" => "INC",
                   "X_Transmission" => "X", "SS_Fraction" => "f")
const TIMEDEP = ("Infectious_Period", "Incubation_Period")

main(args) = begin
    dir = args[1]
    cfg = TOML.parsefile(joinpath(dir, "config.toml"))
    model = get(cfg, "model", "LBDP")
    lines = collect(eachline(joinpath(dir, "trees.tsv")))
    header = split(lines[1], '\t')
    col(name) = findfirst(==(name), header)::Int
    ip, inw, itree = col("p_sample"), col("newick"), col("tree")
    rows = [split(l, '\t') for l in lines[2:end]]
    f = joinpath(dir, "cblv.h5")
    X, S, resc, tree = h5read(f, "cblv"), h5read(f, "ss"), h5read(f, "rescale"), h5read(f, "tree")
    names = h5open(ff -> haskey(ff, "target_names") ? read(ff["target_names"]) : ["R_0", "Infectious_Period"], f)
    pars = [SHORT[n] for n in names]
    N = size(X, 2)
    ntest = min(10_000, N ÷ 20)
    ite = (N - ntest + 1):N
    byid = Dict(parse(Int, r[itree]) => r for r in rows)
    test = [byid[tree[k]] for k in ite]
    gs = [parse_newick(String(r[inw])) for r in test]
    ps = [parse(Float64, r[ip]) for r in test]
    truth = DataFrame(tree = tree[ite])
    for (n, s) in zip(names, pars)
        truth[!, "$(s)_true"] = [parse(Float64, r[col(n)]) for r in test]
    end
    truth.p_sample = ps
    truth.ntips = nsample.(gs)
    println("run $dir: model $model, stop=$(get(cfg, "stop", "time")), $N trees, $ntest test trees"); flush(stdout)
    est = Dict{String,DataFrame}()
    exact = haskey(MLE_MODELS, model)            # LBDP, BDEI, BDSS: all have an exact likelihood
    mm = get(MLE_MODELS, model, nothing)

    ## exact log-likelihood of test tree k at target values v (in the order of `pars`), each kept in its valid range
    lik(k, v) = model == "LBDP" ? loglik(gs[k], ps[k], max(v[1], 1e-6), max(v[2], 1e-6)) : begin
        w = [n == "SS_Fraction" ? clamp(x, 1e-6, 1 - 1e-6) : max(x, 1e-6) for (n, x) in zip(names, v)]
        mtbd_exact(gs[k], mm.table; θ = mm.rates(w, ps[k]), founder = founder_probs(mm, w))
    end

    ## BD: MLE from the exact likelihood, with profile intervals, threaded over trees.
    if model == "LBDP"
        t = @elapsed begin
            mles = Vector{Any}(undef, ntest)
            Threads.@threads :dynamic for k in 1:ntest
                mles[k] = mle_row(gs[k], ps[k])
            end
        end
        @printf("exact MLE + profile CIs: %.1f s for %d trees (%.2f s per tree on %d threads)\n",
                t, ntest, t / ntest, Threads.nthreads()); flush(stdout)
        ℓ̂ = [m.ll for m in mles]
        est["MLE"] = DataFrame(R0 = [m.R0 for m in mles], IP = [m.IP for m in mles],
                               R0_lo = [m.R0_lo for m in mles], R0_hi = [m.R0_hi for m in mles],
                               IP_lo = [m.IP_lo for m in mles], IP_hi = [m.IP_hi for m in mles],
                               R0_open = [m.R0_open for m in mles], IP_open = [m.IP_open for m in mles])
        est["MLE (clipped)"] = DataFrame(R0 = clamp.(est["MLE"].R0, BOX.R0...), IP = clamp.(est["MLE"].IP, BOX.IP...))
    elseif exact
        ## BDEI, BDSS: MLE from mtbd_exact, threaded over trees; no profile intervals (three or four parameters)
        t = @elapsed begin
            fits = Vector{Any}(undef, ntest)
            done = Threads.Atomic{Int}(0); t0 = time()
            Threads.@threads :dynamic for k in 1:ntest
                fits[k] = exact_mle_mtbd(gs[k], mm, ps[k])
                (Threads.atomic_add!(done, 1) + 1) % 50 == 0 &&
                    (@printf("  exact MLE: %d of %d trees, %.0f s\n", done[], ntest, time() - t0); flush(stdout))
            end
        end
        @printf("exact MLE: %.1f s for %d trees (%.2f s per tree on %d threads)\n",
                t, ntest, t / ntest, Threads.nthreads()); flush(stdout)
        ℓ̂ = [f[2] for f in fits]
        est["MLE"] = DataFrame([s => [f[1][j] for f in fits] for (j, s) in enumerate(pars)])
        est["MLE (clipped)"] = DataFrame([s => clamp.(est["MLE"][!, s], mm.box[n]...) for (n, s) in zip(names, pars)])
    end

    ## our networks
    for net in ("cnn", "ffnn"), arch in ("released", "paper", "python")
        file = joinpath(dir, "$(net)_$(arch)_test.csv")
        isfile(file) || continue
        d = CSV.read(file, DataFrame)
        d.tree == truth.tree || error("$file: test trees differ from cblv.h5's split")
        est["$(uppercase(net)) $arch"] = DataFrame([s => d[!, "$(n)_pred"] for (n, s) in zip(names, pars)])
    end

    ## phylodeep's pretrained models, as a reference (BD's are named BD_*; the FFNN needs its exported scaler)
    models = get(ENV, "PHYLODEEP_MODELS", joinpath(homedir(),
        "Desktop/Projects/efficiency/DeepLearningPhylodeep/venv/lib/python3.12/site-packages/phylodeep/pretrained_models/models"))
    tag = (model == "LBDP" ? "BD" : model) * "_" * uppercase(cfg["size"])
    unresc(y) = DataFrame([s => Float64.(y[j, :]) .* (names[j] in TIMEDEP ? resc[ite] : 1.0) for (j, s) in enumerate(pars)])
    cnnfile = joinpath(models, "$(tag)_CNN.h5")
    if isfile(cnnfile)
        m = cnn_cblv(; nout = length(names)); θ, st = Lux.setup(Xoshiro(0), m)
        y, _ = m(X[:, ite], load_keras(m, θ, cnnfile), Lux.testmode(st))
        est["phylodeep CNN (pretrained)"] = unresc(y)
    end
    ffnnfile = joinpath(models, "$(tag)_FFNN.h5")
    scaler = joinpath(@__DIR__, "phylodeep_scalers", "$(tag)_FFNN.csv")
    if isfile(ffnnfile) && isfile(scaler)
        μ, σ = [parse.(Float64, split(l, ',')) for l in eachline(scaler)]
        m = ffnn_ss(; nout = length(names)); θ, st = Lux.setup(Xoshiro(0), m)
        y, _ = m(Float32.((S[:, ite] .- μ) ./ σ), load_keras(m, θ, ffnnfile), Lux.testmode(st))
        est["phylodeep FFNN (pretrained)"] = unresc(y)
    end

    ## per-tree table, with the exact log-likelihood at each estimate and at the true values
    if exact
        ℓtrue = zeros(ntest)
        Threads.@threads :dynamic for k in 1:ntest
            ℓtrue[k] = lik(k, [truth[k, "$(s)_true"] for s in pars])
        end
    end
    long = DataFrame()
    for (name, e) in est
        d = hcat(truth, DataFrame(method = fill(name, ntest)))
        for s in pars
            d[!, "$(s)_est"] = e[!, s]
        end
        if exact
            ℓe = zeros(ntest)
            Threads.@threads :dynamic for k in 1:ntest
                ℓe[k] = lik(k, [e[k, s] for s in pars])
            end
            d.ll_est = ℓe; d.ll_true = ℓtrue; d.ll_max = ℓ̂; d.ll_vs_true = ℓe .- ℓtrue
        end
        if name == "MLE" && model == "LBDP"
            d = hcat(d, e[:, [:R0_lo, :R0_hi, :IP_lo, :IP_hi, :R0_open, :IP_open]])
        end
        long = vcat(long, d; cols = :union)
    end
    CSV.write(joinpath(dir, "comparison.csv"), long)

    ## metrics; the paired z-test is against the clipped MLE when there is one, otherwise phylodeep's pretrained FFNN
    refname = exact ? "MLE (clipped)" : haskey(est, "phylodeep FFNN (pretrained)") ? "phylodeep FFNN (pretrained)" : ""
    re(e, t) = (e .- t) ./ t
    met = DataFrame()
    order = sort(collect(keys(est)); by = n -> (startswith(n, "MLE") ? 0 : startswith(n, "phylodeep") ? 2 : 1, n))
    for name in order, par in pars
        e = est[name][!, par]
        tv = truth[!, "$(par)_true"]
        r = re(e, tv)
        cov = name == "MLE" && model == "LBDP" ?
              mean((est["MLE"][!, Symbol(par, :_lo)] .≤ tv) .& (tv .≤ est["MLE"][!, Symbol(par, :_hi)])) : NaN
        z, pz = (refname == "" || name == refname) ? (NaN, NaN) : ztest(abs.(r), abs.(re(est[refname][!, par], tv)))
        ll = exact ? long[long.method .== name, :ll_vs_true] : [NaN]
        push!(met, (method = name, parameter = par, n = length(e), MRE = mean(abs.(r)), MRB = mean(r),
                    RMSE = sqrt(mean((e .- tv) .^ 2)), coverage95 = cov, mean_ll_vs_true = mean(ll),
                    share_ll_above_true = mean(ll .> 0), reference = refname,
                    z_vs_reference = z, p_vs_reference = pz))
    end
    CSV.write(joinpath(dir, "metrics.csv"), met)
    show(stdout, met[:, [:method, :parameter, :MRE, :MRB, :RMSE, :coverage95, :mean_ll_vs_true, :p_vs_reference]];
         allrows = true); println()
end

main(ARGS)
