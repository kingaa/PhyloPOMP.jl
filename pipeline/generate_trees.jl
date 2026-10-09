# M5: training genealogies for the deep-learning step: BD (the `@mgp` LBDP table), BDEI or BDSS.
#
#   julia -t 16 --project=pipeline pipeline/generate_trees.jl <small|large> <n_trees> <out_dir> [seed]
#                                                               [stop=time|tips] [model=BD|BDEI|BDSS]
#
# Per tree, drawn once and independently, uniformly (user decision), γ = 1/IP, removal splits into sampling χ = γp
# and death μ = γ(1 − p):
#   BD   (LBDP): R0 ~ U(1,5), IP ~ U(1,10), p ~ U(0.01,1); λ = R0·γ.
#   BDEI: R0 ~ U(1,5), IP ~ U(1,10), incubation period = IP × U(0.2,5), p; σ = 1/incubation period, λ = R0·γ.
#         Voznica's text calls the ratio ε/γ, but phylodeep's pretrained BDEI networks are unbiased only on the ratio
#         incubation period / IP = γ/ε drawn uniformly (bias of their IP and incubation estimates: −13.6% and +22.2%
#         under ε/γ ~ U(0.2,5), −1.5% and −1.5% under IP × U(0.2,5); 1,500 trees each, 2026-10-08).
#   BDSS: R0 ~ U(1,5), IP ~ U(1,10), X_SS ~ U(3,10), f_SS ~ U(0.05,0.2), p (Voznica, Methods Extended p. 30):
#         normal hosts transmit at total rate B = R0·γ/(f·X + 1 − f), superspreaders at X·B; a new host is a
#         superspreader with probability f: λ_nn = (1−f)B, λ_ns = fB, λ_sn = X(1−f)B, λ_ss = XfB.
#   Ranges and BDSS relations as Voznica's Supplementary Table 4 (efficiency/DeepLearningPhylodeep/
#   supplementary_info-1-17.pdf, PDF p. 5). For BDEI that table, like the Methods text, gives the incubation factor as
#   ε/γ ~ U(0.2,5); the generator uses its inverse, IP × U(0.2,5), because phylodeep's pretrained BDEI networks are
#   unbiased only on that (above).
# One infectious founder (BDEI: I; BDSS: a superspreader with probability f_SS, as Voznica). stop=time is BD only.
#
# stop=time (default; port of generate_data_phylopomp.R, efficiency/DeepLearningPhylodeep/PhyloDeepArrayRun):
#   target tips ~ U(min, max), real; run to time T = log(target)/r with r = λ − μ − χ, or 999 if r ≤ 0.01; a draw
#   whose T exceeds 500 is skipped. Up to 100 simulations per draw, each stopped after 1 s; the first with
#   min ≤ tips ≤ max is kept, otherwise a new draw. A run is also stopped once it has more than max samples, which
#   R would reject anyway. The tree ends at its last sample; column `time` is T.
# stop=tips (Voznica et al. 2022, Methods Extended p. 31): tree size t ~ U{min, …, max}; each simulation stops at the
#   t-th sample, so the tree ends at its last sample; if the epidemic dies out first, restart, up to 100 times, then
#   discard the parameter set and draw a new one. No timeout (Voznica has none). Column `time` is the stop time.
#   Voznica drew parameters by Latin hypercube; here they stay iid uniform.
#
# Writes <out_dir>/trees.tsv (tree, R_0, Infectious_Period, p_sample, ntips, attempts, newick, time) and
# <out_dir>/config.toml. Each tree has its own random stream, seeded from (seed, tree), so the output does
# not depend on the number of threads, and a stopped run resumes where it left off (with the same stop rule).

using PhyloPOMP
using Random: Xoshiro
using TOML, Dates

struct Stop <: Exception end

const SIZES = Dict("small" => (50, 199), "large" => (200, 500))
const MODELS = Dict("BD" => PhyloPOMP.LBDP, "BDEI" => PhyloPOMP.BDEI, "BDSS" => PhyloPOMP.BDSS)
const TARGETS = Dict("BD" => ["R_0", "Infectious_Period"],
                     "BDEI" => ["R_0", "Infectious_Period", "Incubation_Period"],
                     "BDSS" => ["R_0", "Infectious_Period", "X_Transmission", "SS_Fraction"])
const STOPS = ("time", "tips")
const MAX_ATTEMPTS = 100
const TIMEOUT = 1.0          # seconds per simulation, stop=time only
const MAX_TIME = 500.0
const BATCH = 1000           # trees per append to trees.tsv
issample(ev) = ev.type == PhyloPOMP.SAMPLE

simulate_model(M, θ, x0, graft, tmax, hook, rng) =
    PhyloPOMP._simulate(M, θ, x0, graft, 0.0, tmax, rng, PhyloPOMP.Unstructured, nothing, hook)
simulate_lbdp(θ, tmax, hook, rng) = simulate_model(PhyloPOMP.LBDP, θ, (n = 1,), [1], tmax, hook, rng)

## stop=time: one simulation to time T; `nothing` when stopped by the timeout or the tip cap.
sim_time(θ, T, maxtips, rng) = begin
    t_start = time_ns()
    nsamp = 0
    hook(G, inv, ev, t, x) = begin
        isnothing(ev) && return
        issample(ev) && (nsamp += 1)
        (nsamp > maxtips || time_ns() - t_start > TIMEOUT * 1e9) && throw(Stop())
    end
    try
        simulate_lbdp(θ, T, hook, rng)
    catch e
        e isa Stop ? nothing : rethrow()
    end
end

## stop=tips: one simulation stopped at its `t`-th sample, finished as `simulate` finishes a run (end time,
## prune!, repair!, check_simulated!); `nothing` when the epidemic dies out first.
sim_tips(θ, t, rng; M = PhyloPOMP.LBDP, x0 = (n = 1,), graft = [1]) = begin
    nsamp = 0
    stopped = Ref{Any}(nothing)
    hook(G, inv, ev, tev, x) = begin
        isnothing(ev) && return
        if issample(ev) && (nsamp += 1) == t
            stopped[] = (G, tev)
            throw(Stop())
        end
    end
    try
        simulate_model(M, θ, x0, graft, Inf, hook, rng)
        return nothing                       # died out before t samples
    catch e
        e isa Stop || rethrow()
    end
    G, tstop = stopped[]
    G.time = tstop
    PhyloPOMP.prune!(G)
    PhyloPOMP.repair!(G)
    PhyloPOMP.check_simulated!(G)
    G
end

## The commit of the code that made the trees; "-dirty" when src/ or pipeline/ has uncommitted or untracked changes.
commit_stamp() = begin
    dir = @__DIR__
    head = strip(read(`git -C $dir rev-parse --short HEAD`, String))
    changed = strip(read(`git -C $dir status --porcelain -- $(joinpath(dir, "..", "src")) $dir`, String))
    isempty(changed) ? head : head * "-dirty"
end

row(i, targets, p, g, attempt, tend) =
    join((i, targets..., p, nsample(g), attempt, only(newick(g; sigdigits = 12)), tend), '\t')   # one founder: one root

## BDEI and BDSS draws (Voznica's parameters), their table parameters, and the founder.
draw(::Val{:BDEI}, rng) = begin
    R0 = 1.0 + 4.0 * rand(rng); ip = 1.0 + 9.0 * rand(rng); inc = ip * (0.2 + 4.8 * rand(rng))
    p = 0.01 + 0.99 * rand(rng)
    γ = 1 / ip
    θ = (σ = 1 / inc, λ = R0 * γ, μ = γ * (1 - p), χ = γ * p)
    (R0, ip, inc), p, θ, (E = 0, I = 1), [0, 1]
end
draw(::Val{:BDSS}, rng) = begin
    R0 = 1.0 + 4.0 * rand(rng); ip = 1.0 + 9.0 * rand(rng); X = 3.0 + 7.0 * rand(rng); f = 0.05 + 0.15 * rand(rng)
    p = 0.01 + 0.99 * rand(rng)
    γ = 1 / ip; B = R0 * γ / (f * X + 1 - f)
    θ = (λ_nn = (1 - f) * B, λ_ns = f * B, λ_sn = X * (1 - f) * B, λ_ss = X * f * B, μ = γ * (1 - p), χ = γ * p)
    super = rand(rng) < f
    (R0, ip, X, f), p, θ, super ? (N = 0, S = 1) : (N = 1, S = 0), super ? [0, 1] : [1, 0]
end

## Tree `i`: draw parameters until one gives a tree. Returns the TSV row and the number of discarded draws.
make_tree(i, seed, mintips, maxtips, stop, model = "BD") = begin
    rng = Xoshiro(hash((seed, i)))
    skipped = 0
    while model != "BD"
        targets, p, θ, x0, graft = draw(Val(Symbol(model)), rng)
        t = rand(rng, mintips:maxtips)
        for attempt in 1:MAX_ATTEMPTS
            g = sim_tips(θ, t, rng; M = MODELS[model], x0 = x0, graft = graft)
            isnothing(g) || return row(i, targets, p, g, attempt, g.time), skipped
        end
        skipped += 1
    end
    while true
        R0 = 1.0 + 4.0 * rand(rng)
        inf_period = 1.0 + 9.0 * rand(rng)
        p_sample = 0.01 + 0.99 * rand(rng)
        γ = 1 / inf_period
        θ = (λ = R0 * γ, μ = γ * (1 - p_sample), ψ = 0.0, χ = γ * p_sample)
        if stop == "tips"
            t = rand(rng, mintips:maxtips)
            for attempt in 1:MAX_ATTEMPTS
                g = sim_tips(θ, t, rng)
                isnothing(g) || return row(i, (R0, inf_period), p_sample, g, attempt, g.time), skipped
            end
        else
            r = θ.λ - θ.μ - θ.χ
            target = mintips + (maxtips - mintips) * rand(rng)
            T = r > 0.01 ? log(target) / r : 999.0
            if T ≤ MAX_TIME
                for attempt in 1:MAX_ATTEMPTS
                    g = sim_time(θ, T, maxtips, rng)
                    isnothing(g) && continue
                    mintips ≤ nsample(g) ≤ maxtips && return row(i, (R0, inf_period), p_sample, g, attempt, T), skipped
                end
            end
        end
        skipped += 1
    end
end

main(args) = begin
    pos = filter(a -> !occursin('=', a), args)
    kw = Dict(split(a, '='; limit = 2) for a in args if occursin('=', a))
    length(pos) ≥ 3 || error("usage: generate_trees.jl <small|large> <n_trees> <out_dir> [seed] [stop=time|tips]")
    size, n, outdir = pos[1], parse(Int, pos[2]), pos[3]
    seed = length(pos) ≥ 4 ? parse(Int, pos[4]) : 42
    stop = String(get(kw, "stop", "time"))
    model = uppercase(String(get(kw, "model", "BD")))
    haskey(SIZES, size) || error("size must be small or large")
    stop ∈ STOPS || error("stop must be time or tips")
    haskey(MODELS, model) || error("model must be BD, BDEI or BDSS")
    model == "BD" || stop == "tips" || error("stop=time is defined for BD only; use stop=tips")
    mintips, maxtips = SIZES[size]
    mkpath(outdir)
    tsv = joinpath(outdir, "trees.tsv")
    cfg = joinpath(outdir, "config.toml")
    header = join(["tree"; TARGETS[model]; "p_sample"; "ntips"; "attempts"; "newick"; "time"], '\t')
    done = isfile(tsv) ? max(countlines(tsv) - 1, 0) : 0
    if done > 0 && isfile(cfg)
        oldcfg = TOML.parsefile(cfg)
        old = get(oldcfg, "stop", "time")
        old == stop || error("$outdir was generated with stop=$old; resume it with the same rule or use a new directory")
        get(oldcfg, "model", "LBDP") in (model, model == "BD" ? "LBDP" : model) ||
            error("$outdir was generated with model=$(oldcfg["model"])")
    end
    done == 0 && write(tsv, header * "\n")
    open(cfg, "w") do io
        TOML.print(io, Dict(
            "task" => "generate", "model" => model == "BD" ? "LBDP" : model, "targets" => TARGETS[model],
            "size" => size, "n_trees" => n, "seed" => seed,
            "stop" => stop, "min_tips" => mintips, "max_tips" => maxtips,
            "priors" => Dict("R_0" => [1.0, 5.0], "Infectious_Period" => [1.0, 10.0], "p_sample" => [0.01, 1.0],
                             "incubation_period_over_IP" => model == "BDEI" ? [0.2, 5.0] : Float64[],
                             "X_Transmission" => model == "BDSS" ? [3.0, 10.0] : Float64[],
                             "SS_Fraction" => model == "BDSS" ? [0.05, 0.2] : Float64[]),
            "max_attempts" => MAX_ATTEMPTS,
            "timeout_seconds" => stop == "time" ? TIMEOUT : 0.0,     # 0: none
            "max_time" => stop == "time" ? MAX_TIME : 0.0,           # 0: none
            "threads" => Threads.nthreads(), "started" => string(now()),
            ## per-tree streams are Xoshiro(hash((seed, tree))): identical only under the same Julia version
            "julia_version" => string(VERSION),
            "phylopomp_commit" => commit_stamp(),
        ))
    end
    println("$model $(size), stop=$stop: trees $(done + 1) to $n into $tsv, $(Threads.nthreads()) threads")
    skipped = Threads.Atomic{Int}(0)
    t0 = time()
    for lo in (done + 1):BATCH:n
        idx = lo:min(lo + BATCH - 1, n)
        rows = Vector{String}(undef, length(idx))
        Threads.@threads :dynamic for k in eachindex(idx)
            rows[k], s = make_tree(idx[k], seed, mintips, maxtips, stop, model)
            Threads.atomic_add!(skipped, s)
        end
        open(tsv, "a") do io
            foreach(r -> println(io, r), rows)
        end
        println("  $(last(idx))/$n written, $(skipped[]) draws discarded, $(round(time() - t0, digits = 1)) s")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
