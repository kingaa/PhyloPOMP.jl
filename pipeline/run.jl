# M8: run the pipeline from one TOML file.
#
#   julia --project=pipeline pipeline/run.jl <config.toml>
#
# The config names a run directory and the tasks to run; each task calls one pipeline script in the pipeline
# environment and writes its printed output to <dir>/<task>.log. A task is skipped, unless `force = true`, when its
# output exists and is not older than its inputs; `generate` also needs all `n_trees` rows in trees.tsv, so an
# interrupted generation resumes. The config is copied into the run directory as run_config.toml.
#
#   [run]
#   dir     = "output/bd_small"     # relative to the PhyloPOMP.jl directory
#   model   = "BD"                  # BD, BDEI or BDSS (the deep-learning tasks)
#   size    = "small"               # small (50-199 tips) or large (200-500)
#   tasks   = ["generate", "encode", "train", "compare"]   # also "fit" (BD: mif against the exact MLE)
#   threads = 16
#   force   = false
#
#   [generate]  n_trees = 10000, seed = 42, stop = "tips"
#   [train]     net = "both", arch = "both", epochs = 200, batch = 500, patience = 15, device = "gpu"
#   [fit]       ntrees = 4
#
# Tasks and what they run:
#   generate  generate_trees.jl <size> <n_trees> <dir> <seed> stop=<stop> model=<model>   -> trees.tsv
#   encode    encode_trees.jl <dir>                                                       -> cblv.h5
#   train     train_dl.jl <dir> net= arch= epochs= batch= patience= device=              -> <net>_<arch>_test.csv
#   compare   compare.jl <dir>                                                            -> metrics.csv
#   fit       fit_lbdp.jl <dir>/trees.tsv <ntrees>   (BD only)                           -> fit_lbdp.csv

using TOML, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const TASKS = ("generate", "encode", "train", "compare", "fit")
const DEFAULTS = Dict(
    "run" => Dict("model" => "BD", "size" => "small", "tasks" => ["generate", "encode", "train", "compare"],
                  "threads" => 16, "force" => false),
    "generate" => Dict("n_trees" => 10000, "seed" => 42, "stop" => "tips"),
    "train" => Dict("net" => "both", "arch" => "both", "epochs" => 200, "batch" => 500, "patience" => 15,
                    "device" => "cpu"),
    "fit" => Dict("ntrees" => 4),
)

"""
    load_config(file) -> Dict

The config with defaults filled in. Throws `ArgumentError` for an unknown model, size, task, stop rule, network,
architecture or device, or when `fit` is asked for a model other than BD.
"""
load_config(file::AbstractString) = check_config(TOML.parsefile(file))

function check_config(c::Dict)
    cfg = Dict(k => merge(v, get(c, k, Dict())) for (k, v) in DEFAULTS)
    haskey(get(c, "run", Dict()), "dir") || throw(ArgumentError("[run] needs dir"))
    r = cfg["run"]
    r["model"] = uppercase(r["model"])
    r["model"] in ("BD", "BDEI", "BDSS") || throw(ArgumentError("model must be BD, BDEI or BDSS"))
    r["size"] in ("small", "large") || throw(ArgumentError("size must be small or large"))
    for t in r["tasks"]
        t in TASKS || throw(ArgumentError("unknown task $t; tasks are $(join(TASKS, ", "))"))
    end
    "fit" in r["tasks"] && r["model"] != "BD" &&
        throw(ArgumentError("fit is for BD only (fit_lbdp.jl): the particle filter is not usable on BDEI and BDSS trees of these sizes yet"))
    r["model"] == "BDSS" && r["size"] == "small" && @warn "phylodeep has BDSS networks for large trees only"
    cfg["generate"]["stop"] in ("tips", "time") || throw(ArgumentError("stop must be tips or time"))
    r["model"] != "BD" && cfg["generate"]["stop"] == "time" && throw(ArgumentError("stop=time is for BD only"))
    tr = cfg["train"]
    tr["net"] in ("cnn", "ffnn", "both") || throw(ArgumentError("net must be cnn, ffnn or both"))
    tr["arch"] in ("released", "paper", "python", "both") || throw(ArgumentError("arch must be released, paper, python or both"))
    tr["device"] in ("cpu", "gpu") || throw(ArgumentError("device must be cpu or gpu"))
    cfg
end

## the command, the file it produces, of each task
function task_command(task, cfg)
    r = cfg["run"]; dir = joinpath(ROOT, r["dir"]); g = cfg["generate"]; tr = cfg["train"]
    script(s) = joinpath(ROOT, "pipeline", s)
    jl = `$(joinpath(Sys.BINDIR, "julia")) -t $(r["threads"]) --project=$(joinpath(ROOT, "pipeline"))`
    if task == "generate"
        `$jl $(script("generate_trees.jl")) $(r["size"]) $(g["n_trees"]) $dir $(g["seed"]) stop=$(g["stop"]) model=$(r["model"])`,
            "trees.tsv"
    elseif task == "encode"
        `$jl $(script("encode_trees.jl")) $dir`, "cblv.h5"
    elseif task == "train"
        last = tr["net"] == "cnn" ? "cnn" : "ffnn"
        lastarch = tr["arch"] == "both" ? "paper" : tr["arch"]
        `$jl $(script("train_dl.jl")) $dir net=$(tr["net"]) arch=$(tr["arch"]) epochs=$(tr["epochs"]) batch=$(tr["batch"]) patience=$(tr["patience"]) device=$(tr["device"])`,
            "$(last)_$(lastarch)_test.csv"
    elseif task == "compare"
        `$jl $(script("compare.jl")) $dir`, "metrics.csv"
    else
        `$jl $(script("fit_lbdp.jl")) $(joinpath(dir, "trees.tsv")) $(cfg["fit"]["ntrees"])`, "fit_lbdp.csv"
    end
end

## the files a task reads, in the run directory
task_inputs(task, cfg) =
    task == "encode" ? ["trees.tsv"] :
    task == "train" ? ["cblv.h5"] :
    task == "compare" ? ["cblv.h5", task_command("train", cfg)[2]] :
    task == "fit" ? ["trees.tsv"] : String[]

"""
    task_done(task, output, dir, cfg) -> Bool

Whether `task` can be skipped: its `output` exists and is not older than any input of the task, and for `generate`
`trees.tsv` holds all `n_trees` rows (an interrupted run leaves fewer; `generate_trees.jl` then resumes the file).
"""
function task_done(task, output, dir, cfg)
    f = joinpath(dir, output)
    isfile(f) || return false
    task == "generate" && return countlines(f) - 1 >= cfg["generate"]["n_trees"]
    all(i -> !isfile(joinpath(dir, i)) || mtime(f) >= mtime(joinpath(dir, i)), task_inputs(task, cfg))
end

function run_config(cfg)
    r = cfg["run"]; dir = joinpath(ROOT, r["dir"])
    mkpath(dir)
    open(io -> TOML.print(io, cfg), joinpath(dir, "run_config.toml"), "w")
    println("run $(r["dir"]): $(r["model"]) $(r["size"]); tasks $(join(r["tasks"], ", "))")
    for task in r["tasks"]
        cmd, output = task_command(task, cfg)
        if !r["force"] && task_done(task, output, dir, cfg)
            println("  $task: skipped, $output is complete")
            continue
        end
        log = joinpath(dir, "$task.log")
        t = @elapsed ok = success(pipeline(cmd; stdout = log, stderr = log))
        ok || error("$task failed after $(round(t, digits = 1)) s; see $log")
        @printf("  %s: done in %.1f s (log %s)\n", task, t, relpath(log, ROOT))
    end
    if isfile(joinpath(dir, "metrics.csv"))
        println("metrics: ", relpath(joinpath(dir, "metrics.csv"), ROOT))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run.jl <config.toml>")
    run_config(load_config(ARGS[1]))
end
