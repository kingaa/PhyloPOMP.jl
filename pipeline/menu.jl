# M8: build a run config with menus, then run it (pipeline/run.jl does the work).
#
#   julia --project=pipeline pipeline/menu.jl
#
# Choices, in order: model family -> model -> tree size -> tasks -> parameters (Enter keeps the default shown) ->
# where to save the TOML -> run now or later. The TOML can be edited and run later with
#   julia --project=pipeline pipeline/run.jl <file>.toml
# The linear family (birth-death models without a susceptible pool) has the deep-learning pipeline. The nonlinear
# models (SEIR, MERS, SI2R, ...) have filters and fits (fit_seir_profile.jl for SEIR) but no deep-learning tasks yet.

using REPL.TerminalMenus, TOML
include(joinpath(@__DIR__, "run.jl"))

## the choices as a config Dict, checked by run.jl's `check_config`; separate from the prompts so it can be tested
function build_config(; dir, model, size, tasks, n_trees = 10000, seed = 42, stop = "tips", net = "both",
                      arch = "both", epochs = 200, batch = 500, patience = 15, device = "cpu", ntrees = 4,
                      threads = 16, force = false)
    check_config(Dict(
        "run" => Dict("dir" => dir, "model" => model, "size" => size, "tasks" => collect(tasks),
                      "threads" => threads, "force" => force),
        "generate" => Dict("n_trees" => n_trees, "seed" => seed, "stop" => stop),
        "train" => Dict("net" => net, "arch" => arch, "epochs" => epochs, "batch" => batch,
                        "patience" => patience, "device" => device),
        "fit" => Dict("ntrees" => ntrees),
    ))
end

pick(title, options) = (println(title); options[request(RadioMenu(options; pagesize = 8))])
function ask(prompt, default)
    print("$prompt [$default]: ")
    s = strip(readline())
    isempty(s) ? default : default isa Integer ? parse(Int, s) : String(s)
end

function interactive()
    family = pick("Model family:", ["linear (birth-death: BD, BDEI, BDSS)", "nonlinear (SEIR, MERS, SI2R, ...)"])
    if startswith(family, "nonlinear")
        println("Nonlinear models have filters and fits but no deep-learning tasks yet. For SEIR run:\n" *
                "  julia -t 16 --project=pipeline pipeline/fit_seir_profile.jl")
        return nothing
    end
    model = pick("Model:", ["BD", "BDEI", "BDSS"])
    size = pick("Tree size:", model == "BDSS" ? ["large"] : ["small", "large"])
    alltasks = model == "BD" ? ["generate", "encode", "train", "compare", "fit"] : ["generate", "encode", "train", "compare"]
    println("Tasks (space to select, d when done):")
    sel = request(MultiSelectMenu(alltasks; selected = 1:4))
    tasks = alltasks[sort(collect(sel))]
    dir = ask("Run directory (relative to PhyloPOMP.jl)", "output/$(lowercase(model))_$(size)")
    kw = Dict{Symbol,Any}(:dir => dir, :model => model, :size => size, :tasks => tasks)
    if "generate" in tasks
        kw[:n_trees] = ask("Number of trees", 10000)
        kw[:seed] = ask("Seed", 42)
        kw[:stop] = model == "BD" ? pick("Stop rule:", ["tips", "time"]) : "tips"
    end
    if "train" in tasks
        kw[:device] = pick("Train on:", ["gpu", "cpu"])
        kw[:arch] = pick("Architecture:", ["both", "released", "paper"])
        kw[:epochs] = ask("Maximum epochs", 200)
        kw[:batch] = ask("Batch size", 500)
    end
    "fit" in tasks && (kw[:ntrees] = ask("Trees to fit by mif", 4))
    cfg = build_config(; kw...)
    file = ask("Save the config as", joinpath(dir, "run.toml"))
    path = joinpath(ROOT, file)
    mkpath(dirname(path))
    open(io -> TOML.print(io, cfg), path, "w")
    println("written $file")
    pick("Run it now?", ["yes", "no"]) == "yes" ? run_config(cfg) :
        println("later: julia --project=pipeline pipeline/run.jl $file")
    cfg
end

if abspath(PROGRAM_FILE) == @__FILE__
    interactive()
end
