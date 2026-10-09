# M6: train CNN-CBLV and FFNN-SS on cblv.h5 from encode_trees.jl, as train_bd_models.py does in Keras.
#
#   julia -t 16 --project=pipeline pipeline/train_dl.jl <run_dir> [key=value ...]
#
# Keys (defaults): net=both arch=both epochs=200 batch=8000 patience=15 lr=0.001 seed=42 device=cpu.
# device=gpu trains on the CUDA GPU (LuxCUDA): data, parameters and states go to the GPU once; batches are drawn by
#   shuffling indices on the CPU. cuDNN is not bit-reproducible, so GPU runs differ slightly between runs; CPU runs
#   are deterministic for a seed.
# net: cnn (CNN-CBLV on `cblv`), ffnn (FFNN-SS on `ss`, standardized over all 99 inputs with the training rows'
#      mean and population SD, as sklearn's StandardScaler in phylodeep and train_bd_models.py), or both.
# arch (user decision 2026-10-08: train both):
#   released  phylodeep's pretrained models: no dropout, ELU on every hidden layer.
#   paper     Methods Extended text: dropout 0.5 in the dense part, ELU (its figure caption) on the 8-unit layer.
#   python    train_bd_models.py: dropout 0.5, linear 8-unit layer (only to compare with earlier Python runs).
#   both      released, then paper.
# phylodeep's released models (0.9 .h5 configs; 0.3 JSONs) have no Dropout and an ELU 8-unit layer; the Methods text
# says dropout 0.5 and contradicts its own caption on the 8-unit layer. How the 2022 models were trained is not
# established (the unattached-Dropout bug in scripts/phylodeep_dropout_check.py is from a different paper's code).
# Targets: the `target_names` of cblv.h5 (BD: R0, infectious period; BDEI: + incubation period; BDSS: + X_SS, f_SS),
# periods divided by the rescale factor, as phylodeep; loss MAPE (Keras "mape"); Adam; early stopping on the
# validation loss, keeping the best parameters. Split: last min(10000, N÷20) trees test, the min(10000, N÷10)
# before them validation, the rest training (train_bd_models.py's split leaves no training trees below ~210k).
# Writes <run_dir>/<net>_<arch>_{params.jld2,history.csv,test.csv}; the ffnn params file holds the scaler.

using Lux, LuxCUDA, Optimisers, Zygote, HDF5, JLD2, Random, Statistics, Printf
include(joinpath(@__DIR__, "dl_models.jl"))

const DEFAULTS = Dict("net" => "both", "arch" => "both", "epochs" => "200", "batch" => "8000", "patience" => "15",
                      "lr" => "0.001", "seed" => "42", "device" => "cpu")
const TIME_DEPENDENT = ("Infectious_Period", "Incubation_Period")   # phylodeep's TIME_DEPENDENT_COLUMNS
const ARCHS = Dict("released" => (dropout = 0.0, last = elu), "paper" => (dropout = 0.5, last = elu),
                   "python" => (dropout = 0.5, last = identity))

## Keras "mape": 100·mean(|y − ŷ| / max(|y|, ε)).
mape(ŷ, y) = 100 * mean(abs.(y .- ŷ) ./ max.(abs.(y), 1f-7))

train(net, arch, d, opt) = begin
    epochs, batch, patience, seed = (parse(Int, opt[k]) for k in ("epochs", "batch", "patience", "seed"))
    lr = parse(Float32, opt["lr"])
    (; Y, T, resc, itr, iva, ite, dir, tree, names) = d
    nout = length(names)
    X, scaler = net == "cnn" ? (d.cblv, nothing) : standardize(d.ss, itr)
    println("== $(uppercase(net)), arch $arch: $(ARCHS[arch])")
    rng = Xoshiro(seed)
    dev = opt["device"] == "gpu" ? gpu_device() : cpu_device()
    opt["device"] == "gpu" && !(dev isa CUDADevice) && error("device=gpu but no functional CUDA GPU")
    model = net == "cnn" ? cnn_cblv(; ARCHS[arch]..., nout) : ffnn_ss(; ARCHS[arch]..., ninput = size(X, 1), nout)
    ps, st = Lux.setup(rng, model)
    ps, st = dev(ps), dev(st)
    Xd, Yd = dev(X), dev(Y)
    state = Training.TrainState(model, ps, st, Adam(lr))
    loss(model, ps, st, (x, y)) = begin
        ŷ, st2 = model(x, ps, st)
        mape(ŷ, y), st2, (;)
    end
    evaluate(ps, idx) = mape(first(model(Xd[:, idx], ps, Lux.testmode(state.states))), Yd[:, idx])
    out(s) = joinpath(dir, "$(net)_$(arch)_$s")

    best, bestps, wait = Inf32, state.parameters, 0
    Xtr, Ytr = Xd[:, itr], Yd[:, itr]            # copied once, not every epoch
    open(out("history.csv"), "w") do io
        println(io, "epoch,train_loss,val_loss,seconds")
        for epoch in 1:epochs
            t = @elapsed begin
                tl = 0.0; nb = 0
                perm = randperm(rng, length(itr))
                for lo in 1:batch:length(perm)
                    idx = perm[lo:min(lo + batch - 1, end)]
                    x, y = Xtr[:, idx], Ytr[:, idx]
                    _, l, _, state = Training.single_train_step!(AutoZygote(), loss, (x, y), state)
                    tl += l; nb += 1
                end
                vl = evaluate(state.parameters, iva)
            end
            @printf(io, "%d,%.6g,%.6g,%.2f\n", epoch, tl / nb, vl, t); flush(io)
            @printf("epoch %3d  train %.3f  val %.3f  (%.1f s)\n", epoch, tl / nb, vl, t)
            if vl < best
                best, bestps, wait = vl, deepcopy(state.parameters), 0
            elseif (wait += 1) ≥ patience
                println("early stop: no improvement for $patience epochs"); break
            end
        end
    end

    ŷ = cpu_device()(first(model(Xd[:, ite], bestps, Lux.testmode(state.states))))
    pred = [names[j] in TIME_DEPENDENT ? ŷ[j, :] .* resc[ite] : ŷ[j, :] for j in 1:nout]
    m(a, b) = 100 * mean(abs.((a .- b) ./ a))
    println("test MAPE ($net $arch): ", join(["$(names[j]) $(round(m(T[j, ite], pred[j]), digits = 2))%" for j in 1:nout], ", "))
    open(out("test.csv"), "w") do io
        println(io, join(["tree"; [[n, "$(n)_pred"] for n in names]...], ','))
        for k in eachindex(ite)
            println(io, join([string(tree[ite[k]]); [string.([T[j, ite[k]], pred[j][k]]) for j in 1:nout]...], ','))
        end
    end
    jldsave(out("params.jld2"); ps = cpu_device()(bestps), st = cpu_device()(state.states), net, arch, opt, scaler)
end

## z-scores with the mean and population SD of the training columns `itr`
standardize(S, itr) = begin
    μ = mean(S[:, itr]; dims = 2)
    σ = std(S[:, itr]; dims = 2, corrected = false)
    σ[σ .== 0] .= 1                          # sklearn leaves constant columns unscaled
    Float32.((S .- μ) ./ σ), (mean = vec(μ), scale = vec(σ))
end

main(args) = begin
    dir = args[1]
    opt = merge(DEFAULTS, Dict(split(a, "=")[1] => split(a, "=")[2] for a in args[2:end]))
    archs = opt["arch"] == "both" ? ["released", "paper"] : [opt["arch"]]
    nets = opt["net"] == "both" ? ["cnn", "ffnn"] : [opt["net"]]
    all(in(("cnn", "ffnn")), nets) || error("net must be cnn, ffnn or both")
    all(in(keys(ARCHS)), archs) || error("arch must be released, paper, python or both")

    f = joinpath(dir, "cblv.h5")
    X = h5read(f, "cblv")
    T = Float32.(h5read(f, "targets"))
    names = h5open(ff -> haskey(ff, "target_names") ? read(ff["target_names"]) : ["R_0", "Infectious_Period"], f)
    resc = Float32.(h5read(f, "rescale"))
    Y = reduce(vcat, [names[j] in TIME_DEPENDENT ? T[j:j, :] ./ resc' : T[j:j, :] for j in eachindex(names)])
    N = size(X, 2)
    ntest, nval = min(10_000, N ÷ 20), min(10_000, N ÷ 10)
    itr, iva, ite = 1:(N - nval - ntest), (N - nval - ntest + 1):(N - ntest), (N - ntest + 1):N
    println("trees $N: train $(length(itr)), validation $(length(iva)), test $(length(ite)); input $(size(X, 1))")
    d = (; cblv = X, ss = h5read(f, "ss"), Y, T, resc, itr, iva, ite, dir, tree = h5read(f, "tree"), names)
    for net in nets, a in archs
        train(net, a, d, opt)
    end
end

main(ARGS)
