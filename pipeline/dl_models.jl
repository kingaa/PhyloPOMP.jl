# M6: the networks of Voznica et al. (2022) in Lux. CNN-CBLV reads the 402- or 1002-entry vector of
# encode_trees.jl; FFNN-SS reads the 99 summary statistics (98 + p_sample). Both output two values:
# R0 and the infectious period divided by the tree's rescale factor.
#
# Defaults reproduce phylodeep's pretrained BD models layer for layer (read from their .h5 configs):
# no dropout, ELU on every hidden layer. train_bd_models.py instead uses dropout 0.5 after the 64-, 32- and
# 16-unit layers and a linear 8-unit layer; pass `dropout = 0.5, last = identity` for that.

using Lux, HDF5      # Lux re-exports NNlib's `elu`

## Keras `Reshape((L ÷ 2, 2))` of a row-major (N, L) batch: pairs entries (1,2), (3,4), … into channels.
keras_pairs(v) = permutedims(reshape(v, 2, size(v, 1) ÷ 2, size(v, 2)), (2, 1, 3))

## The dense head shared by both networks: nin → 64 → 32 → 16 → 8 → nout (2 for BD, 3 for BDEI, 4 for BDSS).
head(nin, dropout, last, nout = 2) = dropout > 0 ?
    (Dense(nin => 64, elu), Dropout(dropout), Dense(64 => 32, elu), Dropout(dropout),
     Dense(32 => 16, elu), Dropout(dropout), Dense(16 => 8, last), Dense(8 => nout)) :
    (Dense(nin => 64, elu), Dense(64 => 32, elu), Dense(32 => 16, elu), Dense(16 => 8, last), Dense(8 => nout))

"""
    cnn_cblv(; dropout = 0.0, last = elu)

CNN-CBLV. Input: a (402 or 1002) × batch matrix from `phylodeep_cblv`.
"""
cnn_cblv(; dropout = 0.0, last = elu, nout = 2) = Chain(
    WrappedFunction(keras_pairs),
    Conv((3,), 2 => 50, elu; cross_correlation = true),
    Conv((10,), 50 => 50, elu; cross_correlation = true),
    MaxPool((10,)),
    Conv((10,), 50 => 80, elu; cross_correlation = true),
    GlobalMeanPool(),
    FlattenLayer(),
    head(80, dropout, last, nout)...,
)

"""
    ffnn_ss(; dropout = 0.0, last = elu, ninput = 99)

FFNN-SS. Input: ninput × batch, standardized as phylodeep's scaler does.
"""
ffnn_ss(; dropout = 0.0, last = elu, ninput = 99, nout = 2) = Chain(head(ninput, dropout, last, nout)...)

"""
    load_keras(model, ps, file) -> ps

Parameters of a Keras .h5 model (phylodeep's pretrained files) in the layout of `ps`, the parameters
of the matching `cnn_cblv()` or `ffnn_ss()` with no dropout. Conv kernels are (out, in, k) as HDF5.jl
reads them and become Lux's (k, in, out); Dense kernels are already (out, in).
"""
load_keras(model, ps, file) = h5open(file) do f
    w = f["model_weights"]
    names = filter(n -> startswith(n, "conv1d") || startswith(n, "dense"), read_attribute(w, "layer_names"))
    weighted = [k for k in keys(ps) if haskey(ps[k], :weight)]
    length(names) == length(weighted) || error("$(length(names)) Keras layers with weights, $(length(weighted)) in the model")
    new = Dict{Symbol,Any}(pairs(ps))
    for (k, n) in zip(weighted, names)
        W = read(w[n][n]["kernel:0"])
        b = read(w[n][n]["bias:0"])
        W = ndims(W) == 3 ? permutedims(W, (3, 2, 1)) : W
        size(W) == size(ps[k].weight) || error("$n: kernel $(size(W)) vs $(size(ps[k].weight))")
        new[k] = (weight = Float32.(W), bias = Float32.(b))
    end
    NamedTuple{keys(ps)}(Tuple(new[k] for k in keys(ps)))
end
