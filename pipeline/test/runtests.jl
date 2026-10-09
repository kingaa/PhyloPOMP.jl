# Pipeline checks against phylodeep (Python reference only).
#   julia --project=pipeline pipeline/test/runtests.jl
#
# data/bd_{small,large}_trees.tsv: 25 and 8 trees from generate_trees.jl (seeds 7 and 11).
# data/phylodeep_cblv_{small,large}.csv: phylodeep's encode_into_most_recent on the same trees
#   (scripts/phylodeep_cblv_oracle.py); CBLV columns, then the rescale factor.
# data/phylodeep_ss_{small,large}.csv: encode_into_summary_statistics (scripts/phylodeep_ss_oracle.py), 99 values + rescale.
# data/phylodeep_bd_small_ffnn_{scaler,output}.csv: BD_SMALL_FFNN's scaler (mean; scale) and Keras outputs on the small rows.
# data/phylodeep_bd_small_cnn_output.csv: Keras outputs of phylodeep's pretrained BD_SMALL_CNN on the small rows
#   (scripts/phylodeep_cnn_oracle.py). The CNN check needs that .h5 (PHYLODEEP_MODELS, or the venv default).

using Test, PhyloPOMP, Lux, Random
include(joinpath(@__DIR__, "..", "encode_trees.jl"))     # also includes sumstats.jl
include(joinpath(@__DIR__, "..", "dl_models.jl"))

data(f) = joinpath(@__DIR__, "data", f)
trees(f) = [split(l, '\t') for l in Iterators.drop(eachline(data(f)), 1)]
csv(f) = [parse.(Float64, split(l, ',')) for l in eachline(data(f))]

@testset "pipeline vs phylodeep" begin
    @testset "CBLV encoding ($size)" for size in ("small", "large")
        rows, ref = trees("bd_$(size)_trees.tsv"), csv("phylodeep_cblv_$size.csv")
        @test length(rows) == length(ref)
        for (r, o) in zip(rows, ref)
            v, s = phylodeep_cblv(parse_newick(String(r[7])), parse(Float64, r[4]))
            @test length(v) == length(o) - 1 == (size == "small" ? 402 : 1002)
            @test v ≈ o[1:end-1] atol = 1e-10
            @test s ≈ o[end] rtol = 1e-12
        end
    end

    @testset "summary statistics ($size)" for size in ("small", "large")
        rows, ref = trees("bd_$(size)_trees.tsv"), csv("phylodeep_ss_$size.csv")
        @test length(rows) == length(ref)
        for (r, o) in zip(rows, ref)
            v, s = phylodeep_sumstats(parse_newick(String(r[7])), parse(Float64, r[4]))
            @test length(v) == 99
            @test v ≈ o[1:end-1] rtol = 1e-10
            @test s ≈ o[end] rtol = 1e-12
        end
    end

    models = get(ENV, "PHYLODEEP_MODELS", joinpath(homedir(),
        "Desktop/Projects/efficiency/DeepLearningPhylodeep/venv/lib/python3.12/site-packages/phylodeep/pretrained_models/models"))
    w = joinpath(models, "BD_SMALL_CNN.h5")
    @testset "Lux CNN-CBLV with phylodeep weights" begin
        if isfile(w)
            X = Float32.(reduce(hcat, [phylodeep_cblv(parse_newick(String(r[7])), parse(Float64, r[4]))[1]
                                      for r in trees("bd_small_trees.tsv")]))
            model = cnn_cblv()
            ps, st = Lux.setup(Xoshiro(0), model)
            y, _ = model(X, load_keras(model, ps, w), Lux.testmode(st))
            ref = reduce(hcat, csv("phylodeep_bd_small_cnn_output.csv"))
            @test y ≈ ref rtol = 1e-5
        else
            @info "$w not found; skipped the CNN check"
        end
    end
    wf = joinpath(models, "BD_SMALL_FFNN.h5")
    @testset "Lux FFNN-SS with phylodeep weights and scaler" begin
        if isfile(wf)
            X = reduce(hcat, [phylodeep_sumstats(parse_newick(String(r[7])), parse(Float64, r[4]))[1]
                              for r in trees("bd_small_trees.tsv")])
            μ, σ = csv("phylodeep_bd_small_ffnn_scaler.csv")
            model = ffnn_ss()
            ps, st = Lux.setup(Xoshiro(0), model)
            y, _ = model(Float32.((X .- μ) ./ σ), load_keras(model, ps, wf), Lux.testmode(st))
            @test y ≈ reduce(hcat, csv("phylodeep_bd_small_ffnn_output.csv")) rtol = 1e-5
        else
            @info "$wf not found; skipped the FFNN check"
        end
    end
end
