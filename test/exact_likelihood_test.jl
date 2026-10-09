"""
Exact likelihoods by dynamic programming over individual hosts (`scripts/exact_dp.jl`), which uses neither the
KLI algebra nor any filter, against `lbdp_exact` and against the generic filter.

1. `lbdp_dp` equals `lbdp_exact` on the LBDP trees of `mgp_filter_reference.tsv`.
2. `seir_dp` reproduces the stored L and L' of `exact_seir_reference.tsv` (N = 5 hosts), and L > L':
   the outcome "the lineage stays with the infector" has positive mass when every infectious host is tracked.
3. The generic SEIR filter (`mgp_filter_pomp`) is unbiased for L: the mean of exp(ll - L) over runs is 1
   within Monte Carlo error, and far from exp(L' - L).
"""
module ExactLikelihoodTest

import ..Main: h1, h2

@info h1("Exact likelihoods by dynamic programming")

using Test
using PhyloPOMP
using Random: seed!

include(joinpath(@__DIR__, "..", "scripts", "exact_dp.jl"))

rows(file) = [split(l, '\t') for l in eachline(joinpath(@__DIR__, file)) if !startswith(l, "#")]

@testset verbose=true "exact likelihoods" begin

    @testset "LBDP: dynamic program equals lbdp_exact" begin
        θ = (λ = 1.5, μ = 0.5, ψ = 0.3, χ = 0.2)
        for r in rows("mgp_filter_reference.tsv")
            r[1] == "lbdp" || continue
            g = parse_newick(String(r[3]); t0 = 0.0, time = parse(Float64, r[4]))
            @test lbdp_dp(g; θ...) ≈ lbdp_exact(g; θ...) atol = 1e-8
        end
    end

    θ = (β = 3.0, σ = 1.5, γ = 0.8, ω = 0.5, ψ = 0.6, χ = 0.0, N = 5.0)
    x0 = (S = 4, E = 0, I = 1, R = 0)
    trees = rows("exact_seir_reference.tsv")
    dp(g; kw...) = seir_dp(g; β = θ.β, σ = θ.σ, γ = θ.γ, ω = θ.ω, ψ = θ.ψ, N = 5, x0 = x0, kw...)

    @testset "SEIR: dynamic program, L and L'" begin
        for r in trees
            g = parse_newick(String(r[3]); t0 = 0.0, time = parse(Float64, r[4]))
            L, Lp = parse(Float64, r[5]), parse(Float64, r[6])
            @test dp(g) ≈ L atol = 1e-8
            @test dp(g; drop_stay = true) ≈ Lp atol = 1e-8
            @test L - Lp > 0.5
        end
    end

    @testset "SEIR: generic filter is unbiased for the exact likelihood" begin
        seed!(20261008)
        for r in trees[2:3]
            g = parse_newick(String(r[3]); t0 = 0.0, time = parse(Float64, r[4]))
            L, Lp = parse(Float64, r[5]), parse(Float64, r[6])
            P = mgp_filter_pomp(g, PhyloPOMP.SEIR; θ = θ, x0 = x0)
            w = [exp(logLik(pfilter(P, Np = 1000)) - L) for _ in 1:30]
            m = sum(w) / length(w)
            se = sqrt(sum((w .- m) .^ 2) / (length(w) - 1) / length(w))
            @info "tree $(r[2]): mean exp(ll - L) = $(round(m, digits = 3)) ± $(round(se, digits = 3)); exp(L' - L) = $(round(exp(Lp - L), digits = 3))"
            @test abs(m - 1) < 4se + 0.02
            @test m - exp(Lp - L) > 10se
        end
    end

    @testset "MERS: exact L, L', and the generic filters (naive, guided) are unbiased" begin
        DM = PhyloPOMP.GuidedMERS.Demes
        θm = (β_cc = 3.0, β_ch = 1.0, β_hc = 1.5, β_hh = 2.0, γ_c = 0.7, γ_h = 0.7, χ_c = 0.6, χ_h = 0.6,
              B_c = 0.0, B_h = 0.0, N_c = 3.0, N_h = 3.0)
        xm = (S_c = 2, I_c = 1, S_h = 3, I_h = 0)
        guide = fsmarkov(DM.Camel => 0.5, DM.Human => 0.5, (DM.Camel, DM.Human) => 0.05)
        seed!(20261010)
        for r in rows("exact_mers_reference.tsv")
            g = parse_newick(String(r[3]); demes = DM, t0 = 0.0, time = parse(Float64, r[4]))
            L, Lp = parse(Float64, r[5]), parse(Float64, r[6])
            @test mers_dp(g; θ = θm, Nc = 3, Nh = 3) ≈ L atol = 1e-8
            @test mers_dp(g; θ = θm, Nc = 3, Nh = 3, drop_stay = true) ≈ Lp atol = 1e-8
            @test L > Lp
            r[2] in ("1", "3") || continue
            for (name, kw) in (("naive", (;)), ("guided", (guide = guide,)))
                P = mgp_filter_pomp(g, PhyloPOMP.MERS; θ = θm, x0 = xm, demeset = DM, kw...)
                w = [exp(logLik(pfilter(P, Np = 1000)) - L) for _ in 1:30]
                m = sum(w) / length(w)
                se = sqrt(sum((w .- m) .^ 2) / (length(w) - 1) / length(w))
                @info "MERS $name, tree $(r[2]): mean exp(ll - L) = $(round(m, digits = 3)) ± $(round(se, digits = 3))"
                @test abs(m - 1) < 4se + 0.02
            end
        end
    end

    @testset "SEIR: guided generic filter is unbiased for the exact likelihood" begin
        D2 = PhyloPOMP.MGPDemes2
        guide = fsmarkov(D2.d1 => 0.1, D2.d2 => 1, (D2.d1, D2.d2) => 1)   # demes (E, I)
        seed!(20261009)
        for r in trees[[1, 3]]
            g = parse_newick(String(r[3]); t0 = 0.0, time = parse(Float64, r[4]))
            L = parse(Float64, r[5])
            P = mgp_filter_pomp(g, PhyloPOMP.SEIR; θ = θ, x0 = x0, guide = guide)
            w = [exp(logLik(pfilter(P, Np = 1000)) - L) for _ in 1:30]
            m = sum(w) / length(w)
            se = sqrt(sum((w .- m) .^ 2) / (length(w) - 1) / length(w))
            @info "guided, tree $(r[2]): mean exp(ll - L) = $(round(m, digits = 3)) ± $(round(se, digits = 3))"
            @test abs(m - 1) < 4se + 0.02
        end
    end

    @testset "SEIR: soft and hard kernels are unbiased for the exact likelihood" begin
        D2 = PhyloPOMP.MGPDemes2
        guide = fsmarkov(D2.d1 => 0.1, D2.d2 => 1, (D2.d1, D2.d2) => 1)
        seed!(20261011)
        r = trees[3]
        g = parse_newick(String(r[3]); t0 = 0.0, time = parse(Float64, r[4]))
        L = parse(Float64, r[5])
        for prop in (:soft, :hard)
            P = mgp_filter_pomp(g, PhyloPOMP.SEIR; θ = θ, x0 = x0, guide = guide, proposal = prop)
            w = [exp(logLik(pfilter(P, Np = 1000)) - L) for _ in 1:30]
            m = sum(w) / length(w)
            se = sqrt(sum((w .- m) .^ 2) / (length(w) - 1) / length(w))
            @info "$prop, tree $(r[2]): mean exp(ll - L) = $(round(m, digits = 3)) ± $(round(se, digits = 3))"
            # the weights are skewed (median below the mean), so allow for a low sample mean
            @test abs(m - 1) < 4se + 0.05
        end
        @test_throws ArgumentError mgp_filter_pomp(g, PhyloPOMP.SEIR; θ = θ, x0 = x0, guide = guide, proposal = :fancy)
    end
end

end # module ExactLikelihoodTest
