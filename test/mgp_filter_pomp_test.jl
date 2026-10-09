"""
`mgp_filter_pomp` (generic singular and regular steps) on the LBDP, SIR, BDSS and SI2R tables, against
values from R phylopomp stored in `mgp_filter_reference.tsv` (no R needed here). The SI2R values come from
upstream phylopomp 0.19.8.1 (see the header of that file).

1. The new tables LBDP, BDEI and BDSS pass `validate_model` and simulate.
2. The Julia port of `lbdp_exact` equals R's value on the same genealogies.
3. The generic LBDP filter estimates the exact log-likelihood.
4. The generic SIR, BDSS and SI2R filters agree with R's `sir_pomp`, `bdss_pomp` and upstream `si2rs_pomp`.
5. Argument checks of `mgp_filter_pomp` and `generic_demeset`.
"""
module MgpFilterPompTest

import ..Main: h1, h2

@info h1("Generic filter pomp object vs. R phylopomp values")

using Test
using PhyloPOMP
using PhyloPOMP: validate_model, Unstructured
using Random: MersenneTwister, seed!

const REF = [split(l, '\t') for l in eachline(joinpath(@__DIR__, "mgp_filter_reference.tsv")) if !startswith(l, "#")]

const SETUP = Dict(
    "lbdp" => (LBDP, (λ = 1.5, μ = 0.5, ψ = 0.3, χ = 0.2), (n = 1,)),
    "sir"  => (SIR, (β = 3.0, γ = 1.0, ψ = 0.3, N = 100.0), (S = 99, I = 1, R = 0)),
    "bdss" => (BDSS, (λ_nn = 1.0, λ_ns = 0.3, λ_sn = 1.5, λ_ss = 2.5, μ = 0.5, χ = 0.4), (N = 1, S = 0)),
    "si2r" => (SI2R, (β = 4.0, κ = 3.0, γ = 1.0, ω = 0.5, ψ = 0.3, η_L = 0.5, η_H = 1.0, N = 100.0, r = 1.0),
               (S = 99, I_L = 1, I_H = 0, R = 0)),
)

function estimate(g, model, θ, x0; Np = 1000, nrep = 10)
    p = mgp_filter_pomp(g, model; θ = θ, x0 = x0)
    logmeanexp([logLik(pfilter(p, Np = Np)) for _ in 1:nrep], se = true)
end

@testset verbose=true "generic filter pomp object" begin

    @testset "new tables" begin
        for (M, θ, x0, graft) in ((LBDP, SETUP["lbdp"][2], (n = 1,), [1]),
                                  (BDEI, (σ = 1.0, λ = 2.0, μ = 0.5, χ = 0.4), (E = 0, I = 1), [0, 1]),
                                  (BDSS, SETUP["bdss"][2], (N = 1, S = 0), [1, 0]))
            @test isempty(validate_model(M))
            g = simulate(M, θ; x0 = x0, graft = graft, tmax = 2.0, rng = MersenneTwister(1))
            @test g isa Genealogy
        end
    end

    @testset "lbdp_exact port equals R" begin
        for f in REF
            f[1] == "lbdp" || continue
            g = parse_newick(String(f[3]); t0 = 0.0, time = parse(Float64, f[4]))
            @test lbdp_exact(g; SETUP["lbdp"][2]...) ≈ parse(Float64, f[7]) atol = 1e-10
        end
    end

    @testset "filters against the exact likelihood and R" begin
        seed!(20261007)
        for f in REF
            M, θ, x0 = SETUP[f[1]]
            g = parse_newick(String(f[3]); t0 = 0.0, time = parse(Float64, f[4]))
            est, se = estimate(g, M, θ, x0)
            if f[1] == "lbdp"
                ref, rse = parse(Float64, f[7]), 0.0
            else
                ref, rse = parse(Float64, f[5]), parse(Float64, f[6])
            end
            z = (est - ref) / sqrt(se^2 + rse^2)
            @info "$(f[1]) tree $(f[2]): Julia $(round(est, digits = 3)) ± $(round(se, digits = 3)), " *
                  "reference $(round(ref, digits = 3)), z = $(round(z, digits = 2))"
            @test abs(z) < 4
        end
    end

    @testset "argument checks" begin
        g = parse_newick(String(REF[1][3]); t0 = 0.0, time = parse(Float64, REF[1][4]))
        @test_throws ArgumentError mgp_filter_pomp(g, LBDP; θ = SETUP["lbdp"][2], x0 = (m = 1,))
        @test_throws ArgumentError mgp_filter_pomp(g, LBDP; θ = SETUP["lbdp"][2], x0 = (n = -1,))
        @test_throws ArgumentError mgp_filter_pomp(g, LBDP; θ = SETUP["lbdp"][2], x0 = (n = 1.5,))
        @test generic_demeset(LBDP) === Unstructured
        @test length(instances(generic_demeset(BDEI).DemeSet)) == 2
        @test length(instances(generic_demeset(PhyloPOMP.MERS).DemeSet)) == 2
        @test_throws ArgumentError mgp_filter_pomp(g, LBDP; θ = SETUP["lbdp"][2], x0 = (n = 2,), maxpop = 1)
    end

    @testset "guide knowledge from the table" begin
        know(M, type) = (v = zeros(2); k = PhyloPOMP._model_knowledge(M, generic_demeset(M).DemeSet);
                         (k(v; deme = missing, type = type, time = 0.0), v))
        ## SEIR and BDEI: samples and branch points in I (the second deme); BDSS: samples and births from both types
        for M in (PhyloPOMP.SEIR, BDEI)
            @test know(M, PhyloPOMP.Sample) == (true, [0.0, 1.0])
            @test know(M, PhyloPOMP.Node) == (true, [0.0, 1.0])
            @test know(M, PhyloPOMP.Root)[1] == false
        end
        @test know(BDSS, PhyloPOMP.Sample)[1] == false
        @test know(BDSS, PhyloPOMP.Node)[1] == false
        ## on a tree without demes, the guide's weights are no longer all 1 (with the default knowledge they were)
        D = PhyloPOMP.MGPDemes2
        g = let rng = MersenneTwister(20260813), g = nothing
            for _ in 1:100
                g = simulate(PhyloPOMP.SEIR, (β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.1, χ = 0.0, N = 100.0);
                             x0 = (S = 99, E = 0, I = 1, R = 0), graft = [0, 1], tmax = 20.0, rng = rng)
                nsample(g) >= 20 && break
            end
            g
        end
        fs = fsmarkov(D.d1 => 0.1, D.d2 => 1, (D.d1, D.d2) => 1)
        weights(gd) = reduce(vcat, [PhyloPOMP.relhaz((gd[k].tbeg + gd[k].tend) / 2, gd, k, D.d2, D.d1,
                                                      collect(1:length(gd[k].alllins)))
                                    for k in 1:length(gd.nodes) if !isempty(gd[k].alllins)])
        @test all(isapprox.(weights(PhyloPOMP.Guide(g, fs)), 1.0; atol = 1e-9))
        @test !all(isapprox.(weights(PhyloPOMP.Guide(g, fs, PhyloPOMP._model_knowledge(PhyloPOMP.SEIR, D.DemeSet))), 1.0;
                             atol = 1e-9))
    end

    @testset "population cap (maxpop)" begin
        ## a 30-tip LBDP tree, filtered at a growth rate r = 4 that, uncapped, would simulate about exp(4·span) hosts
        g = let rng = MersenneTwister(11), g = nothing
            for _ in 1:200
                g = simulate(LBDP, (λ = 1.0, μ = 0.2, ψ = 0.0, χ = 0.3); x0 = (n = 1,), graft = [1], tmax = 12.0, rng = rng)
                nsample(g) >= 30 && break
            end
            g
        end
        @test nsample(g) >= 30
        fast = (λ = 5.0, μ = 0.5, ψ = 0.0, χ = 0.5)
        Gcap = mgp_filter_pomp(g, LBDP; θ = fast, x0 = (n = 1,), maxpop = 100 * nsample(g))
        t = @elapsed ll = logLik(pfilter(Gcap, Np = 50))
        @test ll == -Inf
        @test t < 60
        ## a cap far above the plausible counts leaves the estimate as without a cap: same random numbers, same value
        θ = (λ = 1.0, μ = 0.2, ψ = 0.0, χ = 0.3)
        seed!(5); a = logLik(pfilter(mgp_filter_pomp(g, LBDP; θ = θ, x0 = (n = 1,)), Np = 200))
        seed!(5); b = logLik(pfilter(mgp_filter_pomp(g, LBDP; θ = θ, x0 = (n = 1,), maxpop = 10^9), Np = 200))
        @test a === b
    end
end

end # module MgpFilterPompTest
