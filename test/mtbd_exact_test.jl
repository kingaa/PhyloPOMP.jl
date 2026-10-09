"""
`mtbd_exact`: the exact likelihood of the linear tables (LBDP, BDEI, BDSS, MTBD) from the multi-type birth-death
equations, and `lbdp_exact` on genealogies long enough to overflow cosh and sinh.

1. LBDP: `mtbd_exact` equals the closed form `lbdp_exact` (n0 = 1), sampled ancestors included (ψ > 0).
2. BDSS with X = 1 (every host infects at rate B whatever its type) equals `lbdp_exact` with λ = B, for either
   founder type and for their mixture.
3. `lbdp_exact` is finite where d(tf - t)/2 passes 710, and equals `mtbd_exact` there.
4. The generic naive filter is unbiased for exp(mtbd_exact) on small BDEI and BDSS genealogies.
5. `mtbd_rates` reads the BDEI table; a table with compartments outside its demes, and a bad founder, throw.
"""
module MTBDExactTest

import ..Main: h1, h2

@info h1("Exact likelihood of linear tables (mtbd_exact)")

using Test
using PhyloPOMP
using Random: Xoshiro, seed!

## a simulated genealogy with between `nmin` and `nmax` samples
function sim_tree(model, θ, x0, graft, tmax, nmin, nmax, rng)
    for _ in 1:10_000
        g = simulate(model, θ; x0, graft, tmax, rng)
        nmin <= nsample(g) <= nmax && return g
    end
    error("no genealogy with $nmin to $nmax samples")
end

@testset verbose=true "mtbd_exact" begin

    @testset "LBDP: equals lbdp_exact, with sampled ancestors" begin
        rng = Xoshiro(1)
        for _ in 1:10
            θ = (λ = 1.5 + rand(rng), μ = 0.3rand(rng), ψ = 0.2 + 0.5rand(rng), χ = 0.3rand(rng))
            g = sim_tree(PhyloPOMP.LBDP, θ, (n = 1,), [1], 3.0, 3, 500, rng)
            @test mtbd_exact(g, PhyloPOMP.LBDP; θ, founder = :n) ≈ lbdp_exact(g; θ..., n0 = 1) rtol = 1e-8
        end
    end

    @testset "BDSS with X = 1 is LBDP with λ = B" begin
        rng = Xoshiro(2)
        B, f, μ, χ = 2.0, 0.15, 0.4, 0.6
        θ = (λ_nn = (1 - f) * B, λ_ns = f * B, λ_sn = (1 - f) * B, λ_ss = f * B, μ, χ)
        for _ in 1:3
            g = sim_tree(PhyloPOMP.BDSS, θ, (N = 1, S = 0), [1, 0], 4.0, 5, 500, rng)
            L = lbdp_exact(g; λ = B, μ, ψ = 0.0, χ, n0 = 1)
            for founder in (:N, :S, [1 - f, f])
                @test mtbd_exact(g, PhyloPOMP.BDSS; θ, founder) ≈ L rtol = 1e-8
            end
        end
    end

    @testset "lbdp_exact on a genealogy where cosh would overflow" begin
        rng = Xoshiro(3)
        θ = (λ = 1.2, μ = 0.2, ψ = 0.0, χ = 0.5)
        g = sim_tree(PhyloPOMP.LBDP, θ, (n = 1,), [1], 6.0, 10, 200, rng)
        θ300 = map(r -> 300r, θ)          # d(tf - t0)/2 near 1500
        L = lbdp_exact(g; θ300..., n0 = 1)
        @test isfinite(L)
        @test L ≈ mtbd_exact(g, PhyloPOMP.LBDP; θ = θ300, founder = :n) rtol = 1e-8
    end

    @testset "generic filter is unbiased for exp(mtbd_exact): BDEI, BDSS" begin
        seed!(20261009)
        rng = Xoshiro(4)
        cases = [(PhyloPOMP.BDEI, (σ = 0.8, λ = 1.5, μ = 0.3, χ = 0.4), (E = 0, I = 1), [0, 1], :I),
                 (PhyloPOMP.BDSS, (λ_nn = 1.0, λ_ns = 0.2, λ_sn = 4.0, λ_ss = 0.8, μ = 0.4, χ = 0.5),
                  (N = 1, S = 0), [1, 0], :N)]
        for (M, θ, x0, graft, founder) in cases
            g = sim_tree(M, θ, x0, graft, 4.0, 4, 10, rng)
            L = mtbd_exact(g, M; θ, founder)
            P = mgp_filter_pomp(g, M; θ, x0, maxpop = 2000)
            w = [exp(logLik(pfilter(P, Np = 2000)) - L) for _ in 1:30]
            m = sum(w) / length(w)
            se = sqrt(sum((w .- m) .^ 2) / (length(w) - 1) / length(w))
            @info "$(M.name), $(nsample(g)) samples: mean exp(ll - L) = $(round(m, digits = 3)) ± $(round(se, digits = 3))"
            @test abs(m - 1) < 4se + 0.02
        end
    end

    @testset "mtbd_rates and argument checks" begin
        R = mtbd_rates(PhyloPOMP.BDEI, (σ = 0.5, λ = 1.5, μ = 0.3, χ = 0.4))
        @test R.b == [0.0 0.0; 1.5 0.0]          # an I host infects a new E host
        @test R.m == [0.0 0.5; 0.0 0.0]          # E becomes I
        @test R.d == [0.0, 0.3] && R.χ == [0.0, 0.4] && R.ψ == [0.0, 0.0]
        θseir = (β = 3.0, σ = 1.5, γ = 0.8, ω = 0.5, ψ = 0.6, χ = 0.0, N = 5.0)
        @test_throws ArgumentError mtbd_rates(PhyloPOMP.SEIR, θseir)
        θ = (λ = 1.5, μ = 0.3, ψ = 0.0, χ = 0.5)
        g = sim_tree(PhyloPOMP.LBDP, θ, (n = 1,), [1], 3.0, 3, 50, Xoshiro(5))
        @test_throws ArgumentError mtbd_exact(g, PhyloPOMP.LBDP; θ, founder = :x)
        θss = (λ_nn = 1.0, λ_ns = 0.2, λ_sn = 4.0, λ_ss = 0.8, μ = 0.4, χ = 0.5)
        @test_throws ArgumentError mtbd_exact(g, PhyloPOMP.BDSS; θ = θss, founder = [0.5, 0.6])
    end
end

end
