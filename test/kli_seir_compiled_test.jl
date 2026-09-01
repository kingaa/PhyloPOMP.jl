"""
M08 acceptance-gate test: Gate 5 (numerical equivalence) for the compiled
SEIR filter (`src/examples/mgp_seir_filter.jl`) against the trusted,
hand-coded `NaiveSEIR.filter_pomp` (`src/examples/seir_naive.jl`, read-only
oracle, never modified).

Approach (per the M08 task's option (b)): both filters are single-particle
(`Np=1`) importance samplers whose log-likelihood is itself a random variable
(depends on the specific silent/untracked regular events the proposal
happens to draw). Exact numerical equivalence is only meaningful if BOTH
filters draw the IDENTICAL sequence of random numbers for the IDENTICAL
genealogy/parameters -- so each comparison resets the GLOBAL RNG
(`Random.seed!(seed)`) to the SAME seed immediately before each `pfilter`
call. `compiled_regular_part!` was written to make the exact same `rcateg`/
`rand()` calls, in the same order, as `NaiveSEIR.regular_part!`
(`src/examples/mgp_seir_filter.jl`'s own docstring), so this is a legitimate
bit-for-bit comparison, not a statistical one.

Over MANY (parameter, genealogy, seed) combinations, every PAIR of finite
log-likelihoods must agree to a tight tolerance (`atol=1e-6, rtol=1e-8` --
loose enough to absorb floating-point summation-order noise from computing
some quantities via `Rational{Int}` Φ_u then converting to `Float64` at a
different point in the expression than `seir_naive.jl`'s own arithmetic, but
tight enough that a REAL formula discrepancy -- e.g. the M07 decay
discrepancy this milestone resolved, or the M03 Q_u regular/singular gating
bug this milestone found and fixed -- would fail it immediately, as both did
during development).
"""
module KliSeirCompiledTest

import ..Main: h1, h2

@info h1("Compiled SEIR filter vs. seir_naive.jl (M08 Gate 5)")

using Test
using PhyloPOMP
using PhyloPOMP.NaiveSEIR
using Random: MersenneTwister, seed!
import PartiallyObservedMarkovProcesses as POMP

function build_genealogy(θ, x0, target_n, tmax, gseed; attempts = 800)
    rng = MersenneTwister(gseed)
    g = nothing
    for _ in 1:attempts
        g = simulate(PhyloPOMP.SEIR, θ; x0 = x0, graft = [0, 1], tmax = tmax, rng = rng)
        nsample(g) == target_n && return g, true
    end
    g, false
end

function compare_ll(g, β, σ, γ, ω, ψ, χ, pop, x0; nseeds)
    p_naive = NaiveSEIR.filter_pomp(
        g; β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, pop = pop,
        S0 = x0.S / pop, E0 = x0.E / pop, I0 = x0.I / pop, R0 = x0.R / pop,
    )
    p_comp = PhyloPOMP.compiled_filter_pomp(
        g; β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, pop = pop,
        S0 = x0.S / pop, E0 = x0.E / pop, I0 = x0.I / pop, R0 = x0.R / pop,
    )
    nfinite = 0
    nmatch = 0
    worst = 0.0
    for s in 1:nseeds
        seed!(s)
        ll1 = logLik(pfilter(p_naive, Np = 1))
        seed!(s)
        ll2 = logLik(pfilter(p_comp, Np = 1))
        if isfinite(ll1) || isfinite(ll2)
            nfinite += 1
            ok = isapprox(ll1, ll2; atol = 1e-6, rtol = 1e-8)
            ok && (nmatch += 1)
            worst = max(worst, abs(ll1 - ll2))
        end
    end
    nfinite, nmatch, worst
end

@testset verbose=true "Compiled SEIR filter (M08)" begin

    @info h2("Isolated regular-part unit check: many (state, ell) draws, " *
             "mixed event types, bit-exact against NaiveSEIR.regular_part!")
    @testset "compiled_regular_part! isolation" begin
        rng = MersenneTwister(20260818)
        nmismatch = 0
        ntried = 0
        for _ in 1:400
            β = 0.5 + 4 * rand(rng); σ = 0.2 + 2 * rand(rng)
            γ = 0.2 + 2 * rand(rng); ω = 0.1 + 1.5 * rand(rng)
            S = rand(rng, 5:80); E = rand(rng, 1:15); I = rand(rng, 1:15); R = rand(rng, 0:10)
            ellE = rand(rng, 0:E); ellI = rand(rng, 0:I)
            seed = rand(rng, 1:10^7)

            cols_n = PhyloPOMP.Coloring(NaiveSEIR.Demes)
            cols_c = PhyloPOMP.Coloring(NaiveSEIR.Demes)
            lin = 1
            for _ in 1:ellE
                push!(cols_n[NaiveSEIR.Expos], lin); push!(cols_c[NaiveSEIR.Expos], lin); lin += 1
            end
            for _ in 1:ellI
                push!(cols_n[NaiveSEIR.Infec], lin); push!(cols_c[NaiveSEIR.Infec], lin); lin += 1
            end

            seed!(seed)
            ll_n, Sn, En, In, Rn = NaiveSEIR.regular_part!(
                cols_n, 0.0, 0.0, 3.0, S, E, I, R;
                β = β, σ = σ, γ = γ, ω = ω, ψ = 0.0, χ = 0.0, pop = 100.0,
            )
            seed!(seed)
            ll_c, Sc, Ec, Ic, Rc = PhyloPOMP.compiled_regular_part!(
                cols_c, 0.0, 0.0, 3.0, S, E, I, R;
                β = β, σ = σ, γ = γ, ω = ω, ψ = 0.0, χ = 0.0, pop = 100.0, model = PhyloPOMP.SEIR,
            )
            ntried += 1
            ok = isapprox(ll_n, ll_c; atol = 1e-8, rtol = 1e-10) &&
                 (Sn, En, In, Rn) == (Sc, Ec, Ic, Rc) && ell(cols_n) == ell(cols_c)
            ok || (nmismatch += 1)
        end
        @info "compiled_regular_part! isolation: $ntried trials, $nmismatch mismatches"
        @test nmismatch == 0
    end

    @info h2("Gate 5: end-to-end log-likelihood, many (params, genealogy, " *
             "seed) combinations, seed-matched Np=1 particle filters")
    @testset "end-to-end log-likelihood sweep" begin
        rng_master = MersenneTwister(20260819)
        total_finite = 0
        total_match = 0
        ncombos = 0
        worst_abs = 0.0
        for combo in 1:20
            β  = 1.0 + 6.0 * rand(rng_master)
            σ  = 0.3 + 3.0 * rand(rng_master)
            γ  = 0.3 + 3.0 * rand(rng_master)
            ω  = 0.1 + 2.0 * rand(rng_master)
            ψ  = 0.02 + 0.3 * rand(rng_master)
            χ  = 0.0
            pop = rand(rng_master, (60, 100, 150))
            x0 = (S = pop - 1, E = 0, I = 1, R = 0)
            θ = (β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, N = Float64(pop))
            target_n = rand(rng_master, 1:4)
            tmax = 2.0 + 4.0 * rand(rng_master)

            g, found = build_genealogy(θ, x0, target_n, tmax, hash((:m08, combo)))
            found || continue
            ncombos += 1

            nfinite, nmatch, worst = compare_ll(g, β, σ, γ, ω, ψ, χ, pop, x0; nseeds = 150)
            total_finite += nfinite
            total_match += nmatch
            worst_abs = max(worst_abs, worst)
            @test nfinite == nmatch
        end
        @info "Gate 5 sweep: $ncombos parameter/genealogy combos, " *
              "$total_finite finite log-likelihood comparisons, " *
              "$total_match exact matches, worst |Δll|=$worst_abs"
        @test ncombos ≥ 15   # most combos should successfully build a genealogy
        @test total_finite ≥ 50
    end

end

end # module KliSeirCompiledTest
