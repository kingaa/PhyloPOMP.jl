"""
Generic singular_update! for MERS.

1. Exact Rational identities of the regular-move weights over an (ℓ, n) grid.
   These are the terms in Aaron King's `mers_guided.jl` (branch f501185 of
   ../PhyloPOMP.jl_aarons): `1 - C(ℓ,2)/C(I,2)` for within-host no-fork,
   `1 - ℓ_h/I_h` for no-move, `(1 - ℓ_c/I_c)/I_h` for cross.
2. Per-outcome Δll of singular_update! against the closed forms in
   `NaiveMERS.singular_part!` and `GuidedMERS.singular_branch!`/`terminal_sample!`/`singular_root!`.
3. Seed-matched Np=200 log-likelihoods of the compiled filter with generic and with
   NaiveMERS singular parts agree on simulated genealogies.
4. The same for SEIR against `NaiveSEIR.singular_part!` (and `GuidedSEIR` terminal_sample!,
   inline_sample!, singular_branch!): fork, destructive and non-destructive tip, sampled ancestor.
"""
module KliSingularUpdateTest

import ..Main: h1, h2

@info h1("Generic singular_update! (MERS)")

using Test
using PhyloPOMP
using PhyloPOMP.NaiveMERS
using PhyloPOMP: Coloring, ell, plant!, Root, Node, Sample, MERS, singular_update!,
    full_transitions, IdentityTransition, InlineSameDemeTransition, CrossDemeTransition
using PhyloPOMP.NaiveSEIR
using Random: MersenneTwister, seed!

phi_of(ts, T) = (i = findfirst(t -> t isa T, ts); isnothing(i) ? 0 // 1 : ts[i].phi)
find_event(name) = MERS.events[findfirst(e -> e.name == name, MERS.events)]
c2(x) = x * (x - 1) // 2

struct MockNode
    name::Int
    type::PhyloPOMP.NodeType
    lineage::Int
    children::Vector{Int}
    deme
end

D = NaiveMERS.Demes.DemeSet
Camel, Human = NaiveMERS.Camel, NaiveMERS.Human

θ = (β_cc = 3.0, β_ch = 0.8, β_hc = 1.2, β_hh = 2.5, γ_c = 1.0, γ_h = 1.0,
     χ_c = 0.5, χ_h = 0.4, B_c = 0.0, B_h = 0.0, N_c = 60.0, N_h = 50.0)

@testset verbose=true "generic singular_update! (MERS)" begin

    @testset "regular-move weights, exact (ℓ, n) grid" begin
        bad = 0; total = 0
        for Ic in 1:8, Ih in 1:8, lc in 0:Ic, lh in 0:Ih
            if Ic ≥ 2
                ts = full_transitions(find_event(:transmission_cc), [lc, lh], [Ic, Ih])
                total += 1
                bad += !(phi_of(ts, IdentityTransition) + lc * phi_of(ts, InlineSameDemeTransition) ==
                         1 - c2(lc) // c2(Ic))
            end
            if Ih ≥ 2
                ts = full_transitions(find_event(:transmission_hh), [lc, lh], [Ic, Ih])
                total += 1
                bad += !(phi_of(ts, IdentityTransition) + lh * phi_of(ts, InlineSameDemeTransition) ==
                         1 - c2(lh) // c2(Ih))
            end
            ts = full_transitions(find_event(:transmission_hc), [lc, lh], [Ic, Ih])
            total += 1
            bad += !(phi_of(ts, IdentityTransition) + lc * phi_of(ts, InlineSameDemeTransition) == 1 - lh // Ih)
            lh ≥ 1 && (total += 1; bad += !(phi_of(ts, CrossDemeTransition) == (1 - lc // Ic) // Ih))
            ts = full_transitions(find_event(:transmission_ch), [lc, lh], [Ic, Ih])
            total += 1
            bad += !(phi_of(ts, IdentityTransition) + lh * phi_of(ts, InlineSameDemeTransition) == 1 - lc // Ic)
            lc ≥ 1 && (total += 1; bad += !(phi_of(ts, CrossDemeTransition) == (1 - lh // Ih) // Ic))
        end
        @info "regular-move identities: $total checks, $bad failures"
        @test total > 10_000
        @test bad == 0
    end

    @testset "per-outcome Δll against the closed forms" begin
        x = (S_c = 40, I_c = 7, S_h = 30, I_h = 5)
        # parent camel lineage 1 with children 2, 3; node 1 is a fork
        geneal = Dict(
            1 => MockNode(1, Node, 1, [2, 3], missing),
            2 => MockNode(2, Sample, 2, Int[], Camel),
            3 => MockNode(3, Sample, 3, Int[], Human),
        )
        λcc = θ.β_cc * x.S_c * x.I_c / θ.N_c
        λhc = θ.β_hc * x.S_h * x.I_c / θ.N_c
        seen = Dict{Tuple,Float64}()
        for s in 1:400
            seed!(s)
            cols = Coloring(NaiveMERS.Demes)
            plant!(cols, Camel, 1); plant!(cols, Camel, 4)   # ℓ = (2, 0), one other camel lineage
            Δ, x1 = singular_update!(cols, geneal, 1, x, θ, MERS)
            kids = (2 ∈ cols[Camel] ? :C : :H, 3 ∈ cols[Camel] ? :C : :H)
            # proposal probability of the outcome
            if kids == (:C, :C)
                p = λcc / (λcc + λhc)
                Δ0 = log(λcc) - log((x.I_c + 1) * x.I_c / 2) - log(p)
                @test x1 == (S_c = 39, I_c = 8, S_h = 30, I_h = 5)
            else
                p = 0.5λhc / (λcc + λhc)
                Δ0 = log(λhc) - log((x.I_c) * (x.I_h + 1)) - log(p)
                @test x1 == (S_c = 40, I_c = 7, S_h = 29, I_h = 6)
            end
            @test Δ ≈ Δ0 atol = 1e-12
            seen[kids] = Δ
        end
        @test length(seen) == 3 || length(seen) == 4   # (C,C), (C,H), (H,C) reachable; (H,H) is not
        @test !haskey(seen, (:H, :H))

        # human parent
        geneal[1] = MockNode(1, Node, 1, [2, 3], missing)
        λhh = θ.β_hh * x.S_h * x.I_h / θ.N_h
        λch = θ.β_ch * x.S_c * x.I_h / θ.N_h
        for s in 1:200
            seed!(s)
            cols = Coloring(NaiveMERS.Demes)
            plant!(cols, Human, 1)
            Δ, x1 = singular_update!(cols, geneal, 1, x, θ, MERS)
            kids = (2 ∈ cols[Camel] ? :C : :H, 3 ∈ cols[Camel] ? :C : :H)
            if kids == (:H, :H)
                p = λhh / (λhh + λch)
                Δ0 = log(λhh) - log((x.I_h + 1) * x.I_h / 2) - log(p)
            else
                p = 0.5λch / (λhh + λch)
                Δ0 = log(λch) - log((x.I_c + 1) * (x.I_h)) - log(p)
            end
            @test Δ ≈ Δ0 atol = 1e-12
        end

        # sample: charge log(χ I) before the decrement
        geneal[2] = MockNode(2, Sample, 2, Int[], Camel)
        cols = Coloring(NaiveMERS.Demes); plant!(cols, Camel, 2)
        Δ, x1 = singular_update!(cols, geneal, 2, x, θ, MERS)
        @test Δ ≈ log(θ.χ_c * x.I_c)
        @test x1 == (S_c = 40, I_c = 6, S_h = 30, I_h = 5)
        @test 2 ∉ cols[Camel]
        # sample whose lineage is in the wrong deme: -Inf, state net unchanged
        cols = Coloring(NaiveMERS.Demes); plant!(cols, Human, 2)
        Δ, x1 = singular_update!(cols, geneal, 2, x, θ, MERS)
        @test Δ == -Inf
        @test x1 == x
        @test 2 ∉ cols[Camel] && 2 ∉ cols[Human]

        # root: charge -log p, p = (n_d - ℓ_d)/Σ(n - ℓ)
        geneal[1] = MockNode(1, Root, 1, [2], missing)
        cols = Coloring(NaiveMERS.Demes); plant!(cols, Camel, 9)
        for s in 1:50
            seed!(s)
            c = copy(cols)
            Δ, x1 = singular_update!(c, geneal, 1, x, θ, MERS)
            w = [x.I_c - 1, x.I_h]
            d = 1 ∈ c[Camel] ? 1 : 2
            @test Δ ≈ -log(w[d] / sum(w))
            @test x1 == x
        end
    end

    @testset "seed-matched log-likelihood, generic vs NaiveMERS singular part" begin
        function build(θ, x0, target_n, tmax, gseed)
            rng = MersenneTwister(gseed)
            for _ in 1:800
                g = simulate(MERS, θ; x0 = x0, graft = [1, 0], tmax = tmax, rng = rng,
                             demeset = NaiveMERS.Demes, samplemap = [Camel, Human])
                nsample(g) == target_n && return g
            end
            nothing
        end
        rm = MersenneTwister(7)
        nfinite = 0; nmatch = 0; ncombo = 0; worst = 0.0
        for combo in 1:6
            pc = rand(rm, (30, 50)); ph = rand(rm, (30, 50))
            th = (β_cc = 1 + 3rand(rm), β_ch = 0.3 + 2rand(rm), β_hc = 0.3 + 2rand(rm), β_hh = 1 + 3rand(rm),
                  γ_c = 0.5 + 1.5rand(rm), γ_h = 0.5 + 1.5rand(rm), χ_c = 0.2 + 0.4rand(rm), χ_h = 0.2 + 0.4rand(rm),
                  B_c = iseven(combo) ? 0.0 : 0.5, B_h = iseven(combo) ? 0.0 : 0.3,
                  N_c = Float64(pc), N_h = Float64(ph))
            x0 = (S_c = pc - 1, I_c = 1, S_h = ph, I_h = 0)
            g = build(th, x0, rand(rm, 3:8), 3.0 + 6rand(rm), hash((:singular, combo)))
            g === nothing && continue
            ncombo += 1
            kw = (Beta_cc = th.β_cc, Beta_ch = th.β_ch, Beta_hc = th.β_hc, Beta_hh = th.β_hh,
                  gamma_c = th.γ_c, gamma_h = th.γ_h, chi_c = th.χ_c, chi_h = th.χ_h,
                  Bc = th.B_c, Bh = th.B_h, Sc0 = x0.S_c / pc, Sh0 = x0.S_h / ph,
                  Ic0 = x0.I_c / pc, Ih0 = x0.I_h / ph, Nc = pc, Nh = ph)
            p1 = PhyloPOMP.mers_compiled_filter_pomp(g; kw...)
            p2 = PhyloPOMP.mers_compiled_filter_pomp(g; kw..., generic_singular = true)
            for s in 1:4
                seed!(s); a = logLik(pfilter(p1, Np = 100))
                seed!(s); b = logLik(pfilter(p2, Np = 100))
                if isfinite(a) || isfinite(b)
                    nfinite += 1
                    ok = isapprox(a, b; atol = 1e-9, rtol = 1e-10)
                    nmatch += ok
                    worst = max(worst, abs(a - b))
                end
            end
        end
        @info "generic vs naive singular part: $ncombo genealogies, $nfinite finite comparisons, $nmatch match, worst |Δ| = $worst"
        @test ncombo ≥ 4
        @test nfinite ≥ 10
        @test nmatch == nfinite
    end
end

@testset verbose=true "generic singular_update! (SEIR)" begin
    Infec, Expos = NaiveSEIR.Infec, NaiveSEIR.Expos
    SE = PhyloPOMP.SEIR
    θs = (β = 3.0, σ = 1.0, γ = 1.0, ω = 0.5, ψ = 0.2, χ = 0.1, N = 100.0)
    x = (S = 30, E = 4, I = 6, R = 5)

    @testset "per-outcome Δll against the closed forms" begin
        geneal = Dict(
            1 => MockNode(1, Node, 1, [2, 3], missing),
            2 => MockNode(2, Sample, 2, Int[], missing),
            3 => MockNode(3, Sample, 3, Int[], missing),
            4 => MockNode(4, Sample, 1, [2], missing),
            5 => MockNode(5, Sample, 1, Int[], missing),
        )
        mk(extra) = (c = Coloring(NaiveSEIR.Demes); plant!(c, Infec, 1);
                     foreach(e -> plant!(c, Infec, 1000 + e), 1:extra); c)
        seen = Set{Any}()
        for extra in 0:2, s in 1:200
            seed!(s); cols = mk(extra)
            Δ, x1 = singular_update!(cols, geneal, 1, x, θs, SE)
            @test Δ ≈ log(θs.β * x.S * x.I / θs.N) - log((x.E + 1) * x.I) - log(0.5) atol = 1e-12
            @test x1 == (S = 29, E = 5, I = 6, R = 5)
            push!(seen, 2 ∈ cols[Expos])

            seed!(s); cols = mk(extra)
            Δ, x1 = singular_update!(cols, geneal, 5, x, θs, SE)
            if x1.I == x.I
                @test Δ ≈ log(θs.ψ * (x.I - extra)) - log(θs.ψ / (θs.ψ + θs.χ)) atol = 1e-12
                push!(seen, :nondestr)
            else
                @test Δ ≈ log(θs.χ * x.I) - log(θs.χ / (θs.ψ + θs.χ)) atol = 1e-12
                @test x1 == (S = 30, E = 4, I = 5, R = 5)
                push!(seen, :destr)
            end

            seed!(s); cols = mk(extra)
            Δ, x1 = singular_update!(cols, geneal, 4, x, θs, SE)
            @test Δ ≈ log(θs.ψ) atol = 1e-12
            @test x1 == x && 2 ∈ cols[Infec] && 1 ∉ cols[Infec]
        end
        @test seen == Set{Any}([true, false, :nondestr, :destr])
        cols = Coloring(NaiveSEIR.Demes); plant!(cols, Expos, 1)
        @test singular_update!(cols, geneal, 5, x, θs, SE)[1] == -Inf     # sample lineage in E
        @test singular_update!(cols, geneal, 1, x, θs, SE)[1] == -Inf     # fork of an E lineage
    end

    @testset "seed-matched log-likelihood, generic vs NaiveSEIR singular part" begin
        function build(θ, x0, target_n, tmax, gseed)
            rng = MersenneTwister(gseed)
            for _ in 1:800
                g = simulate(SE, θ; x0 = x0, graft = [0, 1], tmax = tmax, rng = rng)
                nsample(g) == target_n && return g
            end
            nothing
        end
        rm = MersenneTwister(5)
        nfinite = 0; nmatch = 0; ncombo = 0; nanc = 0; worst = 0.0
        for combo in 1:8
            β = 1 + 6rand(rm); σ = 0.3 + 3rand(rm); γ = 0.3 + 3rand(rm); ω = 0.1 + 2rand(rm); ψ = 0.05 + 0.4rand(rm)
            χ = combo % 3 == 0 ? 0.0 : 0.05 + 0.3rand(rm)
            pop = rand(rm, (60, 100)); x0 = (S = pop - 1, E = 0, I = 1, R = 0)
            th = (β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, N = Float64(pop))
            g = build(th, x0, rand(rm, 2:7), 2.0 + 4rand(rm), hash((:seir, combo)))
            g === nothing && continue
            ncombo += 1
            nanc += count(i -> g[i].type == Sample && length(g[i].children) == 1, eachindex(g))
            kw = (β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, pop = pop,
                  S0 = x0.S / pop, E0 = 0.0, I0 = 1 / pop, R0 = 0.0)
            p1 = PhyloPOMP.compiled_filter_pomp(g; kw...)
            p2 = PhyloPOMP.compiled_filter_pomp(g; kw..., generic_singular = true)
            for s in 1:4
                seed!(s); a = logLik(pfilter(p1, Np = 100))
                seed!(s); b = logLik(pfilter(p2, Np = 100))
                if isfinite(a) || isfinite(b)
                    nfinite += 1
                    nmatch += isapprox(a, b; atol = 1e-8, rtol = 1e-10)
                    worst = max(worst, abs(a - b))
                end
            end
        end
        @info "SEIR generic vs naive singular part: $ncombo genealogies ($nanc sampled ancestors), $nfinite finite comparisons, $nmatch match, worst |Δ| = $worst"
        @test ncombo ≥ 5
        @test nanc ≥ 1
        @test nfinite ≥ 20
        @test nmatch == nfinite
    end
end

end # module KliSingularUpdateTest
