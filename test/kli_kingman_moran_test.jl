"""
Kingman/Moran special case (KLI 4.3.1, Eqs. 17-20).

TCC and THH put both fork slots in one deme, so the fork weight is 1/C(N,2)
and the any-pair rate is alpha*C(ell,2)/C(N,2). Expected values come from
brute-force pair counts, not the compiler.
"""
module KliKingmanMoranTest

import ..Main: h1, h2

@info h1("Kingman coalescent / Moran special case (KLI 4.3.1)")

using Test
using PhyloPOMP
using PhyloPOMP: production_slots, kli_binomial_ratio, full_transitions,
    IdentityTransition, InlineSameDemeTransition, ForkTransition
using Random: MersenneTwister

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

# Expected values, computed without the compiler.

"""Count unordered pairs of `1:N` by double loop."""
function brute_pairs(N::Int)
    c = 0
    for i in 1:N, j in (i+1):N
        c += 1
    end
    c
end

"""Count pairs of `1:N` lying inside the first `ell` elements."""
function brute_tracked_pairs(N::Int, ell::Int)
    c = 0
    for i in 1:N, j in (i+1):N
        (i <= ell && j <= ell) && (c += 1)
    end
    c
end

"""C(ell,2)/C(N,2) as an exact Rational; 0 when N < 2."""
function moran_pair_prob(N::Integer, ell::Integer)
    @assert 0 <= ell <= N "moran_pair_prob: ell=$ell must satisfy 0 <= ell <= N=$N"
    denom = binomial(N, 2)
    denom == 0 ? zero(Rational{Int}) : binomial(ell, 2) // denom
end

@testset verbose=true "Kingman/Moran special case" begin

    @info h2("Part 0: brute-force ground truth for C(N,2), C(ell,2)/C(N,2)")
    @testset "brute-force pair counts match N(N-1)/2 and Base.binomial" begin
        for N in 0:14
            @test brute_pairs(N) == N * (N - 1) ÷ 2
            @test brute_pairs(N) == binomial(N, 2)
            for ell in 0:N
                @test brute_tracked_pairs(N, ell) == ell * (ell - 1) ÷ 2
                @test brute_tracked_pairs(N, ell) == binomial(ell, 2)
                expected_ratio = N < 2 ? 0 // 1 : brute_tracked_pairs(N, ell) // brute_pairs(N)
                @test moran_pair_prob(N, ell) == expected_ratio
            end
        end
    end

    @info h2("Part 1: TCC/THH Fork phi_u(s=full) == 1/C(N,2) (single specific pair)")
    @testset "TCC (r=(2,0), Camel deme): phi_fork == 1/C(I_C,2), H-deme irrelevant" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        @test production_slots(tcc) == [2, 0]
        @test tcc.from == 1  # Camel deme (I_c) is the fork's ancestral/source deme

        # Anchor values: 1/C(5,2), 1/C(2,2), 1/C(100,2).
        @test kli_binomial_ratio(tcc, [2, 0], [2, 0], [5, 3]) == 1 // 10   # 1/C(5,2)
        @test kli_binomial_ratio(tcc, [2, 0], [2, 0], [2, 0]) == 1 // 1   # 1/C(2,2): guaranteed coalescence
        @test kli_binomial_ratio(tcc, [2, 0], [0, 0], [100, 0]) == 1 // 4950  # 1/C(100,2)

        rng = MersenneTwister(20260818)
        ntrials = 200
        for _ in 1:ntrials
            I_C = rand(rng, 0:400)
            ell_C = rand(rng, 0:I_C)
            I_H = rand(rng, 0:400)
            ell_H = rand(rng, 0:I_H)

            phi_fork = kli_binomial_ratio(production_slots(tcc), [2, 0], [ell_C, ell_H], [I_C, I_H])
            expected = I_C >= 2 ? 1 // binomial(I_C, 2) : zero(Rational{Int})
            @test phi_fork == expected

            # The H-deme factor is 1, so phi_fork must not depend on the other deme.
            I_H2 = rand(rng, 0:400)
            ell_H2 = rand(rng, 0:I_H2)
            phi_fork2 = kli_binomial_ratio(production_slots(tcc), [2, 0], [ell_C, ell_H2], [I_C, I_H2])
            @test phi_fork2 == phi_fork
        end
    end

    @testset "THH (r=(0,2), Human deme): phi_fork == 1/C(I_H,2), mirror of TCC" begin
        thh = find_event(PhyloPOMP.MERS, :transmission_hh)
        @test production_slots(thh) == [0, 2]
        @test thh.from == 2  # Human deme (I_h)

        @test kli_binomial_ratio(thh, [0, 2], [0, 1], [3, 4]) == 1 // 6   # 1/C(4,2)

        rng = MersenneTwister(20260819)
        ntrials = 200
        for _ in 1:ntrials
            I_H = rand(rng, 0:400)
            ell_H = rand(rng, 0:I_H)
            I_C = rand(rng, 0:400)
            ell_C = rand(rng, 0:I_C)

            phi_fork = kli_binomial_ratio(production_slots(thh), [0, 2], [ell_C, ell_H], [I_C, I_H])
            expected = I_H >= 2 ? 1 // binomial(I_H, 2) : zero(Rational{Int})
            @test phi_fork == expected

            I_C2 = rand(rng, 0:400)
            ell_C2 = rand(rng, 0:I_C2)
            phi_fork2 = kli_binomial_ratio(production_slots(thh), [0, 2], [ell_C2, ell_H], [I_C2, I_H])
            @test phi_fork2 == phi_fork
        end
    end

    # Same weight via full_transitions; Identity + ell*Inline + C(ell,2)*Fork sums to 1.
    @info h2("Part 2: full_transitions' ForkTransition + Chu-Vandermonde cross-check")
    @testset "TCC via full_transitions" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        rng = MersenneTwister(555)
        ntrials = 100
        for _ in 1:ntrials
            I_C = rand(rng, 2:200)
            ell_C = rand(rng, 2:I_C)   # >= 2 so the Fork saturation is enumerable
            I_H = rand(rng, 0:50)
            ell_H = rand(rng, 0:I_H)

            ts = full_transitions(tcc, [ell_C, ell_H], [I_C, I_H])
            fork_t = only(filter(t -> t isa ForkTransition, ts))
            @test fork_t.ancestral_deme == 1
            @test fork_t.slot_demes == [1, 1]   # BOTH slots land in the same (Camel) deme
            @test fork_t.phi == 1 // binomial(I_C, 2)

            id_t = only(filter(t -> t isa IdentityTransition, ts))
            inl_t = only(filter(t -> t isa InlineSameDemeTransition, ts))
            @test binomial(ell_C, 0) * id_t.phi + binomial(ell_C, 1) * inl_t.phi +
                  binomial(ell_C, 2) * fork_t.phi == 1
        end
    end

    # Any-pair rate from the compiled hazard and fork weight equals alpha*C(ell,2)/C(N,2).
    @info h2("Part 3: compiled any-pair rate == alpha * C(ell,2)/C(N,2)")
    @testset "TCC: alpha_TCC * C(ell_C,2) * phi_fork == alpha_TCC * moran_pair_prob(I_C,ell_C)" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        rng = MersenneTwister(20260820)
        ntrials = 150
        worst_diff = 0 // 1
        for _ in 1:ntrials
            I_C = rand(rng, 2:500)
            ell_C = rand(rng, 0:I_C)
            I_H = rand(rng, 0:200)
            ell_H = rand(rng, 0:I_H)

            # Positive random parameters keep alpha nonzero.
            S_c = (rand(rng, 1:2000) // rand(rng, 1:37))
            N_c = (rand(rng, 1:5000) // rand(rng, 1:19))
            beta_cc = (rand(rng, 1:97) // rand(rng, 1:23))

            x = (S_c = S_c, I_c = I_C, S_h = zero(S_c), I_h = I_H)
            θ = (β_cc = beta_cc, β_ch = 0 // 1, β_hc = 0 // 1, β_hh = 0 // 1,
                 γ_c = 0 // 1, γ_h = 0 // 1, χ_c = 0 // 1, χ_h = 0 // 1,
                 B_c = 0 // 1, B_h = 0 // 1, N_c = N_c, N_h = 1 // 1)

            alpha_tcc = tcc.hazard(x, θ)
            @test alpha_tcc == beta_cc * S_c * I_C / N_c  # alpha_TCC matches the declared MERS rate
            @test alpha_tcc > 0

            phi_fork = kli_binomial_ratio(production_slots(tcc), [2, 0], [ell_C, ell_H], [I_C, I_H])

            lhs = alpha_tcc * binomial(ell_C, 2) * phi_fork
            rhs = alpha_tcc * moran_pair_prob(I_C, ell_C)
            @test lhs == rhs
            worst_diff = max(worst_diff, abs(lhs - rhs))
        end
        @test worst_diff == 0 // 1
    end

    @testset "THH: alpha_THH * C(ell_H,2) * phi_fork == alpha_THH * moran_pair_prob(I_H,ell_H)" begin
        thh = find_event(PhyloPOMP.MERS, :transmission_hh)
        rng = MersenneTwister(20260821)
        ntrials = 150
        worst_diff = 0 // 1
        for _ in 1:ntrials
            I_H = rand(rng, 2:500)
            ell_H = rand(rng, 0:I_H)
            I_C = rand(rng, 0:200)
            ell_C = rand(rng, 0:I_C)

            S_h = (rand(rng, 1:2000) // rand(rng, 1:37))
            N_h = (rand(rng, 1:5000) // rand(rng, 1:19))
            beta_hh = (rand(rng, 1:97) // rand(rng, 1:23))

            x = (S_c = zero(S_h), I_c = I_C, S_h = S_h, I_h = I_H)
            θ = (β_cc = 0 // 1, β_ch = 0 // 1, β_hc = 0 // 1, β_hh = beta_hh,
                 γ_c = 0 // 1, γ_h = 0 // 1, χ_c = 0 // 1, χ_h = 0 // 1,
                 B_c = 0 // 1, B_h = 0 // 1, N_c = 1 // 1, N_h = N_h)

            alpha_thh = thh.hazard(x, θ)
            @test alpha_thh == beta_hh * S_h * I_H / N_h
            @test alpha_thh > 0

            phi_fork = kli_binomial_ratio(production_slots(thh), [0, 2], [ell_C, ell_H], [I_C, I_H])

            lhs = alpha_thh * binomial(ell_H, 2) * phi_fork
            rhs = alpha_thh * moran_pair_prob(I_H, ell_H)
            @test lhs == rhs
            worst_diff = max(worst_diff, abs(lhs - rhs))
        end
        @test worst_diff == 0 // 1
    end

    # Boundary states: fewer than 2 individuals gives 0; ell = N = 2 gives 1.
    @info h2("Part 4: degenerate/boundary states behave sanely (0, 1 individuals; ell=N)")
    @testset "boundary states" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        # I_C=0 or 1: fewer than 2 individuals, no pair can ever coalesce.
        @test kli_binomial_ratio(tcc, [2, 0], [0, 0], [0, 0]) == 0 // 1
        @test kli_binomial_ratio(tcc, [2, 0], [1, 0], [1, 0]) == 0 // 1
        @test moran_pair_prob(0, 0) == 0 // 1
        @test moran_pair_prob(1, 0) == 0 // 1
        @test moran_pair_prob(1, 1) == 0 // 1
        # ell_C == I_C == 2: the only pair must coalesce.
        @test kli_binomial_ratio(tcc, [2, 0], [2, 0], [2, 0]) == 1 // 1
        @test moran_pair_prob(2, 2) == 1 // 1
    end

end

end # module KliKingmanMoranTest
