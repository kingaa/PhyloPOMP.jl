"""
Reference-value tests for the generic KLI production-slot math added in
M02 (`src/examples/mgp_phi.jl`): `production_slots`, `enumerate_saturations`,
`kli_binomial_ratio` (φ_u).

Every numeric case here is hand-computed in the M02 task spec / handoff
(`handoffs/M02_kli_phi.md`) using exact rational arithmetic, then asserted
against the generic function -- not just checked for "looks reasonable".
Several cases additionally cross-check the generic output against real
`Event`s pulled from the actual `MERS`/`SEIR` `MGPModel`s (not hand-typed
literals) so that no saturation list is hard-coded per model, per the M02
acceptance gate. The MERS cross-checks are verified against this repo's
own worked derivation, `src/examples/mers_filter_suite.tex`:
  - Step B (enumerate saturations): lines 443-449
  - Step C (binomial ratio formula): lines 450-457
  - TCC (within-camel birth, r=(2,0)):  lines 472-480, table rows 554-556
  - THH (within-human birth, r=(0,2)):  lines 482-485, table rows 557-559
  - THC (camel-to-human spillover, r=(1,1)): lines 487-498, table rows 560-563
  - TCH (human-to-camel spillover, r=(1,1)): lines 500-509, table rows 564-567
  - RC/RH/SC/SH (r=(0,0)):  lines 511-539, table rows 568-571
"""
module KliPhiTest

import ..Main: h1, h2

@info h1("KLI φ_u: production slots, saturation enumeration, binomial ratio")

using Test
using PhyloPOMP
using PhyloPOMP: production_slots, enumerate_saturations, kli_binomial_ratio

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

@testset verbose=true "KLI φ_u (M02)" begin

    @info h2("r = (0): trivial saturation, φ = Q always")
    @testset "r=(0)" begin
        r, ℓ, n = [0], [3], [7]
        S = enumerate_saturations(r, ℓ)
        @test S == [[0]]
        @test kli_binomial_ratio(r, [0], ℓ, n) == 1 // 1
        @test kli_binomial_ratio(r, [0], ℓ, n; Q = 0) == 0 // 1
        @test kli_binomial_ratio(r, [0], ℓ, n; Q = 1 // 3) == 1 // 3

        # cross-check against a real r=(0,0) event: SEIR recovery (DEATH).
        recovery = find_event(PhyloPOMP.SEIR, :recovery)
        @test production_slots(recovery) == [0, 0]
        @test enumerate_saturations(recovery, [1, 2]) == [[0, 0]]
        @test kli_binomial_ratio(recovery, [0, 0], [1, 2], [3, 5]) == 1 // 1
    end

    @info h2("r = (1): single deme, single production")
    @testset "r=(1)" begin
        # n=5, ℓ=2: φ(0) = C(3,1)/C(5,1) = 3/5, φ(1) = C(3,0)/C(5,1) = 1/5.
        r, ℓ, n = [1], [2], [5]
        S = enumerate_saturations(r, ℓ)
        @test S == [[0], [1]]
        @test kli_binomial_ratio(r, [0], ℓ, n) == 3 // 5
        @test kli_binomial_ratio(r, [1], ℓ, n) == 1 // 5
        # Chu-Vandermonde sanity: sum_s C(ℓ,s)*φ(s) == 1.
        @test binomial(2, 0) * (3 // 5) + binomial(2, 1) * (1 // 5) == 1

        # ℓ=0: saturation collapses to {0} only.
        @test enumerate_saturations([1], [0]) == [[0]]
        @test kli_binomial_ratio([1], [0], [0], [5]) == 1 // 1  # C(5,1)/C(5,1)

        # cross-check against SEIR progression (MIGRATION, r=(0,1) -- the
        # I-deme slot is the "r=1" component; E-deme slot is a fixed r=0).
        progression = find_event(PhyloPOMP.SEIR, :progression)
        @test production_slots(progression) == [0, 1]
        Sp = enumerate_saturations(progression, [3, 2])  # ℓ_E=3 (irrelevant, r_E=0), ℓ_I=2
        @test Sp == [[0, 0], [0, 1]]
        @test kli_binomial_ratio(progression, [0, 0], [3, 2], [9, 5]) == 3 // 5
        @test kli_binomial_ratio(progression, [0, 1], [3, 2], [9, 5]) == 1 // 5
    end

    @info h2("r = (2): single deme, two productions")
    @testset "r=(2)" begin
        # n=5, ℓ=2: φ(0)=C(3,2)/C(5,2)=3/10, φ(1)=C(3,1)/C(5,2)=3/10,
        # φ(2)=C(3,0)/C(5,2)=1/10.
        r, ℓ, n = [2], [2], [5]
        S = enumerate_saturations(r, ℓ)
        @test S == [[0], [1], [2]]
        @test kli_binomial_ratio(r, [0], ℓ, n) == 3 // 10
        @test kli_binomial_ratio(r, [1], ℓ, n) == 3 // 10
        @test kli_binomial_ratio(r, [2], ℓ, n) == 1 // 10
        @test binomial(2, 0) * (3 // 10) + binomial(2, 1) * (3 // 10) +
              binomial(2, 2) * (1 // 10) == 1

        # cross-check against MERS transmission_cc (TCC, r=(2,0)), I_C=5,
        # ℓ_C=2, H-deme irrelevant (r_H=0). mers_filter_suite.tex:472-480,
        # 554-556: s_C=0 -> (I_C-ℓ_C)(I_C-ℓ_C-1)/(I_C(I_C-1)) = 3*2/(5*4) =
        # 3/10; s_C=1 -> 2(I_C-ℓ_C)/(I_C(I_C-1)) = 6/20 = 3/10;
        # s_C=2 -> 2/(I_C(I_C-1)) = 2/20 = 1/10.
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        @test production_slots(tcc) == [2, 0]
        Stcc = enumerate_saturations(tcc, [2, 0])
        @test Stcc == [[0, 0], [1, 0], [2, 0]]
        @test kli_binomial_ratio(tcc, [0, 0], [2, 0], [5, 3]) == 3 // 10
        @test kli_binomial_ratio(tcc, [1, 0], [2, 0], [5, 3]) == 3 // 10
        @test kli_binomial_ratio(tcc, [2, 0], [2, 0], [5, 3]) == 1 // 10

        # mirror: transmission_hh (THH, r=(0,2)), I_H=4, ℓ_H=1.
        # mers_filter_suite.tex:482-485, 557-559: s_H=0 -> (3*2)/(4*3)=1/2;
        # s_H=1 -> 2*3/(4*3)=1/2; s_H=2 -> 2/(4*3)=1/6.
        thh = find_event(PhyloPOMP.MERS, :transmission_hh)
        @test production_slots(thh) == [0, 2]
        @test kli_binomial_ratio(thh, [0, 0], [0, 1], [3, 4]) == 1 // 2
        @test kli_binomial_ratio(thh, [0, 1], [0, 1], [3, 4]) == 1 // 2
        @test kli_binomial_ratio(thh, [0, 2], [0, 1], [3, 4]) == 1 // 6
    end

    @info h2("r = (1,1): two demes, one production each (cross-deme fork)")
    @testset "r=(1,1)" begin
        # n=(5,4), ℓ=(2,1): the 2x2 product saturation space.
        r, ℓ, n = [1, 1], [2, 1], [5, 4]
        S = enumerate_saturations(r, ℓ)
        @test Set(S) == Set([[0, 0], [1, 0], [0, 1], [1, 1]])
        @test kli_binomial_ratio(r, [0, 0], ℓ, n) == 9 // 20
        @test kli_binomial_ratio(r, [1, 0], ℓ, n) == 3 // 20
        @test kli_binomial_ratio(r, [0, 1], ℓ, n) == 3 // 20
        @test kli_binomial_ratio(r, [1, 1], ℓ, n) == 1 // 20

        # cross-check against MERS transmission_hc (THC, camel-to-human
        # spillover). mers_filter_suite.tex:487-498, 560-563:
        # (s_C,s_H)=(0,0) -> (I_C-ℓ_C)(I_H-ℓ_H)/(I_CI_H) = 3*3/20 = 9/20;
        # (0,1) -> (I_C-ℓ_C)/(I_CI_H) = 3/20; (1,0) -> (I_H-ℓ_H)/(I_CI_H) =
        # 3/20; (1,1) -> 1/(I_CI_H) = 1/20.
        thc = find_event(PhyloPOMP.MERS, :transmission_hc)
        @test production_slots(thc) == [1, 1]
        Sthc = enumerate_saturations(thc, [2, 1])
        @test Set(Sthc) == Set([[0, 0], [1, 0], [0, 1], [1, 1]])
        @test kli_binomial_ratio(thc, [0, 0], [2, 1], [5, 4]) == 9 // 20
        @test kli_binomial_ratio(thc, [0, 1], [2, 1], [5, 4]) == 3 // 20
        @test kli_binomial_ratio(thc, [1, 0], [2, 1], [5, 4]) == 3 // 20
        @test kli_binomial_ratio(thc, [1, 1], [2, 1], [5, 4]) == 1 // 20

        # mirror: transmission_ch (TCH, human-to-camel spillover), same
        # (n,ℓ). mers_filter_suite.tex:500-509, 564-567 -- algebraically
        # identical product formula for r=(1,1), independently re-derived
        # in the tex with C/H roles swapped; the generic function must
        # reproduce the same numeric values since it only consumes r,s,ℓ,n.
        tch = find_event(PhyloPOMP.MERS, :transmission_ch)
        @test production_slots(tch) == [1, 1]
        @test kli_binomial_ratio(tch, [0, 0], [2, 1], [5, 4]) == 9 // 20
        @test kli_binomial_ratio(tch, [1, 0], [2, 1], [5, 4]) == 3 // 20
        @test kli_binomial_ratio(tch, [0, 1], [2, 1], [5, 4]) == 3 // 20
        @test kli_binomial_ratio(tch, [1, 1], [2, 1], [5, 4]) == 1 // 20
    end

    @info h2("boundary: ℓ_d = 0 collapses that deme's saturation range")
    @testset "boundary ℓ=0" begin
        # one deme zero, one deme nonzero.
        S = enumerate_saturations([1, 1], [0, 3])
        @test Set(S) == Set([[0, 0], [0, 1]])
        @test all(s[1] == 0 for s in S)

        # all demes zero: saturation forced to all-zero, φ = Q (worked by
        # hand: s_d must be 0, so C(n_d-0, r_d-0)/C(n_d, r_d) =
        # C(n_d,r_d)/C(n_d,r_d) = 1 for every d, so φ_u = Q * 1 = Q).
        r, ℓ, n = [1, 1], [0, 0], [5, 4]
        S0 = enumerate_saturations(r, ℓ)
        @test S0 == [[0, 0]]
        @test kli_binomial_ratio(r, [0, 0], ℓ, n) == 1 // 1
        @test kli_binomial_ratio(r, [0, 0], ℓ, n; Q = 0) == 0 // 1
    end

    @info h2("boundary: ℓ_d = n_d (fully saturated deme)")
    @testset "boundary ℓ=n" begin
        # single deme, r=1, ℓ=n=5: s=0 forced-impossible (φ=0), s=1 forced
        # (C(0,0)/C(5,1) = 1/5) -- by hand: C(n-ℓ, r-s) = C(0, 1-s), which
        # is 1 iff s == r (=1) and 0 otherwise.
        r, ℓ, n = [1], [5], [5]
        S = enumerate_saturations(r, ℓ)
        @test S == [[0], [1]]
        @test kli_binomial_ratio(r, [0], ℓ, n) == 0 // 1
        @test kli_binomial_ratio(r, [1], ℓ, n) == 1 // 5
    end

    @info h2("boundary: n_d < r_d (deme cannot produce r_d individuals)")
    @testset "boundary n<r" begin
        # r=3 requested from a deme with only n=2 total individuals:
        # C(n,r) = C(2,3) = 0 (b>a) -> denominator zero -> this
        # implementation's documented defensive convention is to return
        # exactly 0 (not throw, not NaN/Inf) for every saturation.
        r, ℓ, n = [3], [1], [2]
        S = enumerate_saturations(r, ℓ)  # min(r,ℓ)=1, so {0,1}; enumeration
        @test S == [[0], [1]]            # itself is still well-defined.
        @test kli_binomial_ratio(r, [0], ℓ, n) == 0 // 1
        @test kli_binomial_ratio(r, [1], ℓ, n) == 0 // 1
    end

    @info h2("enumerate_saturations is always correctly bounded")
    @testset "enumeration bounds (no impossible saturations generated)" begin
        cases = [
            ([0], [0]), ([0], [5]), ([1], [0]), ([1], [1]), ([1], [10]),
            ([2], [1]), ([2], [2]), ([3], [1]),
            ([1, 1], [0, 0]), ([1, 1], [1, 0]), ([1, 1], [0, 1]),
            ([2, 3], [1, 2]), ([0, 2], [4, 1]),
        ]
        for (r, ℓ) in cases
            S = enumerate_saturations(r, ℓ)
            @test !isempty(S)
            for s in S
                @test length(s) == length(r)
                for d in eachindex(r)
                    @test 0 <= s[d] <= min(r[d], ℓ[d])
                end
            end
            # every combinatorially-expected saturation is present exactly once
            expected = length(S) == prod(min(r[d], ℓ[d]) + 1 for d in eachindex(r))
            @test expected
            @test allunique(S)
        end
    end

    @info h2("kli_binomial_ratio is defensive against an out-of-bounds s")
    @testset "defensive φ_u on s outside [0, min(r,ℓ)]" begin
        # s=2 > r=1: numerator argument r-s = -1 < 0 -> safe_binomial -> 0.
        @test kli_binomial_ratio([1], [2], [5], [10]) == 0 // 1
        # s=-1 < 0: r-s = 2 is NOT negative (it's *larger* than r), so the
        # b<0 rule alone would not catch it -- verified this is a real trap
        # (safe_binomial(5,2)/safe_binomial(10,1) = 10/10 = 1 != 0) before
        # adding an explicit `any(s .< 0)` guard to kli_binomial_ratio; now
        # it correctly returns 0 for any negative-count saturation.
        @test kli_binomial_ratio([1], [-1], [5], [10]) == 0 // 1
    end

    @info h2("r=(0,0) events (RC/RH/SC/SH): φ=1 always, no fork")
    @testset "RC/RH/SC/SH (r=(0,0))" begin
        # mers_filter_suite.tex:511-539, table rows 568-571: r=(0,0),
        # s=(0,0) forced, φ=1 (the empty product), for removal and
        # sample/death alike (the Supp indicator is a separate, not-yet-
        # implemented compatibility condition -- out of scope for M02).
        for name in (:removal_c, :removal_h, :sampling_c, :sampling_h)
            ev = find_event(PhyloPOMP.MERS, name)
            @test production_slots(ev) == [0, 0]
            @test enumerate_saturations(ev, [3, 2]) == [[0, 0]]
            @test kli_binomial_ratio(ev, [0, 0], [3, 2], [9, 6]) == 1 // 1
        end
    end

    @info h2("mismatched-length arguments throw ArgumentError")
    @testset "argument validation" begin
        @test_throws ArgumentError enumerate_saturations([1, 1], [1])
        @test_throws ArgumentError kli_binomial_ratio([1, 1], [0, 0], [1], [5, 5])
        @test_throws ArgumentError kli_binomial_ratio([1, 1], [0, 0], [1, 1], [5])
    end

end

end
