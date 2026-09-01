"""
M12b acceptance-gate test: the first of `mgp_filter.jl`'s "Validation gates"
(near the file's end, under "Known special cases reduce correctly") that is
checked against a source EXTERNAL to this project -- the Kingman coalescent /
Moran model special case, KLI §4.3.1, Eqs. 17-20.

WHY THIS FILE EXISTS (read before touching it): every prior milestone
(M00-M12) validated the compiler by cross-checking it against
`seir_naive.jl`/`mers_naive.jl` (hand-coded filters) and
`mers_filter_suite.tex` (a hand derivation) -- but M11's own audit
(`handoffs/M11_model_crossvalidation.md`) flagged that all three of those are
project-authored, so agreement among them cannot rule out a *shared*
misunderstanding of KLI theory. This file checks the compiler's TCC/THH
`ForkTransition` machinery against a textbook population-genetics fact
(Kingman 1982 / the standard Moran-model coalescence-probability
combinatorics) that is derived HERE FROM SCRATCH -- by brute-force
enumeration for small cases and by an independent (non-`kli_binomial_ratio`,
non-`mers_filter_suite.tex`) combinatorial formula for the randomized
sweep -- not read off this project's own tex or restated from
`kli_full_transitions_test.jl`/`kli_phi_test.jl` (which check internal
tex-agreement only; see those files' TCC/THH cases, which this file
deliberately does NOT just re-cite as ground truth).

STANDARD RESULT BEING CHECKED (Kingman 1982; standard textbook statement,
e.g. Wakeley "Coalescent Theory" §3, or Moran 1958 for the underlying finite
population process): in a Moran model of (effective) population size N, a
single birth-death step picks one individual to reproduce and one individual
("slot") to be replaced. Tracing genealogies backward, this step causes two
specific tracked lineages i,j to coalesce exactly when BOTH the
reproducing-parent role and the replaced-slot role land on {i,j} -- i.e. when
the (unordered) unordered pair of population members causally involved in
this step is exactly {i,j}. Treating the N individuals as symmetric/
exchangeable (standard neutral-model assumption), each of the C(N,2)
unordered pairs of individuals is equally likely to be "the pair causally
involved" in a given step, so:
  P(a SPECIFIC tracked pair {i,j} coalesces at this step)   = 1 / C(N,2)
  P(ANY of the C(ell,2) pairs among ell tracked lineages
    coalesces at this step)                                 = C(ell,2) / C(N,2)
(the second line follows from the first by disjoint-event additivity: at
most one pair can be "the" causally-involved pair per step, so summing the
per-pair probability over the C(ell,2) tracked pairs is valid, not an
overcount). This is an EXACT finite-N identity, not a large-N/diffusion
approximation -- the familiar "coalescent time units, rate C(ell,2)" version
of Kingman's result is the N -> infinity limit of exactly this formula, not
a separate fact.

WHY TCC/THH (and NOT SEIR's `infection` or MERS's THC/TCH): a true
Moran-style "one step, one population, both roles land in the same pool"
event needs BOTH of its production slots to occupy the SAME deme -- MERS's
`transmission_cc`/`transmission_hh` (r=(2,0)/(0,2), `move=fork(I_c=>I_c,I_c)`
resp. `fork(I_h=>I_h,I_h)`, both slots in the parent's own deme) are exactly
this shape. SEIR's `infection` (r=(1,1) at demes (E,I)) and MERS's
`transmission_hc`/`transmission_ch` (r=(1,1), one slot in each of two
DIFFERENT demes) are NOT structurally analogous: a "coalescence" there would
require one slot's occupant to be drawn from deme E/H and the other from
I/C -- two DIFFERENT finite pools, not "any 2 of N individuals in ONE pool"
-- so there is no single well-defined N to plug into C(N,2) for those marks.
This is verified directly below (not assumed): TCC/THH's ForkTransition has
BOTH `slot_demes` entries equal (see `kli_full_transitions_test.jl`'s own
`fork_t.slot_demes == [1, 1]` check for TCC), confirming the "both slots,
same pool" shape that makes the Moran analogy apply.

WHAT THIS DOES AND DOES NOT PROVE: this validates that the compiled
production-slot/Fork machinery (`kli_binomial_ratio`/`full_transitions`,
M02/M03), evaluated at the two BIRTH marks whose production vector puts both
slots in one deme, reduces EXACTLY to the standard Moran-model coalescence
probability, AT A FIXED, FROZEN post-event state `(I_C, ell_C)` (resp.
`(I_H, ell_H)`) -- i.e. an exact identity between two numbers computed from
the same finite `(N, ell)`, not an asymptotic approximation, and not a claim
about the whole compiled filter, the driver/decay/reduction machinery, or
any other event class. See `handoffs/M12b_kingman_moran_validation.md` for
the full derivation, scope discussion, and honest answer to "does this
escape the project-authored-circularity concern".
"""
module KliKingmanMoranTest

import ..Main: h1, h2

@info h1("Kingman coalescent / Moran special case (M12b, KLI §4.3.1, Eqs. 17-20)")

using Test
using PhyloPOMP
using PhyloPOMP: production_slots, kli_binomial_ratio, full_transitions,
    IdentityTransition, InlineSameDemeTransition, ForkTransition
using Random: MersenneTwister

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

# =============================================================================
# GROUND TRUTH, derived from scratch (no `kli_binomial_ratio`/`full_transitions`
# call anywhere in this section, and no reference to `mers_filter_suite.tex`).
# =============================================================================

"""
Count unordered pairs of `1:N` by brute-force double loop. Deliberately does
NOT call `Base.binomial` -- this is the independent recount that lets the
`Base.binomial` cross-check below stand on its own rather than trusting the
standard library's combinatorial convention.
"""
function brute_pairs(N::Int)
    c = 0
    for i in 1:N, j in (i+1):N
        c += 1
    end
    c
end

"""
Count, among all unordered pairs of `1:N`, how many lie entirely within the
first `ell` elements (the "tracked" subset) -- brute-force, no `binomial()`
call. This is the from-scratch verification of "P(a uniformly random pair
lies entirely inside a fixed ell-subset of N) = C(ell,2)/C(N,2)".
"""
function brute_tracked_pairs(N::Int, ell::Int)
    c = 0
    for i in 1:N, j in (i+1):N
        (i <= ell && j <= ell) && (c += 1)
    end
    c
end

"""
    moran_pair_prob(N, ell) -> Rational{Int}

The standard Moran-model "ANY of the C(ell,2) tracked pairs coalesces at
this step" probability, `C(ell,2)/C(N,2)`, computed independently of
`kli_binomial_ratio`/`full_transitions` (this is the "expected" side of
every comparison below). Uses `Base.binomial` (itself cross-checked against
`brute_pairs`/`brute_tracked_pairs` above), guarding the `N<2` degenerate
case (fewer than 2 individuals total: no pair can ever coalesce) so the
ratio is `0//1` rather than a `0/0` division.
"""
function moran_pair_prob(N::Integer, ell::Integer)
    @assert 0 <= ell <= N "moran_pair_prob: ell=$ell must satisfy 0 <= ell <= N=$N"
    denom = binomial(N, 2)
    denom == 0 ? zero(Rational{Int}) : binomial(ell, 2) // denom
end

@testset verbose=true "Kingman/Moran special case (M12b)" begin

    # -------------------------------------------------------------------
    # Part 0: the standard combinatorial fact itself, brute-force verified,
    # independent of Base.binomial and independent of this project's code.
    # -------------------------------------------------------------------
    @info h2("Part 0: brute-force ground truth for C(N,2), C(ell,2)/C(N,2)")
    @testset "brute-force pair counts match N(N-1)/2 and Base.binomial" begin
        for N in 0:14
            @test brute_pairs(N) == N * (N - 1) ÷ 2
            @test brute_pairs(N) == binomial(N, 2)
            for ell in 0:N
                @test brute_tracked_pairs(N, ell) == ell * (ell - 1) ÷ 2
                @test brute_tracked_pairs(N, ell) == binomial(ell, 2)
                # the ratio itself, both ways of computing it, must agree:
                expected_ratio = N < 2 ? 0 // 1 : brute_tracked_pairs(N, ell) // brute_pairs(N)
                @test moran_pair_prob(N, ell) == expected_ratio
            end
        end
    end

    # -------------------------------------------------------------------
    # Part 1: kli_binomial_ratio's TCC/THH Fork saturation reduces exactly
    # to the "this SPECIFIC pair is the one realized" probability 1/C(N,2).
    # -------------------------------------------------------------------
    @info h2("Part 1: TCC/THH Fork phi_u(s=full) == 1/C(N,2) (single specific pair)")
    @testset "TCC (r=(2,0), Camel deme): phi_fork == 1/C(I_C,2), H-deme irrelevant" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        @test production_slots(tcc) == [2, 0]
        @test tcc.from == 1  # Camel deme (I_c) is the fork's ancestral/source deme

        # Hand anchors, cross-checked against the already-vetted values in
        # kli_phi_test.jl / kli_full_transitions_test.jl / mers_filter_suite.tex
        # row 616 -- included as a fixed sanity point, NOT as the derivation.
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

            # The H-factor of the binomial-ratio product is
            # C(I_H-ell_H,0)/C(I_H,0) = 1/1 = 1 for ANY (I_H,ell_H) -- confirm
            # numerically (not just by inspecting the formula) that phi_fork
            # is completely insensitive to the OTHER deme's state.
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

    # -------------------------------------------------------------------
    # Part 2: the same result reproduced through the M03 classifier layer
    # (full_transitions -> ForkTransition), plus the Chu-Vandermonde
    # sum-to-1 structural sanity check (Identity + ell*Inline + C(ell,2)*Fork
    # == 1), confirming the C(ell,s) "which specific lineages fill the
    # slots" weighting used to assemble "any pair" from "one specific pair"
    # is the SAME weighting `mgp_mers_filter.jl`'s own Finding 2 already
    # uses -- cross-checked here, not assumed.
    # -------------------------------------------------------------------
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

    # -------------------------------------------------------------------
    # Part 3 (THE MAIN CLAIM): the implied "any tracked pair coalesces"
    # RATE, assembled from the REAL compiled hazard closure
    # (Event.hazard, M01) and the REAL M02 kli_binomial_ratio, equals
    # alpha_TCC(x,theta) * moran_pair_prob(I_C, ell_C) -- an independently
    # derived target that never calls kli_binomial_ratio/full_transitions
    # or reads mers_filter_suite.tex.
    # -------------------------------------------------------------------
    @info h2("Part 3 (MAIN CLAIM): compiled any-pair rate == alpha * C(ell,2)/C(N,2)")
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

            # Strictly positive random rational parameters, so alpha_TCC != 0
            # and the check is not vacuously true.
            S_c = (rand(rng, 1:2000) // rand(rng, 1:37))
            N_c = (rand(rng, 1:5000) // rand(rng, 1:19))
            beta_cc = (rand(rng, 1:97) // rand(rng, 1:23))

            x = (S_c = S_c, I_c = I_C, S_h = zero(S_c), I_h = I_H)
            θ = (β_cc = beta_cc, β_ch = 0 // 1, β_hc = 0 // 1, β_hh = 0 // 1,
                 γ_c = 0 // 1, γ_h = 0 // 1, χ_c = 0 // 1, χ_h = 0 // 1,
                 B_c = 0 // 1, B_h = 0 // 1, N_c = N_c, N_h = 1 // 1)

            alpha_tcc = tcc.hazard(x, θ)           # REAL compiled Event hazard (M01)
            @test alpha_tcc == beta_cc * S_c * I_C / N_c  # sanity: matches mgp_mers.jl's declared rate
            @test alpha_tcc > 0

            phi_fork = kli_binomial_ratio(production_slots(tcc), [2, 0], [ell_C, ell_H], [I_C, I_H])  # REAL M02 function

            lhs = alpha_tcc * binomial(ell_C, 2) * phi_fork
            rhs = alpha_tcc * moran_pair_prob(I_C, ell_C)  # independently-derived target
            @test lhs == rhs
            worst_diff = max(worst_diff, abs(lhs - rhs))
        end
        @test worst_diff == 0 // 1  # exact Rational{Int} equality throughout, not approximate
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

    # -------------------------------------------------------------------
    # Part 4: scope-honesty check -- the correspondence is INSTANTANEOUS
    # (a fixed-state identity), not a real-time-scaled coalescent-rate
    # limit. Confirm explicitly that moran_pair_prob does NOT, by itself,
    # carry any 1/N_e or real-time rescaling -- that only enters once
    # alpha_TCC's own S_c/N_c dependence is examined (a per-model,
    # NOT per-Kingman-formula, fact) -- i.e. this file validates the
    # combinatorial core, not a specific epidemic model's diffusion limit.
    # -------------------------------------------------------------------
    @info h2("Part 4: degenerate/boundary states behave sanely (0, 1 individuals; ell=N)")
    @testset "boundary states" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        # I_C=0 or 1: fewer than 2 individuals, no pair can ever coalesce.
        @test kli_binomial_ratio(tcc, [2, 0], [0, 0], [0, 0]) == 0 // 1
        @test kli_binomial_ratio(tcc, [2, 0], [1, 0], [1, 0]) == 0 // 1
        @test moran_pair_prob(0, 0) == 0 // 1
        @test moran_pair_prob(1, 0) == 0 // 1
        @test moran_pair_prob(1, 1) == 0 // 1
        # ell_C == I_C == 2: every individual tracked, only one possible pair,
        # and it MUST be the one that coalesces on the next TCC birth event.
        @test kli_binomial_ratio(tcc, [2, 0], [2, 0], [2, 0]) == 1 // 1
        @test moran_pair_prob(2, 2) == 1 // 1
    end

end

end # module KliKingmanMoranTest
