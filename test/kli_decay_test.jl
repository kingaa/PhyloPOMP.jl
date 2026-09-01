"""
M07 Part A acceptance-gate tests: `decay_contribution` / `total_decay`
(`src/examples/mgp_decay.jl`), the generic lambda(t,x,y) (KLI Eq. 47 /
Appendix B Eq. B2) generalized from `mers_filter_suite.tex`'s MERS-specific
closed form (the "Decay lambda" subsection, and the "RC and RH"/"SC and SH"
event-specific derivations) to an arbitrary DEATH/SAMPLE `Event`.

Every concrete instance below is hand-computed independently (not copied
from an earlier milestone's table, since this is the first milestone to
touch DEATH/SAMPLE events at all -- M03 explicitly scoped them out of
`full_transitions`), using exact `Rational{Int}` arithmetic throughout
(`x`/`θ` NamedTuples with `Rational{Int}` fields, so `event.hazard(x,θ)`
itself returns an exact `Rational{Int}`, per this milestone's "exact
arithmetic, not floating point" instruction).
"""
module KliDecayTest

import ..Main: h1, h2

@info h1("Decay lambda: decay_contribution / total_decay (M07 Part A)")

using Test
using PhyloPOMP
using PhyloPOMP: Event, EventType, BIRTH, MIGRATION, DEATH, SAMPLE, NEUTRAL,
    DecayContribution, decay_contribution, total_decay

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

@testset verbose=true "Decay lambda (M07 Part A)" begin

    @info h2("MERS: concrete (I_C,I_H,ell_C,ell_H) exercising both the " *
             "above-threshold and at-threshold indicator branches")
    @testset "MERS concrete instance" begin
        M = PhyloPOMP.MERS
        removal_c   = find_event(M, :removal_c)
        removal_h   = find_event(M, :removal_h)
        sampling_c  = find_event(M, :sampling_c)
        sampling_h  = find_event(M, :sampling_h)

        # gamma_c=1/2, gamma_h=2/5, chi_c=1/10, chi_h=3/10 -- exact Rationals.
        x = (S_c = 10 // 1, I_c = 5 // 1, S_h = 10 // 1, I_h = 3 // 1)
        θ = (β_cc = 1 // 1, β_ch = 1 // 1, β_hc = 1 // 1, β_hh = 1 // 1,
             γ_c = 1 // 2, γ_h = 2 // 5, χ_c = 1 // 10, χ_h = 3 // 10,
             B_c = 1 // 1, B_h = 1 // 1, N_c = 100 // 1, N_h = 100 // 1)

        # I_C=5, ell_C=2: strictly above threshold (5 > 2) -> RC gate = 0.
        # I_H=3, ell_H=3: exactly at threshold (3 == 3)     -> RH gate = 1.
        ℓ = [2, 3]
        n = [5, 3]

        # Hand-computed against mers_filter_suite.tex's boxed "Decay lambda"
        # formula (lines 833-836):
        #   lambda = chi_C I_C + chi_H I_H
        #          + gamma_C I_C*1{I_C<=ell_C} + gamma_H I_H*1{I_H<=ell_H}
        #          = (1/10)*5 + (3/10)*3 + (1/2)*5*0 + (2/5)*3*1
        #          = 1/2 + 9/10 + 0 + 6/5 = 13/5
        dc_rc = decay_contribution(removal_c, x, θ, ℓ, n)
        @test dc_rc.alpha  == 5 // 2          # gamma_C*I_C = (1/2)*5
        @test dc_rc.gate   == 0 // 1          # I_C=5 > ell_C=2 -> no decay
        @test dc_rc.value  == 0 // 1
        @test dc_rc.reason == :sub_threshold_removal

        dc_rh = decay_contribution(removal_h, x, θ, ℓ, n)
        @test dc_rh.alpha  == 6 // 5          # gamma_H*I_H = (2/5)*3
        @test dc_rh.gate   == 1 // 1          # I_H=3 == ell_H=3 -> full decay
        @test dc_rh.value  == 6 // 5
        @test dc_rh.reason == :sub_threshold_removal

        dc_sc = decay_contribution(sampling_c, x, θ, ℓ, n)
        @test dc_sc.alpha  == 1 // 2          # chi_C*I_C = (1/10)*5
        @test dc_sc.gate   == 1 // 1          # SAMPLE: unconditional
        @test dc_sc.value  == 1 // 2
        @test dc_sc.reason == :sampling_hazard_full

        dc_sh = decay_contribution(sampling_h, x, θ, ℓ, n)
        @test dc_sh.alpha  == 9 // 10         # chi_H*I_H = (3/10)*3
        @test dc_sh.gate   == 1 // 1
        @test dc_sh.value  == 9 // 10
        @test dc_sh.reason == :sampling_hazard_full

        # Total lambda over the whole model: DEATH/SAMPLE terms only, the
        # four BIRTH marks and the four NEUTRAL demography marks contribute
        # nothing.
        @test total_decay(M, x, θ, ℓ, n) == 13 // 5
        @test total_decay(M, x, θ, ℓ, n) ==
            dc_rc.value + dc_rh.value + dc_sc.value + dc_sh.value
    end

    @info h2("MERS: boundary I_C == ell_C exactly, and I_C == ell_C+1 " *
             "(just above threshold)")
    @testset "MERS boundary" begin
        M = PhyloPOMP.MERS
        removal_c = find_event(M, :removal_c)
        x = (S_c = 1 // 1, I_c = 4 // 1, S_h = 1 // 1, I_h = 1 // 1)
        θ = (β_cc = 1 // 1, β_ch = 1 // 1, β_hc = 1 // 1, β_hh = 1 // 1,
             γ_c = 3 // 4, γ_h = 1 // 1, χ_c = 1 // 1, χ_h = 1 // 1,
             B_c = 1 // 1, B_h = 1 // 1, N_c = 1 // 1, N_h = 1 // 1)

        # I_C == ell_C == 4 exactly: gate = 1, value = gamma_C*I_C = 3.
        dc_at = decay_contribution(removal_c, x, θ, [4, 0], [4, 1])
        @test dc_at.gate  == 1 // 1
        @test dc_at.value == 3 // 1

        # I_C = ell_C+1 = 4 (ell_C=3): just above threshold -> gate = 0.
        dc_above = decay_contribution(removal_c, x, θ, [3, 0], [4, 1])
        @test dc_above.gate  == 0 // 1
        @test dc_above.value == 0 // 1
    end

    @info h2("SEIR: recovery (DEATH) + sampling (SAMPLE, r=(0,1)) match " *
             "the same lambda shape as MERS despite the differing r")
    @testset "SEIR concrete instance" begin
        S = PhyloPOMP.SEIR
        recovery = find_event(S, :recovery)
        sampling = find_event(S, :sampling)
        @test recovery.type == DEATH && recovery.r == [0, 0]
        @test sampling.type == SAMPLE && sampling.r == [0, 1]   # M03's finding

        # gamma=3/4, psi=1/5. Demes are (E,I); recovery/sampling both act on
        # I (deme 2).
        θ = (β = 1 // 1, σ = 1 // 1, γ = 3 // 4, ω = 1 // 1, ψ = 1 // 5,
             χ = 0 // 1, N = 100 // 1)

        # Above-threshold case: I=6, ell_I=3 (I > ell_I) -> recovery gate=0.
        x1 = (S = 1 // 1, E = 1 // 1, I = 6 // 1, R = 1 // 1)
        ℓ1, n1 = [1, 3], [1, 6]
        dc_rec1 = decay_contribution(recovery, x1, θ, ℓ1, n1)
        @test dc_rec1.alpha  == 9 // 2        # gamma*I = (3/4)*6
        @test dc_rec1.gate   == 0 // 1
        @test dc_rec1.value  == 0 // 1
        dc_smp1 = decay_contribution(sampling, x1, θ, ℓ1, n1)
        @test dc_smp1.alpha  == 6 // 5        # psi*I = (1/5)*6, UNCONDITIONAL
        @test dc_smp1.gate   == 1 // 1
        @test dc_smp1.value  == 6 // 5
        @test total_decay(S, x1, θ, ℓ1, n1) == 6 // 5   # 0 (recovery) + 6/5 (sampling)

        # At-threshold boundary: I=4, ell_I=4 (I == ell_I) -> recovery gate=1.
        x2 = (S = 1 // 1, E = 1 // 1, I = 4 // 1, R = 1 // 1)
        ℓ2, n2 = [1, 4], [1, 4]
        dc_rec2 = decay_contribution(recovery, x2, θ, ℓ2, n2)
        @test dc_rec2.gate  == 1 // 1
        @test dc_rec2.value == 3 // 1          # gamma*I = (3/4)*4
        dc_smp2 = decay_contribution(sampling, x2, θ, ℓ2, n2)
        @test dc_smp2.value == 4 // 5          # psi*I = (1/5)*4, still unconditional
        @test total_decay(S, x2, θ, ℓ2, n2) == 3 // 1 + 4 // 5

        # Just-above-threshold: I=5, ell_I=4 (I = ell_I+1) -> gate=0.
        x3 = (S = 1 // 1, E = 1 // 1, I = 5 // 1, R = 1 // 1)
        ℓ3, n3 = [1, 4], [1, 5]
        dc_rec3 = decay_contribution(recovery, x3, θ, ℓ3, n3)
        @test dc_rec3.gate  == 0 // 1
        @test dc_rec3.value == 0 // 1
    end

    @info h2("Scope boundaries: BIRTH/MIGRATION/NEUTRAL rejected or skipped")
    @testset "scope boundaries" begin
        S = PhyloPOMP.SEIR
        M = PhyloPOMP.MERS
        infection   = find_event(S, :infection)    # BIRTH
        progression = find_event(S, :progression)  # MIGRATION
        waning      = find_event(S, :waning)       # NEUTRAL
        θ = (β = 1 // 1, σ = 1 // 1, γ = 1 // 1, ω = 1 // 1, ψ = 1 // 1,
             χ = 0 // 1, N = 1 // 1)
        x = (S = 1 // 1, E = 1 // 1, I = 1 // 1, R = 1 // 1)

        @test_throws ArgumentError decay_contribution(infection, x, θ, [0, 0], [1, 1])
        @test_throws ArgumentError decay_contribution(progression, x, θ, [0, 0], [1, 1])
        @test_throws ArgumentError decay_contribution(waning, x, θ, [0, 0], [1, 1])

        # total_decay silently (correctly) skips BIRTH/MIGRATION/NEUTRAL --
        # confirm the MERS birth_c/birth_h/death_c/death_h NEUTRAL marks
        # contribute nothing, by comparing total_decay against a manual sum
        # over ONLY the DEATH/SAMPLE events.
        x_mers = (S_c = 1 // 1, I_c = 4 // 1, S_h = 1 // 1, I_h = 4 // 1)
        θ_mers = (β_cc = 1 // 1, β_ch = 1 // 1, β_hc = 1 // 1, β_hh = 1 // 1,
                  γ_c = 1 // 1, γ_h = 1 // 1, χ_c = 1 // 1, χ_h = 1 // 1,
                  B_c = 5 // 1, B_h = 5 // 1, N_c = 1 // 1, N_h = 1 // 1)
        ℓm, nm = [2, 2], [4, 4]
        manual = sum(decay_contribution(e, x_mers, θ_mers, ℓm, nm).value
                     for e in M.events if e.type in (DEATH, SAMPLE))
        @test total_decay(M, x_mers, θ_mers, ℓm, nm) == manual
        @test manual > 0 // 1   # sanity: not vacuously true
    end

    @info h2("DEATH invariant guard: ell[d] > n[d] throws")
    @testset "invariant guard" begin
        S = PhyloPOMP.SEIR
        recovery = find_event(S, :recovery)
        θ = (β = 1 // 1, σ = 1 // 1, γ = 1 // 1, ω = 1 // 1, ψ = 1 // 1,
             χ = 0 // 1, N = 1 // 1)
        x = (S = 1 // 1, E = 1 // 1, I = 3 // 1, R = 1 // 1)
        # ell_I=5 > n_I=3 violates the model invariant ell<=n.
        @test_throws ArgumentError decay_contribution(recovery, x, θ, [0, 5], [1, 3])
    end

end

end # module KliDecayTest
