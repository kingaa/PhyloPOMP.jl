"""
SI2R model verification test: runs the M02–M07 compiler pipeline on the
SI2R (superspreading) model — a model that has NEVER been seen by the
compiler before — and verifies every φ_u(s), Φ_u grouping, Chu–Vandermonde
identity, and decay term against the independently hand-derived values in
`si2r_model.qmd`.

WHY THIS TEST EXISTS: milestones M08–M12b validated the compiler against
SEIR and MERS, both of which were used during development. SI2R was
developed independently (in `/home/dislam/Desktop/Projects/efficiency/si2r/`)
and its filter equation hand-derived separately. If the compiler produces
the correct results for SI2R without having been trained on it, that is
genuine evidence that the framework generalizes — not circular validation.

REFERENCE: si2r_model.qmd, "Designing the filter" table (18 lines) and
"Filter equation" section (decay, driver, boost).
"""
module SI2RTest

using Test
using Random
using PhyloPOMP: SI2R, Event, BIRTH, DEATH, MIGRATION, SAMPLE, NEUTRAL,
    enumerate_saturations, kli_binomial_ratio, full_transitions,
    reduce_event_indicator, total_decay,
    IdentityTransition, InlineSameDemeTransition, CrossDemeTransition,
    ForkTransition, production_slots

by_kind(rts) = Dict(rt.kind => rt for rt in rts)

# ---------------------------------------------------------------------------
# Helper: look up an event by name
# ---------------------------------------------------------------------------
function si2r_event(name::Symbol)
    idx = findfirst(e -> e.name == name, SI2R.events)
    isnothing(idx) && error("SI2R event :$name not found")
    SI2R.events[idx]
end

# ---------------------------------------------------------------------------
# Helper: hand-derived binomial ratio formula from the QMD
# All formulas use exact Rational{Int}.
# ---------------------------------------------------------------------------
function qmd_binomial_ratio(n_L, ℓ_L, r_L, s_L, n_H, ℓ_H, r_H, s_H)
    # φ = C(n_L - ℓ_L, r_L - s_L)/C(n_L, r_L) * C(n_H - ℓ_H, r_H - s_H)/C(n_H, r_H)
    function safe_binom(a, b)
        (a < 0 || b < 0 || b > a) && return Rational{Int}(0)
        Rational{Int}(binomial(a, b))
    end
    d_L = safe_binom(n_L, r_L)
    d_H = safe_binom(n_H, r_H)
    (d_L == 0 || d_H == 0) && return Rational{Int}(0)
    safe_binom(n_L - ℓ_L, r_L - s_L) // d_L * safe_binom(n_H - ℓ_H, r_H - s_H) // d_H
end

@testset "SI2R Model Verification" begin

    # ===================================================================
    # Part 0: Structural checks — model declaration matches QMD
    # ===================================================================
    @testset "Part 0: Model structure" begin
        @test length(SI2R.events) == 9
        @test SI2R.name == :SI2R
        @test SI2R.demes == [:I_L, :I_H]

        # TL: BIRTH, from=1 (I_L), r=(2,0)
        tl = si2r_event(:TL)
        @test tl.type == BIRTH
        @test tl.from == 1
        @test tl.r == [2, 0]
        @test tl.regular == true

        # TH: BIRTH, from=2 (I_H), r=(1,1)
        th = si2r_event(:TH)
        @test th.type == BIRTH
        @test th.from == 2
        @test th.r == [1, 1]
        @test th.regular == true
        # New case always enters I_L (low-rate), not I_H -- the parent's
        # deme stays I_H. QMD: alpha_TH indicator{I_L'=I_L+1,S'=S-1}.
        @test th.Δ == [:S => -1, :I_L => +1]

        # L: MIGRATION, from=1 (I_L), r=(0,1)
        l = si2r_event(:L)
        @test l.type == MIGRATION
        @test l.from == 1
        @test l.r == [0, 1]
        @test l.regular == true

        # H: MIGRATION, from=2 (I_H), r=(1,0)
        h = si2r_event(:H)
        @test h.type == MIGRATION
        @test h.from == 2
        @test h.r == [1, 0]
        @test h.regular == true

        # RL: DEATH, from=1 (I_L), r=(0,0)  -- QMD line 12
        rl = si2r_event(:RL)
        @test rl.type == DEATH
        @test rl.from == 1
        @test rl.r == [0, 0]
        @test rl.regular == true

        # RH: DEATH, from=2 (I_H), r=(0,0)  -- QMD line 13
        rh = si2r_event(:RH)
        @test rh.type == DEATH
        @test rh.from == 2
        @test rh.r == [0, 0]
        @test rh.regular == true

        # W: NEUTRAL, r=(0,0)  -- QMD line 14
        w = si2r_event(:W)
        @test w.type == NEUTRAL
        @test w.r == [0, 0]
        @test w.regular == true

        # SL: SAMPLE, from=1, r=(1,0), singular  -- QMD lines 15-16
        sl = si2r_event(:SL)
        @test sl.type == SAMPLE
        @test sl.from == 1
        @test sl.r == [1, 0]
        @test sl.regular == false

        # SH: SAMPLE, from=2, r=(0,1), singular  -- QMD lines 17-18
        sh = si2r_event(:SH)
        @test sh.type == SAMPLE
        @test sh.from == 2
        @test sh.r == [0, 1]
        @test sh.regular == false
    end

    # ===================================================================
    # Part 0b: hazards alpha_u(t,x,x') for all nine marks at a concrete
    # state, checked against the QMD's alpha_u column verbatim. Nothing
    # elsewhere in this file exercises TL/TH/L/H/W's hazard closures --
    # phi_u/Phi_u/Chu-Vandermonde are all hazard-free -- so a wrong
    # hazard (e.g. TH using I_L instead of I_H) would otherwise pass
    # every other test in this file undetected.
    # ===================================================================
    @testset "Part 0b: hazards vs QMD alpha_u" begin
        S, I_L, I_H, R = 30, 12, 7, 15
        NN = S + I_L + I_H + R
        x = (S=S, I_L=I_L, I_H=I_H, R=R)
        θ = (β=1//4, κ=3//2, γ=1//5, ω=1//10, ψ=3//10, η_L=2//5, η_H=1//5, N=NN)

        @test si2r_event(:TL).hazard(x, θ) == θ.β * x.S * x.I_L / θ.N
        @test si2r_event(:TH).hazard(x, θ) == θ.κ * θ.β * x.S * x.I_H / θ.N
        @test si2r_event(:L).hazard(x, θ) == θ.η_L * x.I_L
        @test si2r_event(:H).hazard(x, θ) == θ.η_H * x.I_H
        @test si2r_event(:RL).hazard(x, θ) == θ.γ * x.I_L
        @test si2r_event(:RH).hazard(x, θ) == θ.γ * x.I_H
        @test si2r_event(:W).hazard(x, θ) == θ.ω * x.R
        @test si2r_event(:SL).hazard(x, θ) == θ.ψ * x.I_L
        @test si2r_event(:SH).hazard(x, θ) == θ.ψ * x.I_H

        # TH must use I_H, not I_L, and vice versa for TL -- the
        # discriminating check a copy-paste-from-MERS error would fail.
        @test si2r_event(:TL).hazard(x, θ) != θ.β * x.S * x.I_H / θ.N
        @test si2r_event(:TH).hazard(x, θ) != θ.κ * θ.β * x.S * x.I_L / θ.N
    end

    # ===================================================================
    # Part 1: φ_u(s) — exact binomial ratio verification
    #
    # For each BIRTH/MIGRATION event, at a concrete state, verify that
    # kli_binomial_ratio produces the same value as the QMD's formula.
    # We sweep many (n, ℓ) states using exact Rational{Int}.
    # ===================================================================
    @testset "Part 1: Binomial ratio φ_u(s) vs QMD" begin
        rng = Random.MersenneTwister(20240820)
        n_trials = 100
        total_assertions = 0

        for trial in 1:n_trials
            # Random state: I_L, I_H in [2, 50], ℓ_L ∈ [0, I_L], ℓ_H ∈ [0, I_H]
            I_L = rand(rng, 2:50)
            I_H = rand(rng, 2:50)
            ℓ_L = rand(rng, 0:I_L)
            ℓ_H = rand(rng, 0:I_H)
            n = [I_L, I_H]
            ℓ = [ℓ_L, ℓ_H]

            # Non-vacuity guard: enumerate_saturations must produce exactly
            # prod(min(r_d,ell_d)+1) saturations for every event checked
            # below -- otherwise a short/empty enumeration would make the
            # `for s in sats` loops silently assert nothing.
            expect_nsats(r) = prod(min(r[d], ℓ[d]) + 1 for d in 1:2)

            # --- TL: r=(2,0), 3 saturations: s=(0,0), (1,0), (2,0) ---
            # QMD lines 1-3
            tl = si2r_event(:TL)
            sats = enumerate_saturations(tl.r, ℓ)
            @test length(sats) == expect_nsats(tl.r)
            for s in sats
                compiler_φ = kli_binomial_ratio(tl, s, ℓ, n)
                qmd_φ = qmd_binomial_ratio(I_L, ℓ_L, 2, s[1], I_H, ℓ_H, 0, s[2])
                @test compiler_φ == qmd_φ
                total_assertions += 1
            end

            # --- TH: r=(1,1), 4 saturations: (0,0), (1,0), (0,1), (1,1) ---
            # QMD lines 4-7
            th = si2r_event(:TH)
            sats = enumerate_saturations(th.r, ℓ)
            @test length(sats) == expect_nsats(th.r)
            for s in sats
                compiler_φ = kli_binomial_ratio(th, s, ℓ, n)
                qmd_φ = qmd_binomial_ratio(I_L, ℓ_L, 1, s[1], I_H, ℓ_H, 1, s[2])
                @test compiler_φ == qmd_φ
                total_assertions += 1
            end

            # --- L: r=(0,1), 2 saturations: (0,0), (0,1) ---
            # QMD lines 8-9
            l = si2r_event(:L)
            sats = enumerate_saturations(l.r, ℓ)
            @test length(sats) == expect_nsats(l.r)
            for s in sats
                compiler_φ = kli_binomial_ratio(l, s, ℓ, n)
                qmd_φ = qmd_binomial_ratio(I_L, ℓ_L, 0, s[1], I_H, ℓ_H, 1, s[2])
                @test compiler_φ == qmd_φ
                total_assertions += 1
            end

            # --- H: r=(1,0), 2 saturations: (0,0), (1,0) ---
            # QMD lines 10-11
            h = si2r_event(:H)
            sats = enumerate_saturations(h.r, ℓ)
            @test length(sats) == expect_nsats(h.r)
            for s in sats
                compiler_φ = kli_binomial_ratio(h, s, ℓ, n)
                qmd_φ = qmd_binomial_ratio(I_L, ℓ_L, 1, s[1], I_H, ℓ_H, 0, s[2])
                @test compiler_φ == qmd_φ
                total_assertions += 1
            end

            # --- RL/RH/W: r=(0,0), 1 saturation s=(0,0), φ=1 trivially ---
            # QMD lines 12-14
            for name in (:RL, :RH, :W)
                ev = si2r_event(name)
                sats = enumerate_saturations(ev.r, ℓ)
                @test length(sats) == 1
                @test kli_binomial_ratio(ev, sats[1], ℓ, n) == 1 // 1
                total_assertions += 1
            end

            # --- SL: r=(1,0), 2 saturations: (0,0), (1,0) --- QMD lines 15-16
            sl = si2r_event(:SL)
            sats = enumerate_saturations(sl.r, ℓ)
            @test length(sats) == expect_nsats(sl.r)
            for s in sats
                compiler_φ = kli_binomial_ratio(sl, s, ℓ, n)
                qmd_φ = qmd_binomial_ratio(I_L, ℓ_L, 1, s[1], I_H, ℓ_H, 0, s[2])
                @test compiler_φ == qmd_φ
                total_assertions += 1
            end

            # --- SH: r=(0,1), 2 saturations: (0,0), (0,1) --- QMD lines 17-18
            sh = si2r_event(:SH)
            sats = enumerate_saturations(sh.r, ℓ)
            @test length(sats) == expect_nsats(sh.r)
            for s in sats
                compiler_φ = kli_binomial_ratio(sh, s, ℓ, n)
                qmd_φ = qmd_binomial_ratio(I_L, ℓ_L, 0, s[1], I_H, ℓ_H, 1, s[2])
                @test compiler_φ == qmd_φ
                total_assertions += 1
            end
        end
        @info "Part 1: $total_assertions φ_u assertions passed"
    end

    # ===================================================================
    # Part 2: Transition classification
    #
    # Verify that full_transitions classifies each saturation correctly:
    #   TL: s=(0,0)→Identity, s=(1,0)→InlineSameDeme, s=(2,0)→Fork
    #   TH: s=(0,0)→Identity, s=(0,1)→InlineSameDeme, s=(1,0)→CrossDeme, s=(1,1)→Fork
    #   L:  s=(0,0)→Identity, s=(0,1)→CrossDeme
    #   H:  s=(0,0)→Identity, s=(1,0)→CrossDeme
    # ===================================================================
    @testset "Part 2: Transition classification" begin
        # Use a specific state where all saturations are feasible
        I_L, I_H = 10, 8
        ℓ_L, ℓ_H = 3, 2
        n = [I_L, I_H]
        ℓ = [ℓ_L, ℓ_H]

        # --- TL: r=(2,0), from=1 (I_L) ---
        tl = si2r_event(:TL)
        trans_tl = full_transitions(tl, ℓ, n)
        @test length(trans_tl) == 3  # QMD lines 1-3

        # s=(0,0): Identity
        @test trans_tl[1] isa IdentityTransition
        @test trans_tl[1].s == [0, 0]

        # s=(1,0): InlineSameDeme (one slot in deme 1 = I_L = from deme)
        @test trans_tl[2] isa InlineSameDemeTransition
        @test trans_tl[2].s == [1, 0]
        @test trans_tl[2].deme == 1  # same deme as from

        # s=(2,0): Fork (both slots filled by tracked lineages, same deme)
        @test trans_tl[3] isa ForkTransition
        @test trans_tl[3].s == [2, 0]

        # --- TH: r=(1,1), from=2 (I_H) ---
        th = si2r_event(:TH)
        trans_th = full_transitions(th, ℓ, n)
        @test length(trans_th) == 4  # QMD lines 4-7

        # s=(0,0): Identity
        @test trans_th[1] isa IdentityTransition
        @test trans_th[1].s == [0, 0]

        # Saturation order: (0,0), (1,0), (0,1), (1,1)
        # s=(1,0): CrossDeme (slot in deme 1 = I_L ≠ from(=2))
        # QMD line 6: swap(HL) = cross-deme H→L
        @test trans_th[2] isa CrossDemeTransition
        @test trans_th[2].s == [1, 0]

        # s=(0,1): InlineSameDeme (slot in deme 2 = I_H = from deme)
        # QMD line 5: swap(HH) = same-deme inline
        @test trans_th[3] isa InlineSameDemeTransition
        @test trans_th[3].s == [0, 1]
        @test trans_th[3].deme == 2  # same as from

        # s=(1,1): Fork (both slots filled)
        # QMD line 7: fork(HLH)
        @test trans_th[4] isa ForkTransition
        @test trans_th[4].s == [1, 1]

        # --- L: r=(0,1), from=1 (I_L) ---
        l = si2r_event(:L)
        trans_l = full_transitions(l, ℓ, n)
        @test length(trans_l) == 2  # QMD lines 8-9

        # s=(0,0): Identity
        @test trans_l[1] isa IdentityTransition

        # s=(0,1): CrossDeme (slot in deme 2 = I_H ≠ from(=1))
        # QMD line 9: swap(LH) = cross-deme L→H
        @test trans_l[2] isa CrossDemeTransition
        @test trans_l[2].s == [0, 1]

        # --- H: r=(1,0), from=2 (I_H) ---
        h = si2r_event(:H)
        trans_h = full_transitions(h, ℓ, n)
        @test length(trans_h) == 2  # QMD lines 10-11

        # s=(0,0): Identity
        @test trans_h[1] isa IdentityTransition

        # s=(1,0): CrossDeme (slot in deme 1 = I_L ≠ from(=2))
        # QMD line 11: swap(HL) = cross-deme H→L
        @test trans_h[2] isa CrossDemeTransition
        @test trans_h[2].s == [1, 0]
    end

    # ===================================================================
    # Part 2b: reduce_event_indicator collapse structure (QMD lines 1+2,
    # 4+5) — the specific claim "lines 1&2 combine" / "lines 4&5 combine"
    # from the QMD's Filter equation section, checked two ways:
    #  (i)  GROUPING: at a concrete state, reduce_event_indicator must
    #       put exactly the Identity+InlineSameDeme pair into one :noop
    #       group (not 3 separate groups), leaving Fork/CrossDeme alone
    #       — mirrors kli_reduce_test.jl's TCC/THC pattern.
    #  (ii) WEIGHTED CLOSED FORM: reduce_event_indicator's Φ_noop is the
    #       UNWEIGHTED sum phi_id + phi_inline (M09 Finding 2) — it is
    #       NOT itself equal to the QMD's collapsed formula, which
    #       additionally weights the InlineSameDeme term by the C(ell,1)
    #       "any of the ell tracked lineages" multiplicity. Reproduce
    #       that weighting explicitly and check against QMD Eq (132)/(135)
    #       in exact Rational{Int}, swept over random states.
    # ===================================================================
    @testset "Part 2b: reduce_event_indicator collapse (lines 1+2, 4+5)" begin
        # --- (i) grouping structure at the Part-2 concrete state ---
        I_L, I_H = 10, 8
        ℓ_L, ℓ_H = 3, 2
        n = [I_L, I_H]
        ℓ = [ℓ_L, ℓ_H]

        tl = si2r_event(:TL)
        rts_tl = reduce_event_indicator(full_transitions(tl, ℓ, n))
        @test length(rts_tl) == 2  # QMD: {1,2} collapse, 3 alone
        g_tl = by_kind(rts_tl)
        @test Set(typeof(t) for t in g_tl[:noop].transitions) ==
              Set([IdentityTransition, InlineSameDemeTransition])
        @test Set(typeof(t) for t in g_tl[:fork].transitions) == Set([ForkTransition])

        th = si2r_event(:TH)
        rts_th = reduce_event_indicator(full_transitions(th, ℓ, n))
        @test length(rts_th) == 3  # QMD: {4,5} collapse, 6 and 7 alone
        g_th = by_kind(rts_th)
        @test Set(typeof(t) for t in g_th[:noop].transitions) ==
              Set([IdentityTransition, InlineSameDemeTransition])
        @test Set(typeof(t) for t in g_th[:cross].transitions) == Set([CrossDemeTransition])
        @test Set(typeof(t) for t in g_th[:fork].transitions) == Set([ForkTransition])

        l = si2r_event(:L)
        rts_l = reduce_event_indicator(full_transitions(l, ℓ, n))
        @test length(rts_l) == 2  # QMD: 8 alone (no InlineSameDeme possible, r_L=0), 9 alone
        g_l = by_kind(rts_l)
        @test Set(typeof(t) for t in g_l[:noop].transitions) == Set([IdentityTransition])
        @test Set(typeof(t) for t in g_l[:cross].transitions) == Set([CrossDemeTransition])

        h = si2r_event(:H)
        rts_h = reduce_event_indicator(full_transitions(h, ℓ, n))
        @test length(rts_h) == 2  # QMD: 10 alone, 11 alone
        g_h = by_kind(rts_h)
        @test Set(typeof(t) for t in g_h[:noop].transitions) == Set([IdentityTransition])
        @test Set(typeof(t) for t in g_h[:cross].transitions) == Set([CrossDemeTransition])

        # --- (ii) weighted closed-form check, swept ---
        # Safe extraction: InlineSameDeme is only reachable when ell>=1 at
        # that deme, so it may be absent from the :noop group entirely
        # (ell=0) -- treat as phi=0 in that case (its C(ell,1)=0 weight
        # would zero it out anyway, but `only(...)` on an empty generator
        # throws, so this must be handled explicitly, not left to `only`).
        function phi_of(group, ::Type{T}) where {T}
            members = [t for t in group.transitions if t isa T]
            isempty(members) ? Rational{Int}(0) : only(members).phi
        end

        rng = Random.MersenneTwister(20240825)
        n_trials = 150
        assertions = 0
        for trial in 1:n_trials
            I_L = rand(rng, 2:60)
            I_H = rand(rng, 2:60)
            ℓ_L = rand(rng, 0:I_L)
            ℓ_H = rand(rng, 0:I_H)
            n = [I_L, I_H]
            ℓ = [ℓ_L, ℓ_H]

            # QMD Eq (132): lines 1&2 -> [1 - C(ℓ_L,2)/C(I_L,2)]
            tl = si2r_event(:TL)
            rts = reduce_event_indicator(full_transitions(tl, ℓ, n))
            g = by_kind(rts)
            φ_id = phi_of(g[:noop], IdentityTransition)
            φ_inl = phi_of(g[:noop], InlineSameDemeTransition)
            weighted = φ_id + ℓ_L * φ_inl
            # binomial(I_L,2) is 0 only for I_L<2, in which case the QMD's
            # RHS gate 1{I_L>=ell_L} degenerates; sampled range is I_L>=2.
            expected = 1 - Rational{Int}(binomial(ℓ_L, 2)) // binomial(I_L, 2)
            @test weighted == expected
            assertions += 1

            # QMD Eq (135): lines 4&5 -> [1 - ℓ_L/I_L]  (H-factor cancels)
            th = si2r_event(:TH)
            rts2 = reduce_event_indicator(full_transitions(th, ℓ, n))
            g2 = by_kind(rts2)
            φ_id2 = phi_of(g2[:noop], IdentityTransition)
            φ_inl2 = phi_of(g2[:noop], InlineSameDemeTransition)
            weighted2 = φ_id2 + ℓ_H * φ_inl2
            expected2 = 1 - ℓ_L // I_L
            @test weighted2 == expected2
            assertions += 1
        end
        @info "Part 2b: $assertions weighted-collapse (QMD Eq 132/135) assertions passed"
    end

    # ===================================================================
    # Part 3: Chu–Vandermonde identity: Σ φ_u(s) * C(ℓ,s) = 1
    # Verified for every BIRTH/MIGRATION event across random states.
    # ===================================================================
    @testset "Part 3: Chu-Vandermonde identity" begin
        rng = Random.MersenneTwister(20240821)
        n_trials = 200
        assertions = 0

        for trial in 1:n_trials
            I_L = rand(rng, 2:100)
            I_H = rand(rng, 2:100)
            ℓ_L = rand(rng, 0:I_L)
            ℓ_H = rand(rng, 0:I_H)
            n = [I_L, I_H]
            ℓ = [ℓ_L, ℓ_H]

            for event_name in (:TL, :TH, :L, :H)
                ev = si2r_event(event_name)
                sats = enumerate_saturations(ev.r, ℓ)
                cv_sum = Rational{Int}(0)
                for s in sats
                    φ = kli_binomial_ratio(ev, s, ℓ, n)
                    cv_sum += φ * prod(binomial(ℓ[d], s[d]) for d in 1:2)
                end
                @test cv_sum == 1 // 1
                assertions += 1
            end
        end
        @info "Part 3: $assertions Chu-Vandermonde assertions passed"
    end

    # ===================================================================
    # Part 4: m-reduction — verify Φ grouping consistency
    #
    # For each event's reduced transitions, verify:
    # (a) Φ = Σ φ(member transitions) — the Φ invariant
    # (b) Each member's φ matches qmd_binomial_ratio — φ correctness
    # (c) Correct group classification (noop/cross/fork)
    # (d) Σ Φ across all groups = Σ φ across all transitions (no loss)
    # ===================================================================
    @testset "Part 4: m-reduction and Φ_u grouping" begin
        rng = Random.MersenneTwister(20240822)
        n_trials = 150
        assertions = 0

        for trial in 1:n_trials
            I_L = rand(rng, 3:80)
            I_H = rand(rng, 3:80)
            ℓ_L = rand(rng, 0:I_L)
            ℓ_H = rand(rng, 0:I_H)
            n = [I_L, I_H]
            ℓ = [ℓ_L, ℓ_H]

            for event_name in (:TL, :TH, :L, :H)
                ev = si2r_event(event_name)
                trans = full_transitions(ev, ℓ, n)
                reduced = reduce_event_indicator(trans)

                # (a) Φ invariant: each group's Φ = sum of member φ values
                for rt in reduced
                    member_sum = sum(t.phi for t in rt.transitions)
                    @test rt.Φ == member_sum
                    assertions += 1
                end

                # (b) Each member transition's φ matches QMD formula
                for rt in reduced
                    for t in rt.transitions
                        r = production_slots(ev)
                        expected_phi = qmd_binomial_ratio(
                            n[1], ℓ[1], r[1], t.s[1],
                            n[2], ℓ[2], r[2], t.s[2])
                        @test t.phi == expected_phi
                        assertions += 1
                    end
                end

                # (c) No loss: sum of Φ across groups = sum of φ across all transitions
                total_Φ = sum(rt.Φ for rt in reduced)
                total_φ = sum(t.phi for t in trans)
                @test total_Φ == total_φ
                assertions += 1
            end
        end
        @info "Part 4: $assertions Φ_u grouping assertions passed"
    end

    # ===================================================================
    # Part 5: Decay λ — matches QMD Eq. (decay)
    #
    # QMD: λ(t,x,d) = [α_SL + α_SH]
    #                 + α_RL * 1{I_L ≤ ℓ_L} + α_RH * 1{I_H ≤ ℓ_H}
    #
    # For α_SL = ψ*I_L, α_SH = ψ*I_H, α_RL = γ*I_L, α_RH = γ*I_H
    # ===================================================================
    @testset "Part 5: Decay λ vs QMD" begin
        rng = Random.MersenneTwister(20240823)
        n_trials = 100
        assertions = 0
        boundary_L_hits = 0  # I_L <= ℓ_L, forces the RL leftover branch
        boundary_H_hits = 0  # I_H <= ℓ_H, forces the RH leftover branch

        for trial in 1:n_trials
            # Random state with rational parameters
            S = Rational{Int}(rand(rng, 1:100))
            I_L = rand(rng, 1:50)
            I_H = rand(rng, 1:50)
            R = Rational{Int}(rand(rng, 0:30))
            ℓ_L = rand(rng, 0:I_L)
            ℓ_H = rand(rng, 0:I_H)
            # Force the I<=ell boundary (decay leftover) to fire on some
            # trials -- with I in 1:50 drawn uniformly it would otherwise
            # only hit by chance (~6%/trial); guarantee non-vacuity.
            trial % 5 == 0 && (ℓ_L = I_L)
            trial % 7 == 0 && (ℓ_H = I_H)
            I_L <= ℓ_L && (boundary_L_hits += 1)
            I_H <= ℓ_H && (boundary_H_hits += 1)

            # SI2R parameters (all positive rationals)
            β = Rational{Int}(rand(rng, 1:10))
            κ = Rational{Int}(rand(rng, 1:5))
            γ = Rational{Int}(rand(rng, 1:10)) // 10
            ω = Rational{Int}(rand(rng, 1:10)) // 10
            ψ = Rational{Int}(rand(rng, 1:10)) // 10
            η_L = Rational{Int}(rand(rng, 1:10)) // 10
            η_H = Rational{Int}(rand(rng, 1:10)) // 10
            NN = Rational{Int}(S + I_L + I_H + R)

            x = (S=S, I_L=Rational{Int}(I_L), I_H=Rational{Int}(I_H), R=R)
            θ = (β=β, κ=κ, γ=γ, ω=ω, ψ=ψ, η_L=η_L, η_H=η_H, N=NN)
            n = [I_L, I_H]
            ℓ = [ℓ_L, ℓ_H]

            # Hand-computed decay from QMD Eq. (decay):
            # λ = ψ*I_L + ψ*I_H + γ*I_L*1{I_L ≤ ℓ_L} + γ*I_H*1{I_H ≤ ℓ_H}
            expected_decay = ψ * I_L + ψ * I_H +
                γ * I_L * (I_L <= ℓ_L ? 1 : 0) +
                γ * I_H * (I_H <= ℓ_H ? 1 : 0)

            compiler_decay = total_decay(SI2R, x, θ, ℓ, n)
            @test compiler_decay == expected_decay
            assertions += 1
        end
        # Non-vacuity guard: the I<=ell leftover branches must actually
        # have fired at least once each, not merely be dead code that
        # happens to agree with the ordinary branch every time it's hit.
        @test boundary_L_hits >= 1
        @test boundary_H_hits >= 1
        @info "Part 5: $assertions decay assertions passed " *
              "($boundary_L_hits I_L<=ℓ_L boundary hits, $boundary_H_hits I_H<=ℓ_H)"
    end

    # ===================================================================
    # Part 6: TL Fork = Kingman/Moran (same as TCC in M12b)
    #
    # TL has r=(2,0): same-deme fork with N=I_L.
    # φ_fork = 1/C(I_L,2), independent of ℓ_L.
    # This is the Kingman/Moran result, now confirmed on a THIRD model.
    # ===================================================================
    @testset "Part 6: TL Fork = Kingman/Moran" begin
        rng = Random.MersenneTwister(20240824)
        n_trials = 200
        assertions = 0

        tl = si2r_event(:TL)
        for trial in 1:n_trials
            I_L = rand(rng, 2:200)
            I_H = rand(rng, 0:100)
            # Fork requires ℓ_L ≥ 2 (need two tracked lineages to fill
            # both production slots)
            ℓ_L = rand(rng, 2:I_L)
            ℓ_H = rand(rng, 0:I_H)
            n = [I_L, I_H]
            ℓ = [ℓ_L, ℓ_H]

            trans = full_transitions(tl, ℓ, n)
            fork_trans = filter(t -> t isa ForkTransition, trans)
            @test length(fork_trans) == 1

            expected_moran = 1 // Rational{Int}(binomial(I_L, 2))
            @test fork_trans[1].phi == expected_moran

            # Independence from ℓ_L, ℓ_H, I_H: verified implicitly by
            # sweeping many different values and getting the same formula.
            assertions += 1
        end
        @info "Part 6: $assertions Moran/Kingman assertions passed (SI2R TL)"
    end

    # ===================================================================
    # Part 6b: TH Fork is an ORDERED pair, 1/(I_L*I_H) -- discriminates
    # from TL's UNORDERED same-deme pair, 1/C(I_L,2) (Part 6). TH's two
    # production slots sit in DIFFERENT demes (L and H), so there is no
    # "which of 2 slots" symmetry to collapse away — QMD boost line 167
    # states this denominator as I_L'*I_H' directly (not C(I_L+I_H,2) or
    # any single-N binomial). This is exactly the shape a compiler that
    # merely pattern-matched MERS's TCC/THH (same-deme, r=(2,0)/(0,2))
    # could get wrong when it meets a cross-deme r=(1,1) fork instead.
    # ===================================================================
    @testset "Part 6b: TH Fork = ordered pair 1/(I_L*I_H)" begin
        rng = Random.MersenneTwister(20240826)
        n_trials = 200
        assertions = 0

        th = si2r_event(:TH)
        for trial in 1:n_trials
            I_L = rand(rng, 1:150)
            I_H = rand(rng, 1:150)
            ℓ_L = rand(rng, 1:I_L)  # need >=1 tracked lineage in each deme
            ℓ_H = rand(rng, 1:I_H)  # to fill both of TH's fork slots
            n = [I_L, I_H]
            ℓ = [ℓ_L, ℓ_H]

            trans = full_transitions(th, ℓ, n)
            fork_trans = filter(t -> t isa ForkTransition, trans)
            @test length(fork_trans) == 1

            expected = 1 // (Rational{Int}(I_L) * Rational{Int}(I_H))
            @test fork_trans[1].phi == expected
            assertions += 1
        end
        @info "Part 6b: $assertions TH ordered-fork assertions passed"
    end

    # ===================================================================
    # Part 7: QMD's SIMPLIFIED closed forms (boost Eq lines 157-162),
    # independent of `qmd_binomial_ratio`'s product-form reimplementation
    # used in Parts 1/3/4. Part 1 compares two implementations of the
    # SAME product formula; this part instead checks `kli_binomial_ratio`
    # against algebraically-simplified expressions the QMD states
    # directly, giving genuinely independent verification content.
    # ===================================================================
    @testset "Part 7: QMD simplified closed forms (boost Eq)" begin
        rng = Random.MersenneTwister(20240827)
        n_trials = 200
        assertions = 0

        for trial in 1:n_trials
            I_L = rand(rng, 2:100)
            I_H = rand(rng, 2:100)
            ℓ_L = rand(rng, 0:I_L)
            ℓ_H = rand(rng, 0:I_H)
            n = [I_L, I_H]
            ℓ = [ℓ_L, ℓ_H]

            # TH cross (s=(1,0), QMD line 6) == (1/I_L)*(1 - ℓ_H/I_H)
            # (boost line 158's 1/I_L' * [1-ℓ_H'/I_H'] factor, restated
            # as a phi_u value rather than a boost/proposal quantity)
            th = si2r_event(:TH)
            φ_th_cross = kli_binomial_ratio(th, [1, 0], ℓ, n)
            @test φ_th_cross == (1 // I_L) * (1 - ℓ_H // I_H)
            assertions += 1

            # L cross (s=(0,1), QMD line 9) == 1/I_H  (boost line 160)
            l = si2r_event(:L)
            φ_l_cross = kli_binomial_ratio(l, [0, 1], ℓ, n)
            @test φ_l_cross == 1 // I_H
            assertions += 1

            # L noop (s=(0,0), QMD line 8) == 1 - ℓ_H/I_H  (boost line 159's
            # π_L^∅ = 1-ℓ_L/I_L is the PROPOSAL, not phi_u itself; phi_u at
            # r=(0,1),s=(0,0) reduces to the H-deme factor alone since
            # r_L=0 contributes an unconditional factor of 1)
            φ_l_noop = kli_binomial_ratio(l, [0, 0], ℓ, n)
            @test φ_l_noop == 1 - ℓ_H // I_H
            assertions += 1

            # H cross (s=(1,0), QMD line 11) == 1/I_L  (boost line 162)
            h = si2r_event(:H)
            φ_h_cross = kli_binomial_ratio(h, [1, 0], ℓ, n)
            @test φ_h_cross == 1 // I_L
            assertions += 1

            # H noop (s=(0,0), QMD line 10) == 1 - ℓ_L/I_L
            φ_h_noop = kli_binomial_ratio(h, [0, 0], ℓ, n)
            @test φ_h_noop == 1 - ℓ_L // I_L
            assertions += 1
        end
        @info "Part 7: $assertions QMD simplified-closed-form assertions passed"
    end

end  # @testset "SI2R Model Verification"

end # module SI2RTest
