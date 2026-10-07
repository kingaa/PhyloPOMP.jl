"""
Generic regular step (`regular_step!` with the naive proposal) for MERS and SEIR.

1. Slot structure: `kli_slots` gives the 12 MERS and 6 SEIR slots of the compiled regular parts.
2. Per-outcome weights of `apply_move!` and selection factors of `kli_select` against the closed forms
   (`1 - C(ℓ,2)/C(n',2)`, `1 - ℓ_other/n'_other`, `(1 - ℓ'_src/n_src)/n'_dst`, `no_move_share`), on an
   (ℓ, n) grid that includes ℓ = n.
3. `kli_decay` against `compiled_decay`.
4. Isolated `regular_step!` against `mers_compiled_regular_part!` and `compiled_regular_part!`: same seed, same
   final ll, state, coloring, and the same next draw of the global RNG.
5. Seed-matched Np=100 log-likelihoods of the compiled filters with `generic_regular = true` (and with both generic
   parts) against the default filters, on simulated genealogies.
6. A coloring whose demeset has a different number of demes than the model is rejected with `ArgumentError`
   by `singular_update!`, `regular_step!`, `kli_select` and `apply_move!`.
"""
module KliRegularStepTest

import ..Main: h1, h2

@info h1("Generic regular step (MERS, SEIR)")

using Test
using PhyloPOMP
using PhyloPOMP.NaiveMERS
using PhyloPOMP.NaiveSEIR
using PhyloPOMP: Coloring, ell, plant!, MERS, Event, BIRTH, MIGRATION, regular_step!, kli_slots,
    kli_select, kli_rates!, kli_decay, apply_move!, apply_pop, compiled_decay, no_move_share
using Random: MersenneTwister, seed!, default_rng

const SE = PhyloPOMP.SEIR
Camel, Human = NaiveMERS.Camel, NaiveMERS.Human
Infec, Expos = NaiveSEIR.Infec, NaiveSEIR.Expos
find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]
c2(x) = x * (x - 1) // 2

function mers_cols(lc, lh)
    cols = Coloring(NaiveMERS.Demes)
    foreach(b -> plant!(cols, Camel, b), 1:lc)
    foreach(b -> plant!(cols, Human, 100 + b), 1:lh)
    cols
end

function seir_cols(lE, lI)
    cols = Coloring(NaiveSEIR.Demes)
    foreach(b -> plant!(cols, Expos, b), 1:lE)
    foreach(b -> plant!(cols, Infec, 100 + b), 1:lI)
    cols
end

@testset verbose=true "generic regular step" begin

    @testset "slot structure" begin
        idx(name) = findfirst(e -> e.name == name, MERS.events)
        @test kli_slots(MERS) == [
            (idx(:transmission_cc), 0), (idx(:transmission_hh), 0),
            (idx(:transmission_hc), 0), (idx(:transmission_hc), 2),
            (idx(:transmission_ch), 0), (idx(:transmission_ch), 1),
            (idx(:removal_c), 0), (idx(:removal_h), 0),
            (idx(:birth_c), 0), (idx(:birth_h), 0), (idx(:death_c), 0), (idx(:death_h), 0)]
        sidx(name) = findfirst(e -> e.name == name, SE.events)
        @test kli_slots(SE) == [(sidx(:infection), 0), (sidx(:infection), 1),
                                (sidx(:progression), 0), (sidx(:progression), 2),
                                (sidx(:recovery), 0), (sidx(:waning), 0)]
        # a birth that does not continue its parent has no naive selection rule here
        orphan = Event(:orphan, Pair{Symbol,Int}[], (x, θ) -> 1.0, [0, 2], BIRTH, 1, [2], true, false)
        @test_throws ArgumentError kli_slots(orphan)
    end

    @info h2("apply_move! weights and kli_select factors against the closed forms, ℓ = n included")
    @testset "per-outcome weights (MERS)" begin
        tcc, thh = find_event(MERS, :transmission_cc), find_event(MERS, :transmission_hh)
        thc, tch = find_event(MERS, :transmission_hc), find_event(MERS, :transmission_ch)
        rc, rh = find_event(MERS, :removal_c), find_event(MERS, :removal_h)
        neutral = [find_event(MERS, s) for s in (:birth_c, :birth_h, :death_c, :death_h)]
        rng = MersenneTwister(1)
        close(w, target) = isapprox(w, log(Float64(target)); atol = 1e-12)
        bad = 0; total = 0; nboundary = 0
        check(ok) = (total += 1; bad += !ok)
        for Ic in 1:6, Ih in 1:6, lc in 0:Ic, lh in 0:Ih
            x = (S_c = 20, I_c = Ic, S_h = 15, I_h = Ih)
            # within-deme births: one slot, no lineage moves, weight 1 - C(ℓ,2)/C(n',2)
            for (ev, l, n1) in ((tcc, lc, Ic + 1), (thh, lh, Ih + 1))
                cols = mers_cols(lc, lh)
                check(kli_select(ev, cols, x, MERS) == [1.0])
                w = apply_move!(cols, ev, 0, apply_pop(x, ev), MERS; rng)
                check(close(w, 1 - c2(l) // c2(n1)) && ell(cols) == (lc, lh))
            end
            # cross-deme births: (src, dst) = (camel, human) for THC, (human, camel) for TCH
            for (ev, ls, ld, ns, nd, dst) in ((thc, lc, lh, Ic, Ih, 2), (tch, lh, lc, Ih, Ic, 1))
                x1 = apply_pop(x, ev)
                cols = mers_cols(lc, lh)
                p = kli_select(ev, cols, x, MERS)
                check(isapprox(p[1], no_move_share(ls, ld, ns, nd); atol = 1e-14) && p[2] == 1 - p[1])
                w0 = apply_move!(cols, ev, 0, x1, MERS; rng)
                T0 = 1 - ld // (nd + 1)
                check(close(w0, T0) && ell(cols) == (lc, lh))
                check(p[1] > 0 && isfinite(w0))
                ls == ns && (nboundary += 1)
                T1 = 0 // 1
                if ls ≥ 1
                    w1 = apply_move!(cols, ev, dst, x1, MERS; rng)
                    T1 = ls * (1 - (ls - 1) // ns) // (nd + 1)
                    moved = dst == 2 ? (lc - 1, lh + 1) : (lc + 1, lh - 1)
                    check(close(w1, T1) && ell(cols) == moved)
                    # the naive selection is the target share, so the full weights e^w/π are equal
                    check(isapprox(exp(w0) / p[1], exp(w1) / p[2]; rtol = 1e-12))
                    check(isapprox(exp(w0) / p[1], Float64(T0 + T1); rtol = 1e-12))
                else
                    check(p[2] == 0)
                end
            end
            # removal: π is the chance that the dying host is untracked; no coloring change
            for (ev, l, n) in ((rc, lc, Ic), (rh, lh, Ih))
                cols = mers_cols(lc, lh)
                check(kli_select(ev, cols, x, MERS) == [n > l ? 1 - l / n : 0.0])
                check(apply_move!(cols, ev, 0, apply_pop(x, ev), MERS; rng) == 0.0)
            end
            for ev in neutral
                cols = mers_cols(lc, lh)
                check(kli_select(ev, cols, x, MERS) == [1.0])
                check(apply_move!(cols, ev, 0, apply_pop(x, ev), MERS; rng) == 0.0)
            end
        end
        @info "MERS per-outcome checks: $total checks, $bad failures ($nboundary states with every source host tracked)"
        @test total > 5_000
        @test nboundary > 50
        @test bad == 0
    end

    @testset "per-outcome weights (SEIR)" begin
        inf, prog = find_event(SE, :infection), find_event(SE, :progression)
        rec = find_event(SE, :recovery)
        rng = MersenneTwister(2)
        close(w, target) = isapprox(w, log(Float64(target)); atol = 1e-12)
        bad = 0; total = 0; nboundary = 0
        check(ok) = (total += 1; bad += !ok)
        for E in 0:6, I in 1:6, lE in 0:E, lI in 0:I
            x = (S = 30, E = E, I = I, R = 4)
            # infection: BIRTH from I into E
            x1 = apply_pop(x, inf)
            cols = seir_cols(lE, lI)
            p = kli_select(inf, cols, x, SE)
            check(isapprox(p[1], no_move_share(lI, lE, I, E); atol = 1e-14) && p[1] > 0)
            w0 = apply_move!(cols, inf, 0, x1, SE; rng)
            check(close(w0, 1 - lE // (E + 1)) && ell(cols) == (lE, lI))
            lI == I && (nboundary += 1)
            if lI ≥ 1
                w1 = apply_move!(cols, inf, 1, x1, SE; rng)
                check(close(w1, lI * (1 - (lI - 1) // I) // (E + 1)) && ell(cols) == (lE + 1, lI - 1))
            end
            # progression: MIGRATION from E to I; the moving host is the parent
            E == 0 && continue
            x1 = apply_pop(x, prog)
            cols = seir_cols(lE, lI)
            p = kli_select(prog, cols, x, SE)
            check(p == [1 - lE / E, lE / E])
            lE == E && check(p[1] == 0)   # every exposed host tracked: the lineage must move
            if lE < E
                w0 = apply_move!(cols, prog, 0, x1, SE; rng)
                check(close(w0, 1 - lI // (I + 1)) && ell(cols) == (lE, lI))
            end
            if lE ≥ 1
                w1 = apply_move!(cols, prog, 2, x1, SE; rng)
                check(close(w1, lE // (I + 1)) && ell(cols) == (lE - 1, lI + 1))
            end
            cols = seir_cols(lE, lI)
            check(kli_select(rec, cols, x, SE) == [I > lI ? 1 - lI / I : 0.0])
        end
        @info "SEIR per-outcome checks: $total checks, $bad failures ($nboundary states with every infectious host tracked)"
        @test total > 1_000
        @test nboundary > 20
        @test bad == 0
    end

    @testset "kli_decay equals compiled_decay" begin
        rng = MersenneTwister(3)
        worst = 0.0
        for _ in 1:300
            Ic, Ih = rand(rng, 0:10), rand(rng, 0:10)
            lc, lh = rand(rng, 0:Ic), rand(rng, 0:Ih)
            x = (S_c = rand(rng, 0:40), I_c = Ic, S_h = rand(rng, 0:40), I_h = Ih)
            θ = (β_cc = 3rand(rng), β_ch = 2rand(rng), β_hc = 2rand(rng), β_hh = 3rand(rng),
                 γ_c = 2rand(rng), γ_h = 2rand(rng), χ_c = rand(rng), χ_h = rand(rng),
                 B_c = rand(rng), B_h = rand(rng), N_c = 50.0, N_h = 60.0)
            cols = mers_cols(lc, lh)
            slots = kli_slots(MERS)
            al, pv = zeros(length(slots)), zeros(length(slots))
            d = kli_rates!(al, pv, slots, cols, x, θ, MERS)
            worst = max(worst, abs(d - compiled_decay(MERS, x, θ, [lc, lh], [Ic, Ih])))

            E, I = rand(rng, 0:10), rand(rng, 0:10)
            lE, lI = rand(rng, 0:E), rand(rng, 0:I)
            x = (S = rand(rng, 0:60), E = E, I = I, R = rand(rng, 0:10))
            θ = (β = 4rand(rng), σ = 2rand(rng), γ = 2rand(rng), ω = rand(rng), ψ = rand(rng), χ = rand(rng), N = 100.0)
            cols = seir_cols(lE, lI)
            slots = kli_slots(SE)
            al, pv = zeros(length(slots)), zeros(length(slots))
            d = kli_rates!(al, pv, slots, cols, x, θ, SE)
            worst = max(worst, abs(d - compiled_decay(SE, x, θ, [lE, lI], [E, I])))
        end
        @info "kli_decay vs compiled_decay: worst |Δ| = $worst"
        @test worst < 1e-12
    end

    @info h2("Isolated regular step, bit-exact against the compiled regular parts")
    @testset "regular_step! vs mers_compiled_regular_part!" begin
        rng = MersenneTwister(20261007)
        nmismatch = 0; nmoved = 0; worst = 0.0
        for _ in 1:400
            Beta_cc = 0.3 + 3rand(rng); Beta_hh = 0.3 + 3rand(rng)
            Beta_hc = 0.2 + 2rand(rng); Beta_ch = 0.2 + 2rand(rng)
            gamma_c = 0.2 + 2rand(rng); gamma_h = 0.2 + 2rand(rng)
            chi_c = 0.05 + 0.5rand(rng); chi_h = 0.05 + 0.5rand(rng)
            Bc = 0.3rand(rng); Bh = 0.3rand(rng)
            Nc = Float64(rand(rng, 30:200)); Nh = Float64(rand(rng, 30:200))
            Sc = rand(rng, 5:60); Ic = rand(rng, 1:15); Sh = rand(rng, 5:60); Ih = rand(rng, 1:15)
            lc = rand(rng, 0:Ic); lh = rand(rng, 0:Ih)
            dt = 0.3 + 2.7rand(rng); seed = rand(rng, 1:10^7)
            kw = (Beta_cc = Beta_cc, Beta_ch = Beta_ch, Beta_hc = Beta_hc, Beta_hh = Beta_hh,
                  gamma_c = gamma_c, gamma_h = gamma_h, chi_c = chi_c, chi_h = chi_h,
                  Bc = Bc, Bh = Bh, Nc = Nc, Nh = Nh)
            θ = (β_cc = Beta_cc, β_ch = Beta_ch, β_hc = Beta_hc, β_hh = Beta_hh, γ_c = gamma_c, γ_h = gamma_h,
                 χ_c = chi_c, χ_h = chi_h, B_c = Bc, B_h = Bh, N_c = Nc, N_h = Nh)
            c1, c2 = mers_cols(lc, lh), mers_cols(lc, lh)
            seed!(seed)
            l1, Sc1, Ic1, Sh1, Ih1 = PhyloPOMP.mers_compiled_regular_part!(
                c1, 0.0, 0.0, dt, Sc, Ic, Sh, Ih; kw..., model = MERS)
            u1 = rand(UInt64)
            seed!(seed)
            l2, x2, _ = regular_step!(c2, 0.0, 0.0, dt, (S_c = Sc, I_c = Ic, S_h = Sh, I_h = Ih), MERS, θ)
            u2 = rand(UInt64)
            ok = isapprox(l1, l2; atol = 1e-8, rtol = 1e-10) && (Sc1, Ic1, Sh1, Ih1) == Tuple(x2) &&
                 c1.cols == c2.cols && u1 == u2
            nmismatch += !ok
            worst = max(worst, abs(l1 - l2))
            nmoved += ell(c1) != (lc, lh)
        end
        @info "MERS isolation: 400 trials, $nmismatch mismatches, worst |Δll| = $worst, a lineage moved in $nmoved"
        @test nmismatch == 0
        @test nmoved ≥ 100
    end

    @testset "regular_step! vs compiled_regular_part!" begin
        rng = MersenneTwister(20261008)
        nmismatch = 0; nmoved = 0; worst = 0.0
        for _ in 1:400
            β = 0.5 + 4rand(rng); σ = 0.2 + 2rand(rng); γ = 0.2 + 2rand(rng); ω = 0.1 + 1.5rand(rng)
            ψ = rand(rng, Bool) ? 0.0 : 0.3rand(rng); χ = rand(rng, Bool) ? 0.0 : 0.3rand(rng)
            S = rand(rng, 5:80); E = rand(rng, 0:15); I = rand(rng, 1:15); R = rand(rng, 0:10)
            lE = rand(rng, 0:E); lI = rand(rng, 0:I)
            dt = 0.3 + 2.7rand(rng); seed = rand(rng, 1:10^7)
            θ = (β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, N = 100.0)
            c1, c2 = seir_cols(lE, lI), seir_cols(lE, lI)
            seed!(seed)
            l1, S1, E1, I1, R1 = PhyloPOMP.compiled_regular_part!(
                c1, 0.0, 0.0, dt, S, E, I, R; β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, pop = 100.0, model = SE)
            u1 = rand(UInt64)
            seed!(seed)
            l2, x2, _ = regular_step!(c2, 0.0, 0.0, dt, (S = S, E = E, I = I, R = R), SE, θ)
            u2 = rand(UInt64)
            ok = isapprox(l1, l2; atol = 1e-8, rtol = 1e-10) && (S1, E1, I1, R1) == Tuple(x2) &&
                 c1.cols == c2.cols && u1 == u2
            nmismatch += !ok
            worst = max(worst, abs(l1 - l2))
            nmoved += ell(c1) != (lE, lI)
        end
        @info "SEIR isolation: 400 trials, $nmismatch mismatches, worst |Δll| = $worst, a lineage moved in $nmoved"
        @test nmismatch == 0
        @test nmoved ≥ 100
    end

    @testset "explicit rng" begin
        θ = (β = 3.0, σ = 1.0, γ = 1.0, ω = 0.5, ψ = 0.2, χ = 0.1, N = 100.0)
        x = (S = 40, E = 5, I = 6, R = 3)
        run(rng) = (c = seir_cols(2, 3); (regular_step!(c, 0.0, 0.0, 2.0, x, SE, θ; rng)..., c.cols))
        @test run(MersenneTwister(9)) == run(MersenneTwister(9))
        # an explicit rng leaves the global RNG alone; default_rng() is the global RNG
        seed!(4); u0 = rand(UInt64)
        seed!(4); run(MersenneTwister(9)); @test rand(UInt64) == u0
        seed!(4); a = (c = seir_cols(2, 3); (regular_step!(c, 0.0, 0.0, 2.0, x, SE, θ)..., c.cols))
        seed!(4); @test run(default_rng()) == a
    end

    @info h2("Seed-matched Np=100 log-likelihoods: generic regular step inside the compiled filters")
    @testset "end-to-end (MERS)" begin
        function build(θ, x0, target_n, tmax, gseed)
            rng = MersenneTwister(gseed)
            for _ in 1:800
                g = simulate(MERS, θ; x0 = x0, graft = [1, 0], tmax = tmax, rng = rng,
                             demeset = NaiveMERS.Demes, samplemap = [Camel, Human])
                nsample(g) == target_n && return g
            end
            nothing
        end
        rm = MersenneTwister(16)
        nfinite = 0; nmatch = 0; nmatch_all = 0; ncombo = 0; worst = 0.0
        for combo in 1:6
            pc = rand(rm, (30, 50)); ph = rand(rm, (30, 50))
            th = (β_cc = 1 + 3rand(rm), β_ch = 0.3 + 2rand(rm), β_hc = 0.3 + 2rand(rm), β_hh = 1 + 3rand(rm),
                  γ_c = 0.5 + 1.5rand(rm), γ_h = 0.5 + 1.5rand(rm), χ_c = 0.2 + 0.4rand(rm), χ_h = 0.2 + 0.4rand(rm),
                  B_c = iseven(combo) ? 0.0 : 0.5, B_h = iseven(combo) ? 0.0 : 0.3,
                  N_c = Float64(pc), N_h = Float64(ph))
            x0 = (S_c = pc - 1, I_c = 1, S_h = ph, I_h = 0)
            g = build(th, x0, rand(rm, 3:8), 3.0 + 6rand(rm), hash((:regular, combo)))
            g === nothing && continue
            ncombo += 1
            kw = (Beta_cc = th.β_cc, Beta_ch = th.β_ch, Beta_hc = th.β_hc, Beta_hh = th.β_hh,
                  gamma_c = th.γ_c, gamma_h = th.γ_h, chi_c = th.χ_c, chi_h = th.χ_h,
                  Bc = th.B_c, Bh = th.B_h, Sc0 = x0.S_c / pc, Sh0 = x0.S_h / ph,
                  Ic0 = x0.I_c / pc, Ih0 = x0.I_h / ph, Nc = pc, Nh = ph)
            p1 = PhyloPOMP.mers_compiled_filter_pomp(g; kw...)
            p2 = PhyloPOMP.mers_compiled_filter_pomp(g; kw..., generic_regular = true)
            p3 = PhyloPOMP.mers_compiled_filter_pomp(g; kw..., generic_regular = true, generic_singular = true)
            for s in 1:4
                seed!(s); a = logLik(pfilter(p1, Np = 100))
                seed!(s); b = logLik(pfilter(p2, Np = 100))
                seed!(s); c = logLik(pfilter(p3, Np = 100))
                if isfinite(a) || isfinite(b) || isfinite(c)
                    nfinite += 1
                    nmatch += isapprox(a, b; atol = 1e-9, rtol = 1e-10)
                    nmatch_all += isapprox(a, c; atol = 1e-9, rtol = 1e-10)
                    worst = max(worst, abs(a - b), abs(a - c))
                end
            end
        end
        @info "MERS end-to-end: $ncombo genealogies, $nfinite finite comparisons, $nmatch match (generic regular), " *
              "$nmatch_all match (generic regular and singular), worst |Δ| = $worst"
        @test ncombo ≥ 4
        @test nfinite ≥ 10
        @test nmatch == nfinite
        @test nmatch_all == nfinite
    end

    @testset "end-to-end (SEIR)" begin
        function build(θ, x0, target_n, tmax, gseed)
            rng = MersenneTwister(gseed)
            for _ in 1:800
                g = simulate(SE, θ; x0 = x0, graft = [0, 1], tmax = tmax, rng = rng)
                nsample(g) == target_n && return g
            end
            nothing
        end
        rm = MersenneTwister(17)
        nfinite = 0; nmatch = 0; nmatch_all = 0; ncombo = 0; worst = 0.0
        for combo in 1:8
            β = 1 + 6rand(rm); σ = 0.3 + 3rand(rm); γ = 0.3 + 3rand(rm); ω = 0.1 + 2rand(rm); ψ = 0.05 + 0.4rand(rm)
            χ = combo % 3 == 0 ? 0.0 : 0.05 + 0.3rand(rm)
            pop = rand(rm, (60, 100)); x0 = (S = pop - 1, E = 0, I = 1, R = 0)
            th = (β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, N = Float64(pop))
            g = build(th, x0, rand(rm, 2:7), 2.0 + 4rand(rm), hash((:seir_regular, combo)))
            g === nothing && continue
            ncombo += 1
            kw = (β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, pop = pop,
                  S0 = x0.S / pop, E0 = 0.0, I0 = 1 / pop, R0 = 0.0)
            p1 = PhyloPOMP.compiled_filter_pomp(g; kw...)
            p2 = PhyloPOMP.compiled_filter_pomp(g; kw..., generic_regular = true)
            p3 = PhyloPOMP.compiled_filter_pomp(g; kw..., generic_regular = true, generic_singular = true)
            for s in 1:4
                seed!(s); a = logLik(pfilter(p1, Np = 100))
                seed!(s); b = logLik(pfilter(p2, Np = 100))
                seed!(s); c = logLik(pfilter(p3, Np = 100))
                if isfinite(a) || isfinite(b) || isfinite(c)
                    nfinite += 1
                    nmatch += isapprox(a, b; atol = 1e-9, rtol = 1e-10)
                    nmatch_all += isapprox(a, c; atol = 1e-9, rtol = 1e-10)
                    worst = max(worst, abs(a - b), abs(a - c))
                end
            end
        end
        @info "SEIR end-to-end: $ncombo genealogies, $nfinite finite comparisons, $nmatch match (generic regular), " *
              "$nmatch_all match (generic regular and singular), worst |Δ| = $worst"
        @test ncombo ≥ 5
        @test nfinite ≥ 20
        @test nmatch == nfinite
        @test nmatch_all == nfinite
    end
    @testset "deme count of the coloring must match the model" begin
        cols1 = Coloring(PhyloPOMP.Unstructured)                  # 1 deme; MERS and SEIR have 2
        x = (S_c = 30, I_c = 5, S_h = 30, I_h = 4)
        θ = (β_cc = 1.0, β_ch = 1.0, β_hc = 1.0, β_hh = 1.0, γ_c = 1.0, γ_h = 1.0, χ_c = 0.5, χ_h = 0.5,
             B_c = 0.0, B_h = 0.0, N_c = 60.0, N_h = 60.0)
        ev = find_event(MERS, :transmission_hc)
        @test_throws ArgumentError PhyloPOMP.singular_update!(cols1, Dict(), 1, x, θ, MERS)
        @test_throws ArgumentError regular_step!(cols1, 0.0, 0.0, 1.0, x, MERS, θ)
        @test_throws ArgumentError kli_select(ev, cols1, x, MERS)
        @test_throws ArgumentError apply_move!(cols1, ev, 2, x, MERS)
        msg = try; kli_select(ev, cols1, x, MERS); "" ; catch e; sprint(showerror, e); end
        @test occursin("1 demes", msg) && occursin("MERS", msg) && occursin("I_c", msg)
        cols2 = mers_cols(1, 1)                                   # the right count passes the check
        @test PhyloPOMP._check_demes(cols2, MERS) === nothing
    end
end

end # module KliRegularStepTest
