"""
Compiled MERS filter vs NaiveMERS: seed-matched Np=1 log-likelihoods must
agree on simulated genealogies. naive_oracle_filter_pomp is
NaiveMERS.filter_pomp taking gen as an argument. Genealogies use
NaiveMERS.Demes so the oracle's deme comparisons see the same enum type.
"""
module KliMersCompiledTest

import ..Main: h1, h2

@info h1("Compiled MERS filter vs. mers_naive.jl")

using Test
using PhyloPOMP
using PhyloPOMP.NaiveMERS
using PhyloPOMP: Name, Prob, Time, Genealogy, Coloring, ell,
    full_transitions, IdentityTransition, InlineSameDemeTransition,
    CrossDemeTransition, reduce_event_indicator
using Random: MersenneTwister, seed!
import PartiallyObservedMarkovProcesses as POMP

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

# Saturations are absent when ell is too small; phi_of returns 0 for them.
phi_of(ts, T) = begin
    i = findfirst(t -> t isa T, ts)
    isnothing(i) ? zero(Rational{Int}) : ts[i].phi
end

# NaiveMERS.filter_pomp with gen as an argument.
naive_oracle_filter_pomp(
    gen::Genealogy;
    Beta_cc = 4.0, Beta_ch = 0.0, Beta_hc = 1.0, Beta_hh = 4.0,
    gamma_c = 1.0, gamma_h = 1.0,
    chi_c = 1.0, chi_h = 0.0,
    Bc = 0.1, Bh = 0.03,
    Sc0 = 1.0, Sh0 = 1.0, Ic0 = 0.01, Ih0 = 0.0,
    Nc = 10000, Nh = 10000,
) = begin
    POMP.pomp(
        params = (
            Beta_cc = Float64(Beta_cc), Beta_ch = Float64(Beta_ch),
            Beta_hc = Float64(Beta_hc), Beta_hh = Float64(Beta_hh),
            gamma_c = Float64(gamma_c), gamma_h = Float64(gamma_h),
            chi_c = Float64(chi_c), chi_h = Float64(chi_h),
            Bc = Float64(Bc), Bh = Float64(Bh),
            Sc0 = Float64(Sc0), Sh0 = Float64(Sh0),
            Ic0 = Float64(Ic0), Ih0 = Float64(Ih0),
            Nc = Float64(Nc), Nh = Float64(Nh),
        ),
        t0 = timezero(gen),
        times = times(gen),
        rinit = function (; Sc0, Sh0, Ic0, Ih0, Nc, Nh, _...)
            fc = Nc / (Sc0 + Ic0)
            fh = Nh / (Sh0 + Ih0)
            (
                node = one(Name),
                ll = zero(Float64),
                cols = Coloring(NaiveMERS.Demes),
                Sc = round(Int64, fc * Sc0),
                Ic = round(Int64, fc * Ic0),
                Sh = round(Int64, fh * Sh0),
                Ih = round(Int64, fh * Ih0),
            )
        end,
        rprocess = POMP.onestep(
            function (; node, ll, cols, geneal, Sc, Ic, Sh, Ih, t, dt, args...)
                cols = copy(cols)
                ll = zero(Float64)
                ll, Sc, Ic, Sh, Ih = NaiveMERS.singular_part!(
                    cols, ll, geneal, node, Sc, Ic, Sh, Ih; args...,
                )
                if isfinite(ll)
                    ll, Sc, Ic, Sh, Ih = NaiveMERS.regular_part!(
                        cols, ll, t, dt, Sc, Ic, Sh, Ih; args...,
                    )
                end
                (; node = node + one(Name), ll, cols, Sc, Ic, Sh, Ih)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,),
    )
end

function build_genealogy(θ, x0, target_n, tmax, gseed; attempts = 800)
    rng = MersenneTwister(gseed)
    g = nothing
    for _ in 1:attempts
        g = simulate(
            PhyloPOMP.MERS, θ; x0 = x0, graft = [1, 0], tmax = tmax, rng = rng,
            demeset = NaiveMERS.Demes, samplemap = [NaiveMERS.Camel, NaiveMERS.Human],
        )
        if nsample(g) == target_n
            return g, true
        end
    end
    g, false
end

function compare_ll(g, θ, pop_c, pop_h, x0; nseeds)
    p_naive = naive_oracle_filter_pomp(
        g;
        Beta_cc = θ.β_cc, Beta_ch = θ.β_ch, Beta_hc = θ.β_hc, Beta_hh = θ.β_hh,
        gamma_c = θ.γ_c, gamma_h = θ.γ_h, chi_c = θ.χ_c, chi_h = θ.χ_h,
        Bc = θ.B_c, Bh = θ.B_h,
        Sc0 = x0.S_c / pop_c, Sh0 = x0.S_h / pop_h,
        Ic0 = x0.I_c / pop_c, Ih0 = x0.I_h / pop_h,
        Nc = pop_c, Nh = pop_h,
    )
    p_comp = PhyloPOMP.mers_compiled_filter_pomp(
        g;
        Beta_cc = θ.β_cc, Beta_ch = θ.β_ch, Beta_hc = θ.β_hc, Beta_hh = θ.β_hh,
        gamma_c = θ.γ_c, gamma_h = θ.γ_h, chi_c = θ.χ_c, chi_h = θ.χ_h,
        Bc = θ.B_c, Bh = θ.B_h,
        Sc0 = x0.S_c / pop_c, Sh0 = x0.S_h / pop_h,
        Ic0 = x0.I_c / pop_c, Ih0 = x0.I_h / pop_h,
        Nc = pop_c, Nh = pop_h,
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

@testset verbose=true "Compiled MERS filter" begin

    @info h2("TCC/THH: naive's drawless no-fork term equals " *
             "Phi_id + ell*Phi_inline, not the plain noop sum")
    @testset "TCC/THH weighted aggregate" begin
        tcc = find_event(PhyloPOMP.MERS, :transmission_cc)
        thh = find_event(PhyloPOMP.MERS, :transmission_hh)

        # camel side: ell_C=2, I_C=5 (post-event)
        ts = full_transitions(tcc, [2, 0], [5, 3])
        Φid = only(filter(t -> t isa IdentityTransition, ts)).phi
        Φinl = only(filter(t -> t isa InlineSameDemeTransition, ts)).phi
        @test Φid == 3 // 10 && Φinl == 3 // 10
        weighted = Φid + 2 * Φinl
        @test weighted == 9 // 10
        naive_raw = 1 - (2 * 1) // (5 * 4)   # naive's 1-ell*(ell-1)/(I*(I-1))
        @test weighted == naive_raw
        # This does NOT equal reduce_event_indicator's plain (unweighted) sum.
        rts = reduce_event_indicator(ts)
        Φnoop = only(filter(r -> r.kind == :noop, rts)).Φ
        @test Φnoop == 3 // 5
        @test Φnoop != weighted

        # Other (ell, I) cases; ell_C=1 has no Fork saturation, so phi_of returns 0.
        # Identity + ell*Inline + C(ell,2)*Fork sums to 1 in each.
        for (ellc, Ic) in ((3, 5), (1, 4), (4, 7))
            ts2 = full_transitions(tcc, [ellc, 0], [Ic, 3])
            Φid2 = phi_of(ts2, IdentityTransition)
            Φinl2 = phi_of(ts2, InlineSameDemeTransition)
            Φfork2 = phi_of(ts2, PhyloPOMP.ForkTransition)
            @test Φid2 + ellc * Φinl2 + binomial(ellc, 2) * Φfork2 == 1
            naive2 = 1 - (ellc * (ellc - 1)) // (Ic * (Ic - 1))
            @test Φid2 + ellc * Φinl2 == naive2
        end

        # THH mirrors TCC exactly (I_H=6, ell_H=2 -- the length(ts)==3
        # instance already hand-verified in kli_full_transitions_test.jl).
        ts3 = full_transitions(thh, [0, 2], [3, 6])
        Φidh = only(filter(t -> t isa IdentityTransition, ts3)).phi
        Φinlh = only(filter(t -> t isa InlineSameDemeTransition, ts3)).phi
        weighted_h = Φidh + 2 * Φinlh
        naive_h = 1 - (2 * 1) // (6 * 5)
        @test weighted_h == naive_h
    end

    @info h2("THC/TCH no-move/cross check: no move = Identity + " *
             "ell*InlineSameDeme, cross = boost(Phi,1/ell)")
    @testset "THC/TCH identity+cross" begin
        thc = find_event(PhyloPOMP.MERS, :transmission_hc)
        ellc, ellh, Ic, Ih = 2, 1, 5, 4
        ts = full_transitions(thc, [ellc, ellh], [Ic, Ih])
        Φid = only(filter(t -> t isa IdentityTransition, ts)).phi
        Φcr = only(filter(t -> t isa CrossDemeTransition, ts)).phi
        Φinl = only(filter(t -> t isa InlineSameDemeTransition, ts)).phi
        @test Φid == (1 - ellc // Ic) * (1 - ellh // Ih)
        # naive's k==3 target: no lineage passes to the new human, whether
        # the parent is untracked (Identity) or tracked (InlineSameDeme)
        @test Φid + ellc * Φinl == 1 - ellh // Ih
        # ...and it stays positive when every camel is tracked
        tsf = full_transitions(thc, [2, 1], [2, 4])
        @test sum(t.phi for t in tsf if t isa IdentityTransition; init=0//1) +
              2 * only(filter(t -> t isa InlineSameDemeTransition, tsf)).phi ==
              1 - 1 // 4
        @test Φcr * ellc == 3 // 10          # boost(Φcr,1/ellc) == Φcr*ellc
    end

    @info h2("Decay: this MERS instance reduces to 18/5")
    @testset "decay instance" begin
        γc, γh = 1 // 2, 2 // 5
        χc, χh = 1 // 10, 3 // 10
        x = (S_c = 0 // 1, I_c = 5 // 1, S_h = 0 // 1, I_h = 3 // 1)
        θ = (β_cc = 0 // 1, β_ch = 0 // 1, β_hc = 0 // 1, β_hh = 0 // 1,
             γ_c = γc, γ_h = γh, χ_c = χc, χ_h = χh,
             B_c = 0 // 1, B_h = 0 // 1, N_c = 1 // 1, N_h = 1 // 1)
        ℓ = [2, 3]   # ell_C=2 (above threshold), ell_H=3 (boundary)
        n = [5, 3]
        λ_tex = total_decay(PhyloPOMP.MERS, x, θ, ℓ, n)
        @test λ_tex == 13 // 5
        # compiled_decay computes in Float64, so compare with isapprox.
        total = compiled_decay(PhyloPOMP.MERS, x, θ, ℓ, n)
        @test isapprox(total, 18 / 5; atol = 1e-12)
    end

    @info h2("Isolated regular-part unit check: many (state, ell, params) " *
             "draws, mixed event realizations, bit-exact against " *
             "NaiveMERS.regular_part!")
    @testset "mers_compiled_regular_part! isolation" begin
        rng = MersenneTwister(20260918)
        nmismatch = 0
        ntried = 0
        for _ in 1:400
            Beta_cc = 0.3 + 3 * rand(rng); Beta_hh = 0.3 + 3 * rand(rng)
            Beta_hc = 0.2 + 2 * rand(rng); Beta_ch = 0.2 + 2 * rand(rng)
            gamma_c = 0.2 + 2 * rand(rng); gamma_h = 0.2 + 2 * rand(rng)
            chi_c = 0.05 + 0.5 * rand(rng); chi_h = 0.05 + 0.5 * rand(rng)
            Bc = 0.3 * rand(rng); Bh = 0.3 * rand(rng)
            Nc = Float64(rand(rng, 30:200)); Nh = Float64(rand(rng, 30:200))
            Sc = rand(rng, 5:60); Ic = rand(rng, 1:15)
            Sh = rand(rng, 5:60); Ih = rand(rng, 1:15)
            ellc = rand(rng, 0:Ic); ellh = rand(rng, 0:Ih)
            dt = 0.3 + 2.7 * rand(rng)
            seed = rand(rng, 1:10^7)

            cols_n = Coloring(NaiveMERS.Demes)
            cols_c = Coloring(NaiveMERS.Demes)
            lin = 1
            for _ in 1:ellc
                push!(cols_n[NaiveMERS.Camel], lin); push!(cols_c[NaiveMERS.Camel], lin); lin += 1
            end
            for _ in 1:ellh
                push!(cols_n[NaiveMERS.Human], lin); push!(cols_c[NaiveMERS.Human], lin); lin += 1
            end

            kwargs = (
                Beta_cc = Beta_cc, Beta_ch = Beta_ch, Beta_hc = Beta_hc, Beta_hh = Beta_hh,
                gamma_c = gamma_c, gamma_h = gamma_h, chi_c = chi_c, chi_h = chi_h,
                Bc = Bc, Bh = Bh, Nc = Nc, Nh = Nh,
            )

            seed!(seed)
            ll_n, Scn, Icn, Shn, Ihn = NaiveMERS.regular_part!(
                cols_n, 0.0, 0.0, dt, Sc, Ic, Sh, Ih; kwargs...,
            )
            seed!(seed)
            ll_c, Scc, Icc, Shc, Ihc = PhyloPOMP.mers_compiled_regular_part!(
                cols_c, 0.0, 0.0, dt, Sc, Ic, Sh, Ih; kwargs..., model = PhyloPOMP.MERS,
            )

            ntried += 1
            ok = isapprox(ll_n, ll_c; atol = 1e-8, rtol = 1e-10) &&
                 (Scn, Icn, Shn, Ihn) == (Scc, Icc, Shc, Ihc) &&
                 ell(cols_n) == ell(cols_c)
            ok || (nmismatch += 1)
        end
        @info "mers_compiled_regular_part! isolation: $ntried trials, $nmismatch mismatches"
        @test nmismatch == 0
    end

    @info h2("End-to-end log-likelihood over many (params, genealogy, seed) " *
             "combinations; beta_hc/beta_ch nonzero")
    @testset "end-to-end log-likelihood sweep" begin
        rng_master = MersenneTwister(20260919)
        total_finite = 0
        total_match = 0
        ncombos = 0
        worst_abs = 0.0
        for combo in 1:15
            β_cc = 1.0 + 3.0 * rand(rng_master)
            β_hh = 1.0 + 3.0 * rand(rng_master)
            β_hc = 0.3 + 2.0 * rand(rng_master)   # never zero
            β_ch = 0.3 + 2.0 * rand(rng_master)   # never zero
            γ_c = 0.5 + 1.5 * rand(rng_master)
            γ_h = 0.5 + 1.5 * rand(rng_master)
            χ_c = 0.2 + 0.4 * rand(rng_master)
            χ_h = 0.2 + 0.4 * rand(rng_master)
            B_c = 0.0
            B_h = 0.0
            pop_c = rand(rng_master, (20, 30, 40))
            pop_h = rand(rng_master, (20, 30, 40))
            x0 = (S_c = pop_c - 1, I_c = 1, S_h = pop_h, I_h = 0)
            θ = (β_cc = β_cc, β_ch = β_ch, β_hc = β_hc, β_hh = β_hh,
                 γ_c = γ_c, γ_h = γ_h, χ_c = χ_c, χ_h = χ_h,
                 B_c = B_c, B_h = B_h, N_c = Float64(pop_c), N_h = Float64(pop_h))
            target_n = rand(rng_master, 1:4)
            tmax = 3.0 + 5.0 * rand(rng_master)

            g, found = build_genealogy(θ, x0, target_n, tmax, hash((:m09, combo)))
            found || continue
            ncombos += 1

            nfinite, nmatch, worst = compare_ll(g, θ, pop_c, pop_h, x0; nseeds = 100)
            total_finite += nfinite
            total_match += nmatch
            worst_abs = max(worst_abs, worst)
            @test nfinite == nmatch
        end
        @info "End-to-end sweep: $ncombos parameter/genealogy combos, " *
              "$total_finite finite log-likelihood comparisons, " *
              "$total_match exact matches, worst |Δll|=$worst_abs"
        @test ncombos ≥ 8
        @test total_finite ≥ 20
    end

end

end # module KliMersCompiledTest
