"""
Randomized property tests over full_transitions, reduce_event_indicator and
kli_binomial_ratio for every BIRTH/MIGRATION event of SEIR and MERS.
Hundreds of random (ell, n) states per property.
"""
module KliPropertiesTest

import ..Main: h1, h2

@info h1("Property-based verification over SEIR/MERS")

using Test
using Random: MersenneTwister
using PhyloPOMP
using PhyloPOMP: Event, EventType, BIRTH, MIGRATION,
    production_slots, enumerate_saturations, kli_binomial_ratio,
    full_transitions, reduce_event_indicator, reduced_transitions,
    IdentityTransition, InlineSameDemeTransition, CrossDemeTransition,
    ForkTransition, KLITransition, ReducedTransition

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

# All BIRTH/MIGRATION events, collected by type so new events are covered.
birth_migration_events() = vcat(
    [(:SEIR, e) for e in PhyloPOMP.SEIR.events if e.type in (BIRTH, MIGRATION)],
    [(:MERS, e) for e in PhyloPOMP.MERS.events if e.type in (BIRTH, MIGRATION)],
)

# Random state with ℓ_d <= n_d by construction.
function random_valid_state(rng, D; nmax = 12)
    n = [rand(rng, 1:nmax) for _ in 1:D]
    ℓ = [rand(rng, 0:n[d]) for d in 1:D]
    return ℓ, n
end

@testset verbose=true "Property-based verification" begin

    @info h2("Property 1: population/genealogy support ℓ_d <= n_d")
    @testset "P1: population support ell <= n" begin
        rng = MersenneTwister(20260819_1)
        events = birth_migration_events()
        ntrials = 400
        for _ in 1:ntrials
            _, event = events[rand(rng, 1:length(events))]
            D = length(production_slots(event))
            ℓ, n = random_valid_state(rng, D)

            # Generator check.
            @test all(ℓ[d] <= n[d] for d in 1:D)

            ts = full_transitions(event, ℓ, n)
            # enumerate_saturations bounds s_d by min(r_d, ℓ_d).
            S = enumerate_saturations(event, ℓ)
            @test all(all(s[d] <= ℓ[d] for d in 1:D) for s in S)
            @test all(all(t.s[d] <= ℓ[d] for d in 1:D) for t in ts)

            # ForkTransition/CrossDemeTransition/InlineSameDemeTransition
            # slot demes never imply an occupied count exceeding ℓ in their
            # own deme either -- spot check via sum(s) <= sum(ℓ).
            @test all(sum(t.s) <= sum(ℓ) for t in ts)

            rts = reduce_event_indicator(ts)
            # A reduced Fork's slot-deme multiset never exceeds ℓ per deme
            # (mirrors the full-level check one level up, after grouping).
            for rt in rts
                if rt.kind == :fork
                    _, _, slot_demes = rt.key
                    counts = Dict{Int,Int}()
                    for d in slot_demes
                        counts[d] = get(counts, d, 0) + 1
                    end
                    for (d, c) in counts
                        @test c <= ℓ[d]
                    end
                end
            end
        end
    end

    @info h2("Property 2: full-to-reduced conservation, Σφ_u == ΣΦ_u")
    @testset "P2: full-to-reduced conservation" begin
        rng = MersenneTwister(20260819_2)
        events = birth_migration_events()
        ntrials = 500
        for _ in 1:ntrials
            _, event = events[rand(rng, 1:length(events))]
            D = length(production_slots(event))
            ℓ, n = random_valid_state(rng, D)

            ts = full_transitions(event, ℓ, n)
            rts = reduce_event_indicator(ts)

            # Both sides sum the same phi_u values, partitioned differently.
            @test sum(t.phi for t in ts) == sum(rt.Φ for rt in rts)

            # Partition property: every full transition appears in exactly
            # one reduced group, none dropped/duplicated.
            ngrouped = sum(length(rt.transitions) for rt in rts)
            @test ngrouped == length(ts)

            # convenience composition agrees with the two-step call.
            rts2 = reduced_transitions(event, ℓ, n)
            @test sum(rt.Φ for rt in rts2) == sum(t.phi for t in ts)
        end
    end

    @info h2("Property 3: proposal support, Φ_u(z) >= 0 always")
    @testset "P3: proposal support Phi_u >= 0" begin
        rng = MersenneTwister(20260819_3)
        events = birth_migration_events()
        ntrials = 500
        for _ in 1:ntrials
            _, event = events[rand(rng, 1:length(events))]
            D = length(production_slots(event))
            ℓ, n = random_valid_state(rng, D; nmax = 20)

            ts = full_transitions(event, ℓ, n)
            # phi_u is a ratio of non-negative binomials.
            @test all(t.phi >= 0 for t in ts)

            rts = reduce_event_indicator(ts)
            # Φ_u is non-negative, so proposal weights are well defined.
            @test all(rt.Φ >= 0 for rt in rts)

            # Check kli_binomial_ratio directly, bypassing classification.
            S = enumerate_saturations(event, ℓ)
            for s in S
                φ = kli_binomial_ratio(event, s, ℓ, n)
                @test φ >= 0
            end
        end
    end

    @info h2("Property 4: impossible transitions have zero compatibility")
    @testset "P4: impossible transitions -> zero, not error/crash" begin
        rng = MersenneTwister(20260819_4)
        events = birth_migration_events()
        ntrials = 400
        for _ in 1:ntrials
            _, event = events[rand(rng, 1:length(events))]
            r = production_slots(event)
            D = length(r)

            case = rand(rng, 1:3)
            if case == 1
                # n_d < r_d for some deme: C(n_d, r_d) = 0, so the WHOLE
                # product must be exactly 0 for every s, not an error.
                n = [rand(rng, 0:1) for _ in 1:D]      # tiny n
                ℓ = [rand(rng, 0:n[d]) for d in 1:D]
                # (only meaningful if some r_d > n_d actually occurs)
                if any(r[d] > n[d] for d in 1:D)
                    S = enumerate_saturations(r, ℓ)
                    for s in S
                        φ = kli_binomial_ratio(event, s, ℓ, n)
                        @test φ == 0
                    end
                    # full_transitions must not throw -- it degrades to an
                    # all-zero-phi vector, not a crash.
                    ts = full_transitions(event, ℓ, n)
                    @test all(t.phi == 0 for t in ts)
                end
            elseif case == 2
                # ell_d > n_d for some deme -- should never be a valid MODEL
                # state, but the machinery must degrade to zero, not crash or
                # return a nonsensical (e.g. negative or >1) answer.
                n = [rand(rng, 1:6) for _ in 1:D]
                ℓ = [n[d] + rand(rng, 1:5) for d in 1:D]   # ell > n, every deme
                # enumerate_saturations ignores n, so it does not error here.
                S = enumerate_saturations(r, ℓ)
                @test S isa Vector
                for s in S
                    φ = kli_binomial_ratio(event, s, ℓ, n)
                    @test φ == 0   # n_d - ell_d < 0 zeroes every factor
                end
                ts = full_transitions(event, ℓ, n)
                @test all(t.phi == 0 for t in ts)
                rts = reduce_event_indicator(ts)
                @test all(rt.Φ == 0 for rt in rts)
            else
                # s_d < 0 is never enumerated; kli_binomial_ratio must still return 0.
                n = [rand(rng, 3:10) for _ in 1:D]
                ℓ = [rand(rng, 0:n[d]) for d in 1:D]
                s = [rand(rng, Bool) ? -rand(rng, 1:3) : rand(rng, 0:min(r[d], ℓ[d])) for d in 1:D]
                if any(sd < 0 for sd in s)
                    φ = kli_binomial_ratio(event, s, ℓ, n)
                    @test φ == 0
                end
                # s_d > r_d (structurally infeasible: requesting more
                # production than the event has slots for) -- caught by
                # safe_binomial's b<0 rule (r_d - s_d < 0).
                s2 = [r[d] + rand(rng, 1:4) for d in 1:D]
                φ2 = kli_binomial_ratio(event, s2, ℓ, n)
                @test φ2 == 0
            end
        end
    end

    @info h2("Property 5: genealogy-operation structural invariants")
    @testset "P5: KLITransition structural invariants (fuzz)" begin
        rng = MersenneTwister(20260819_5)
        events = birth_migration_events()
        ntrials = 600
        for _ in 1:ntrials
            _, event = events[rand(rng, 1:length(events))]
            D = length(production_slots(event))
            ℓ, n = random_valid_state(rng, D; nmax = 10)

            ts = full_transitions(event, ℓ, n)
            for t in ts
                if t isa IdentityTransition
                    @test sum(t.s) == 0
                elseif t isa InlineSameDemeTransition
                    @test sum(t.s) == 1
                    # the occupied slot's deme equals the ancestral deme
                    # (event.from), by construction/classification rule.
                    @test t.deme == event.from
                elseif t isa CrossDemeTransition
                    @test sum(t.s) == 1
                    @test t.ancestral_deme != t.target_deme
                    @test t.ancestral_deme == event.from
                elseif t isa ForkTransition
                    @test sum(t.s) >= 2
                    @test length(t.slot_demes) == sum(t.s)
                    @test t.ancestral_deme == event.from
                else
                    @test false  # unreachable: exhaustive over concrete subtypes
                end
                # phi is always the same value production would recompute
                # directly -- cross-check the two code paths agree.
                @test t.phi == kli_binomial_ratio(event, t.s, ℓ, n)
            end

            # Reduced-level structural invariants: every ReducedTransition's
            # `kind` matches its `key`'s first element, and noop/cross/fork
            # keys never collide across kinds.
            rts = reduce_event_indicator(ts)
            for rt in rts
                @test rt.kind == rt.key[1]
                if rt.kind == :cross
                    _, d_from, d_to = rt.key
                    @test d_from != d_to
                elseif rt.kind == :fork
                    _, d_anc, slot_demes = rt.key
                    @test length(slot_demes) >= 2
                end
            end
            keys = [rt.key for rt in rts]
            @test length(keys) == length(unique(keys))  # no duplicate groups
        end
    end

end

end # module
