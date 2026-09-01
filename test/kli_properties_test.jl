"""
M12 acceptance-gate tests: generic, RANDOMIZED property tests over the KLI
compiler machinery built in M02-M09 (`src/examples/mgp_phi.jl`,
`mgp_transitions.jl`, `mgp_reduce.jl`), exercised against the two models
that actually exist in this repository (`PhyloPOMP.SEIR`, `PhyloPOMP.MERS`).

Unlike `kli_phi_test.jl`/`kli_full_transitions_test.jl`/`kli_reduce_test.jl`
(M02-M04's own acceptance tests, which check a handful of hand-picked
`(ℓ,n)` instances term-by-term against `mers_filter_suite.tex`), this file
implements broad RANDOMIZED-TRIAL sweeps (hundreds of trials per property,
varying `(ℓ,n,s)` programmatically) of five structural invariants the
compiler must satisfy for EVERY state, not just the hand-picked ones. This
generalizes the pattern M08/M09 already used for their own 400-trial
isolated-unit checks and 289-2903-trial Gate-5 finite-comparison sweeps.

No QuickCheck-style property-testing library exists in this repo's
`Project.toml`/`Manifest.toml` (checked directly -- neither `test/Project.toml`
nor the root `Project.toml` list anything resembling one); per this
milestone's own instruction, randomized-trial loops are hand-written here
using Julia's standard `Test`/`Random` machinery, matching the pattern
`kli_seir_compiled_test.jl`/`kli_mers_compiled_test.jl` already established
(`MersenneTwister` seeded once per testset, `rand(rng, ...)` draws).
"""
module KliPropertiesTest

import ..Main: h1, h2

@info h1("Property-based verification over SEIR/MERS (M12)")

using Test
using Random: MersenneTwister
using PhyloPOMP
using PhyloPOMP: Event, EventType, BIRTH, MIGRATION,
    production_slots, enumerate_saturations, kli_binomial_ratio,
    full_transitions, reduce_event_indicator, reduced_transitions,
    IdentityTransition, InlineSameDemeTransition, CrossDemeTransition,
    ForkTransition, KLITransition, ReducedTransition

find_event(model, name) = model.events[findfirst(e -> e.name == name, model.events)]

# Every BIRTH/MIGRATION event in both models that actually exist in this
# repository -- SEIR's `infection` (BIRTH)/`progression` (MIGRATION) and
# MERS's `transmission_cc`/`transmission_hh`/`transmission_hc`/`transmission_ch`
# (all BIRTH, via `move=fork(...)`). Collected generically (filtered by
# `event.type`), not hand-listed by name, so this automatically covers any
# future BIRTH/MIGRATION event added to either model without editing this
# file.
birth_migration_events() = vcat(
    [(:SEIR, e) for e in PhyloPOMP.SEIR.events if e.type in (BIRTH, MIGRATION)],
    [(:MERS, e) for e in PhyloPOMP.MERS.events if e.type in (BIRTH, MIGRATION)],
)

# Random valid (ℓ, n) state for an event with D lineage-carrying demes:
# n_d drawn first (bounded, but comfortably larger than any r_d in SEIR/MERS,
# both of which have r_d <= 2), then ℓ_d drawn uniformly in [0, n_d] so the
# MODEL invariant ℓ_d <= n_d holds by construction (this is itself asserted
# below, as a sanity check on the generator, not just assumed).
function random_valid_state(rng, D; nmax = 12)
    n = [rand(rng, 1:nmax) for _ in 1:D]
    ℓ = [rand(rng, 0:n[d]) for d in 1:D]
    return ℓ, n
end

@testset verbose=true "Property-based verification (M12)" begin

    # =========================================================================
    @info h2("Property 1: population/genealogy support ℓ_d <= n_d")
    @testset "P1: population support ell <= n" begin
        rng = MersenneTwister(20260819_1)
        events = birth_migration_events()
        ntrials = 400
        for _ in 1:ntrials
            _, event = events[rand(rng, 1:length(events))]
            D = length(production_slots(event))
            ℓ, n = random_valid_state(rng, D)

            # Sanity check on our OWN generator, not the compiler: every
            # trial actually satisfies the model invariant it's supposed to.
            @test all(ℓ[d] <= n[d] for d in 1:D)

            ts = full_transitions(event, ℓ, n)
            # No saturation this milestone's machinery ever considers can
            # imply s_d > ℓ_d (enumerate_saturations bounds s_d by
            # min(r_d, ℓ_d) structurally) -- confirmed on the actual
            # enumerated saturations, not just the construction rule.
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

    # =========================================================================
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

            # Exact Rational{Int} equality -- both sides are sums of the same
            # underlying phi_u values, just partitioned differently
            # (M04's `reduce_event_indicator` docstring invariant, checked
            # here broadly rather than at a handful of hand-picked instances).
            @test sum(t.phi for t in ts) == sum(rt.Φ for rt in rts)

            # Partition property: every full transition appears in exactly
            # one reduced group's provenance, none dropped/duplicated.
            total_provenance = sum(length(rt.transitions) for rt in rts)
            @test total_provenance == length(ts)

            # convenience composition agrees with the two-step call.
            rts2 = reduced_transitions(event, ℓ, n)
            @test sum(rt.Φ for rt in rts2) == sum(t.phi for t in ts)
        end
    end

    # =========================================================================
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
            # phi_u itself is never negative -- a product of
            # safe_binomial(...)/safe_binomial(...) ratios, each a ratio of
            # two non-negative integers by construction (safe_binomial
            # returns 0 or a genuine non-negative binomial coefficient).
            @test all(t.phi >= 0 for t in ts)

            rts = reduce_event_indicator(ts)
            # Φ_u is a sum of non-negative phi_u values within each group,
            # so it is never negative either -- the quantity the compiled
            # SEIR/MERS proposals (mgp_seir_filter.jl/mgp_mers_filter.jl)
            # actually assign probability mass proportional to (directly, or
            # via the collapsed :noop group's weighted variant, M09 Finding
            # 2) must never be negative for the proposal to be well-defined.
            @test all(rt.Φ >= 0 for rt in rts)

            # Every individual saturation's phi_u is also checked directly
            # via kli_binomial_ratio (bypassing full_transitions'
            # classification), to isolate the two code paths.
            S = enumerate_saturations(event, ℓ)
            for s in S
                φ = kli_binomial_ratio(event, s, ℓ, n)
                @test φ >= 0
            end
        end
    end

    # =========================================================================
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
                n = [rand(rng, 0:1) for _ in 1:D]      # deliberately tiny n
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
                # enumerate_saturations does not consult n at all (it only
                # bounds by min(r_d, ell_d)), so it does not error either --
                # confirmed explicitly, not assumed.
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
                # s_d < 0 fed directly to kli_binomial_ratio (never produced
                # by enumerate_saturations, but a defensive-convention check
                # per mgp_phi.jl's documented, tested design: M02 found and
                # fixed a real bug here).
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

    # =========================================================================
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
