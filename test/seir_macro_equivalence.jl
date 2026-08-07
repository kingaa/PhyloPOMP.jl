"""
Executable equivalence checks between the `@mgp`-generated `SEIR` event table
and the model encoded by `NaiveSEIR`.

The proved scope is the model layer: compartments, lineage demes, population
hazards, population increments, event classifications, regular driver rates,
and inter-event decay. The parameter correspondence is `N = pop`, and exact
equivalence requires `χ = 0` because the current macro table contains the
non-destructive `ψ` sampling mark but no destructive `χ` sampling mark.

This file does not claim full-filter equivalence: the KLI move and singular
weight functions in `mgp_filter.jl` are still explicit stubs.
"""
module SEIRMacroEquivalenceTest

import ..Main: h1, h2

@info h1("macro-generated SEIR model equivalence")

using Test
using Random: MersenneTwister, rand
using PhyloPOMP

const Existing = PhyloPOMP.NaiveSEIR
const MacroModel = PhyloPOMP.SEIR
const ExplicitTable = PhyloPOMP.SEIR_REFERENCE

event_signature(ev) = (
    ev.name,
    ev.Δ,
    ev.r,
    ev.type,
    ev.from,
    ev.into,
    ev.regular,
    ev.observed,
)

function coloring(ellE, ellI)
    cols = Coloring(Existing.Demes)
    for lineage in 1:ellE
        plant!(cols, Existing.Expos, lineage)
    end
    for lineage in 1:ellI
        plant!(cols, Existing.Infec, ellE + lineage)
    end
    cols
end

function existing_rates(x, θ, ellE, ellI; χ = 0.0)
    alpha = zeros(6)
    pi = zeros(6)
    decay = Existing.event_rates!(
        alpha, pi, coloring(ellE, ellI),
        x.S, x.E, x.I, x.R;
        β = θ.β, σ = θ.σ, γ = θ.γ, ω = θ.ω,
        ψ = θ.ψ, χ, pop = θ.N,
    )
    alpha .* pi, decay
end

function macro_rates(x, θ, ellE, ellI)
    hazards = [ev.hazard(x, θ) for ev in MacroModel.events]
    pIoff = x.I > 0 ? 1 - ellI/x.I : 0.0
    pIon = x.I > 0 ? ellI/x.I : 0.0
    pEoff = x.E > 0 ? 1 - ellE/x.E : 0.0
    pEon = x.E > 0 ? ellE/x.E : 0.0

    driver = [
        hazards[1]*pIoff,
        hazards[1]*pIon,
        hazards[2]*pEoff,
        hazards[2]*pEon,
        hazards[3]*pIoff,
        hazards[4],
    ]
    decay = hazards[5] + hazards[3]*(1-pIoff)
    driver, decay
end

@testset verbose=true "macro-generated SEIR model equivalence" begin

    @info h2("event-table structure")
    @test MacroModel.compartments == [:S, :E, :I, :R]
    @test MacroModel.demes == [:E, :I]
    @test map(event_signature, MacroModel.events) ==
        map(event_signature, ExplicitTable.events)

    expected = [
        (:infection,   [:S=>-1, :E=>1], [1, 1], PhyloPOMP.BIRTH,     2, [1],   true,  false),
        (:progression, [:E=>-1, :I=>1], [0, 1], PhyloPOMP.MIGRATION, 1, [2],   true,  false),
        (:recovery,    [:I=>-1, :R=>1], [0, 0], PhyloPOMP.DEATH,     2, Int[], true,  false),
        (:waning,      [:R=>-1, :S=>1], [0, 0], PhyloPOMP.NEUTRAL,   0, Int[], true,  false),
        (:sampling,    Pair{Symbol,Int}[], [0, 1], PhyloPOMP.SAMPLE, 2, Int[], false, true),
    ]
    @test map(event_signature, MacroModel.events) == expected

    x = (S = 90, E = 4, I = 5, R = 1)
    θ = (β = 4.0, σ = 0.8, γ = 0.5, ω = 0.2,
         ψ = 0.03, χ = 0.0, N = 100.0)
    @test [ev.hazard(x, θ) for ev in MacroModel.events] ≈
        [θ.β*x.S*x.I/θ.N, θ.σ*x.E, θ.γ*x.I, θ.ω*x.R, θ.ψ*x.I]

    expected_states = [
        (S = 89, E = 5, I = 5, R = 1),
        (S = 90, E = 3, I = 6, R = 1),
        (S = 90, E = 4, I = 4, R = 2),
        (S = 91, E = 4, I = 5, R = 0),
        x,
    ]
    @test [PhyloPOMP.apply_pop(x, ev) for ev in MacroModel.events] == expected_states

    @info h2("randomized driver/decay equivalence vs NaiveSEIR")
    rng = MersenneTwister(20260715)
    for _ in 1:1000
        E = rand(rng, 0:100)
        I = rand(rng, 0:100)
        x = (S = rand(rng, 0:100), E, I, R = rand(rng, 0:100))
        ellE = rand(rng, 0:E)
        ellI = rand(rng, 0:I)
        θ = (
            β = 5rand(rng), σ = 5rand(rng),
            γ = 5rand(rng), ω = 5rand(rng),
            ψ = rand(rng), χ = 0.0, N = 1.0 + rand(rng, 0:500),
        )
        old_driver, old_decay = existing_rates(x, θ, ellE, ellI)
        new_driver, new_decay = macro_rates(x, θ, ellE, ellI)
        @test new_driver ≈ old_driver
        @test new_decay ≈ old_decay

        χ = rand(rng)
        _, destructive_decay = existing_rates(x, θ, ellE, ellI; χ)
        @test destructive_decay - new_decay ≈ χ*I
    end
end

end
