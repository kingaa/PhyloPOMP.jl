module GenealogyPlotTest

import ..Main: h1, h2

@info h1("genealogy plot")

using Test
using Random: MersenneTwister
using PhyloPOMP
using PhyloPOMP: treelayout, MTBDDemes
import CairoMakie
using CairoMakie: Figure, Axis

## A root, then a branch point with one tip and one further branch point; three tips.
small_tree() = parse_newick("(((t1:0.5,(t2:1.0,t3:1.0)B:1.0)A:1.0)R:0);", t0 = 0.0, time = 4.0)

@testset verbose=true "genealogy plot" begin

    @info h2("layout: tips are 1:ntip, internal nodes sit between their children")
    g = small_tree()
    L = treelayout(g; ladderize = false)
    tips_ = tips(g)
    @test sort(L.y[tips_]) == collect(1.0:length(tips_))
    for i ∈ nodes(g)
        k = Int.(g[i].children)
        @test L.y[i] == (minimum(L.y[k]) + maximum(L.y[k])) / 2
    end
    @test length(L.branches) == length(g) - length(roots(g))
    @test length(L.connectors) == count(i -> !isempty(g[i].children), 1:length(g))
    for (x0, x1, y) ∈ L.branches
        @test x0 < x1
    end

    @info h2("layout: ladderize puts the smaller clade first (lowest y)")
    Lf = treelayout(g; ladderize = true)
    n = only(i for i ∈ nodes(g) if g[i].slate == minimum(g[j].slate for j ∈ nodes(g)))
    small, big = sort(Int.(g[n].children); by = c -> length(g[c].children))
    @test Lf.y[small] < Lf.y[big]

    @info h2("layout: a forest keeps each tree's tips together")
    gf = simulate(PhyloPOMP.SIR, (β = 4.0, γ = 1.0, ψ = 0.3, N = 100.0);
                  x0 = (S = 96, I = 4, R = 0), graft = [4], tmax = 3.0, rng = MersenneTwister(2))
    Lg = treelayout(gf)
    @test sort(Lg.y[tips(gf)]) == collect(1.0:length(tips(gf)))

    @info h2("layout: simulated trees")
    for s ∈ 1:20
        gs = simulate(PhyloPOMP.SIR, (β = 4.0, γ = 1.0, ψ = 0.3, N = 100.0);
                      x0 = (S = 99, I = 1, R = 0), graft = [1], tmax = 3.0, rng = MersenneTwister(s))
        nsample(gs) == 0 && continue
        Ls = treelayout(gs)
        @test sort(Ls.y[tips(gs)]) == collect(1.0:length(tips(gs)))
    end

    @info h2("genealogyplot returns a Figure (needs a Makie backend)")
    @test genealogyplot(g) isa Figure
    @test genealogyplot(g; ladderize = false, points = false, title = "t") isa Figure
    gd = simulate(PhyloPOMP.MTBD,
                  (lambda11 = 1.2, lambda12 = 0.3, lambda21 = 0.2, lambda22 = 0.9, m12 = 0.2, m21 = 0.1,
                   mu1 = 0.5, mu2 = 0.5, psi1 = 0.3, psi2 = 0.3, r1 = 0.7, r2 = 1.0);
                  x0 = (I1 = 1, I2 = 0), graft = [1, 0], tmax = 6.0, rng = MersenneTwister(4),
                  demeset = MTBDDemes, samplemap = [MTBDDemes.I1, MTBDDemes.I2])
    @test genealogyplot(gd) isa Figure
    @test genealogyplot(simulate(PhyloPOMP.SIR, (β = 4.0, γ = 1.0, ψ = 0.3, N = 100.0);
        x0 = (S = 99, I = 1, R = 0), graft = [1], tmax = 3.0, rng = MersenneTwister(1))) isa Figure

end

end
