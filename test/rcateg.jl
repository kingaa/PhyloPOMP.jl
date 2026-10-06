module RCategTest

import ..Main: h1, h2

@info h1("testing random draw from categorical distribution")

using PhyloPOMP
using Tally
using Random: seed!
import Random
using Test

## RNG with a pinned uniform draw, for the edge-case tests below.
struct PinnedRNG <: Random.AbstractRNG; u::Float64; end
Random.rand(r::PinnedRNG, ::Type{Float64}) = r.u

@testset verbose=true "rcateg tests" begin

    seed!(263261083)

    x = [rcateg([1.0, 2.0, 3.0])[1] for _ in 1:10000]
    t=tally(x)
    @test 1.9 < t[2]/t[1] < 2.1
    @test 2.9 < t[3]/t[1] < 3.1
    
    @marks U trans recov wane
    using .U: trans, recov, wane, T as UType
    
    x = [rcateg([1.0, 2.0, 3.0],UType)[1] for _ in 1:10000]
    t = tally(x)
    @test 1.9 < t[recov]/t[trans] < 2.1
    @test 2.9 < t[wane]/t[trans] < 3.1
    k,s = rcateg([0,0,0],UType)
    @test ismissing(k)

    b = BitSet([33, 45, 101])
    x = [rcateg([1.0, 2.0, 3.0],b)[1] for _ in 1:10000]
    t = tally(x)
    @test 1.9 < t[45]/t[33] < 2.1
    @test 2.9 < t[101]/t[33] < 3.1
    k,s = rcateg([0,0,0],b)
    @test ismissing(k)

    v = [33, 45, 101]
    x = [rcateg([1.0, 2.0, 3.0],v)[1] for _ in 1:10000]
    t = tally(x)
    @test 1.9 < t[45]/t[33] < 2.1
    @test 2.9 < t[101]/t[33] < 3.1
    k,s = rcateg([0,0,0],v)
    @test ismissing(k)

    k,s,p = rcateg([0,0,0],true)
    @test k==0 && s==0 && p==0
    k,s = rcateg([0,0,0])
    @test k==0 && s==0 && p==0

    @info h2("edge cases: rounding overshoot and r == 0 never select a zero weight")
    ## A pinned RNG hits the measure-zero branches: u -> 1 overshoots the last
    ## weight by rounding; u == 0 lands on a zero-weight first category.
    for (p, rng, want) ∈ [
            ([0.0, 0.0, 1.0], PinnedRNG(0.0), 3),
            ([0.0, 2.0, 1.0], PinnedRNG(0.0), 2),
            ([1.0, 0.0],      PinnedRNG(prevfloat(1.0)), 1),
            ([0.5, 0.5, 0.0], PinnedRNG(prevfloat(1.0)), 2),
        ]
        k, = rcateg(p; rng)
        @test k == want
    end
    ## fuzz: with u -> 1 the selected category must exist and have positive weight
    seed!(7)
    for _ ∈ 1:2000
        n = rand(1:6)
        p = [rand() < 0.3 ? 0.0 : rand() * 10.0^rand(-8:8) for _ ∈ 1:n]
        sum(p) > 0 || continue
        k, = rcateg(p; rng = PinnedRNG(prevfloat(1.0)))
        @test 1 ≤ k ≤ n && p[k] > 0
    end

end

end
