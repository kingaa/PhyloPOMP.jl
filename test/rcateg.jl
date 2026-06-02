module RCategTest

import ..Main: h1, h2

@info h1("testing random draw from categorical distribution")

using PhyloPOMP
using Tally
using Random: seed!
using Test

@testset verbose=true "rcateg tests" begin

    seed!(263260083)

    x = [rcateg([1.0, 2.0, 3.0])[1] for _ in 1:10000]
    t=tally(x)
    @test 1.9 < t[2]/t[1] < 2.1
    @test 2.9 < t[3]/t[1] < 3.1
    
    @marks U trans recov wane
    x = [rcateg([1.0, 2.0, 3.0],U)[1] for _ in 1:10000]
    t = tally(x)
    @test 1.9 < t[U.recov]/t[U.trans] < 2.1
    @test 2.9 < t[U.wane]/t[U.trans] < 3.1
    
end

end
