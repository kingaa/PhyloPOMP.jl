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
    using .U: trans, recov, wane, T as UType
    
    x = [rcateg([1.0, 2.0, 3.0],UType)[1] for _ in 1:10000]
    t = tally(x)
    @test 1.9 < t[recov]/t[trans] < 2.1
    @test 2.9 < t[wane]/t[trans] < 3.1

    k,s,p = rcateg([0,0,0],true)
    @test k==1 && s==0 && p==0
    k,s = rcateg([0,0,0])
    @test k==1 && s==0 && p==0
    
end

end
