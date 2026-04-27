using PhyloPOMP
using PhyloPOMP.FSMarkov
using Test

@info "testing finite-state Markov processes"
@testset verbose=true "finite-state Markov" begin

    @demes SI3R I1 I2 I3
    using .SI3R: I1, I2, I3
    p = fsmarkov(I1=>1,I2=>1,I3=>1);
    @test sum(statdist(p))==1
    @test sum(abs.(generator(p)))==0
    @test all(sum(generator(p),dims=1).==0)
    @test occursin(r"with generator",sprint(show,p))

    @test_throws UndefVarError fsmarkov()
    @test_throws r"unspecified stationary probability" fsmarkov(I1=>1,I3=>1)
    @test_throws r"unspecified stationary probability" fsmarkov(I1=>1,I3=>1,I3=>1)
    @test_throws r"non-positive stationary probability" fsmarkov(I1=>1,I3=>1,I2=>-1)

    p = fsmarkov(Float32,I1=>1,I3=>2,(I1,I2)=>3,(I2,I3)=>1,I2=>1);
    @test generator(p) isa Matrix{Float32}
    @test forward_action(p,0) isa Matrix{Float32}

    p = fsmarkov(I1=>1,I3=>2,(I1,I2)=>3,(I2,I3)=>1,I2=>1);
    @test sum(statdist(p))==1
    @test sum(abs.(generator(p)))==4.5
    @test all(sum(generator(p),dims=1).==0)
    @test_throws r"size mismatch" forward_action(p,10,[1,2])
    @test maximum(abs.(forward_action(p,0,[1,2,3])-[1,2,3]))+100-100==0
    @test maximum(abs.(forward_action(p,1e5,[1,1,1])-[0.75,0.75,1.5]))+100-100==0

    @test_throws r"negative conductance" fsmarkov(I1=>1,I3=>2,(I1,I2)=>3,I2=>1,(I2,I3)=>-1)
    @test_throws r"double specification of conductance" fsmarkov(I2=>1,I3=>2,(I1,I2)=>3,I1=>1,(I2,I1)=>3)
    @test_throws r"cannot specify conductance of a deme to itself" fsmarkov(I1=>1,I2=>1,I3=>2,(I1,I2)=>3,(I2,I2)=>3)

end
