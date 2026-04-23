using PhyloPOMP
using Test

@info "testing finite-state Markov processes"
@testset verbose=true "finite-state Markov" begin

    @demes SI3R I1 I2 I3

    g = fsmarkov(I1=>1,I2=>1,I3=>1);
    @test sum(g.statdist)==1
    @test sum(abs.(g.generator))==0
    @test all(sum(g.generator,dims=1).==0)

    @test_throws UndefVarError fsmarkov()
    @test_throws r"unspecified stationary probability" fsmarkov(I1=>1,I3=>1)
    @test_throws r"unspecified stationary probability" fsmarkov(I1=>1,I3=>1,I3=>1)
    @test_throws r"non-positive stationary probability" fsmarkov(I1=>1,I3=>1,I2=>-1)

    g = fsmarkov(I1=>1,I3=>2,(I1,I2)=>3,(I2,I3)=>1,I2=>1);
    @test sum(g.statdist)==1
    @test sum(abs.(g.generator))==4.5
    @test all(sum(g.generator,dims=1).==0)

    @test_throws r"negative conductance" fsmarkov(I1=>1,I3=>2,(I1,I2)=>3,I2=>1,(I2,I3)=>-1)
    @test_throws r"double specification of conductance" fsmarkov(I2=>1,I3=>2,(I1,I2)=>3,I1=>1,(I2,I1)=>3)
    @test_throws r"cannot specify conductance of a deme to itself" fsmarkov(I1=>1,I2=>1,I3=>2,(I1,I2)=>3,(I2,I2)=>3)

    g = fsmarkov(Float32,I1=>1,I3=>2,(I1,I2)=>3,(I2,I3)=>1,I2=>1);

end
