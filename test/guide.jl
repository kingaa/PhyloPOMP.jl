using PhyloPOMP
using Test

@info "testing guides"
@testset verbose=true "guides" begin

    x = [
        "((((((:0.197946,[&&PhyloPOMP type=sample deme=I]:0.221360):0.160234,[&&PhyloPOMP type=sample deme=I]:0.572431)[&&PhyloPOMP type=branch deme=I]:0.192588)[&&PhyloPOMP type=sample deme=I]:0.033812):0.0913934,(([&&PhyloPOMP type=sample deme=I]:0.102642e-2,(([&&PhyloPOMP type=sample deme=I]:0.100027)[&&PhyloPOMP type=sample deme=I]:0.020000)[&&PhyloPOMP deme=I type=sample]:0.414031):0.161720,(([&&PhyloPOMP deme=I type=sample]:8.22165e-2,[&&PhyloPOMP deme=E type=node]:0.164473):0.015591)[&&PhyloPOMP deme=I]:0.449640):0.179927):0.0297577):0;"
    ];

    g = parse_newick(x,demes=SEIR,t0=0,time=1);
    m = fsmarkov(I=>0.5,E=>0.5,(I,E)=>5);
    z = guide(g,m)
    @test length(g)==length(z)
    @test occursin(r"^<guide:.*>$"s,sprint(show,z))
    @test occursin(r"^<t ∈ \[[\d\.]+,[\d\.]+\]:.*>$"s,sprint(show,z[19]))
    @test length(z[12:15])==4
    t = z[7].tend
    h1 = relhaz(z,t-0.1,7)
    h2 = relhaz(z,t-0.01,7)
    h3 = relhaz(z,t-0.001,7)
    @test_throws r"cannot evaluate at t" relhaz(z,t+0.1,7)
    @test all(keys(h1).===keys(h2))
    @test h1[(2,E=>I)]>4.0
    @test h2[(2,E=>I)]>40.0
    @test h3[(2,E=>I)]>400.0
    @test h2[(2,E=>I)]*h2[(2,I=>E)]==1
    @test isempty(relhaz(z,0.1,length(z)))

end
