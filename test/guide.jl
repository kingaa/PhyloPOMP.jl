using PhyloPOMP
using Test

@info "testing guides"
@testset verbose=true "filter guides" begin

    x = [
        "((((((:0.197946,[&&PhyloPOMP type=sample deme=I]:0.221360):0.160234,[&&PhyloPOMP type=sample deme=I]:0.572431)[&&PhyloPOMP type=branch deme=I]:0.192588)[&&PhyloPOMP type=sample deme=I]:0.033812):0.0913934,(([&&PhyloPOMP type=sample deme=I]:0.102642e-2,(([&&PhyloPOMP type=sample deme=I]:0.100027)[&&PhyloPOMP type=sample deme=I]:0.020000)[&&PhyloPOMP deme=I type=sample]:0.414031):0.161720,(([&&PhyloPOMP deme=I type=sample]:8.22165e-2,[&&PhyloPOMP deme=E type=node]:0.164473):0.015591)[&&PhyloPOMP deme=I]:0.449640):0.179927):0.0297577):0;"
    ];

    @demes SEIRDemes E I
    g = parse_newick(x,demes=SEIRDemes,t0=0,time=1);
    z = guide(g,I=>0.5,E=>0.5,(I,E)=>5);
    @test length(g)==length(z)
    @test occursin(r"^<guide:.*>$"s,sprint(show,z))
    @test occursin(r"t ∈ \[[\d\.]+,[\d\.]+\]"s,sprint(show,z[19]))
    @test length(z[12:15])==4
    @test all(
        map(eachindex(g),eachindex(z)) do i,j
            length(g[i].children)==length(z[i].chillins)
        end
    )
    t = z[7].tend
    h1 = relhaz(z,t-0.1,E,7)
    h2 = relhaz(z,t-0.01,E,7)
    h3 = relhaz(z,t-0.001,E,7)
    @test all(keys(h1).===keys(h2))
    @test h1[(2,I)]>4.0
    @test h2[(2,I)]>40.0
    @test h3[(2,I)]>400.0
    @test isempty(relhaz(z,0.1,I,length(z)))
    @test_throws r"cannot evaluate at t" relhaz(z,t+0.1,I,7)

end
