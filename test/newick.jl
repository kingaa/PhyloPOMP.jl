using PhyloPOMP
using Test

@info "testing Newick parser"
@testset verbose=true "Newick parser" begin

    x1 = readlines("MERS_274_sCoal_phylopomp.nwk");
    @enum MERSDemes Camel Human
    @test_throws r"final time from data .* exceeds" parse_newick(x1,Val(MERSDemes),0.0,1)
    @time g1 = parse_newick(x1,Val(MERSDemes),0.0,6);
    @test isa(g1,Genealogy)
    @test g1.time == 6

    @enum Strains Strain1 Strain2 Strain3
    x2 = readlines("B.1.617.all.nwk");
    @time g2 = parse_newick(x2,Val(Strains),1.0);
    @test isa(g2,Genealogy)
    @test g2.t0 == 1.0

    x3 = [
        "((50:0.195305,(64:0.155799,81:0.000000)55:0.308953)41:1.296678)1:0.000000;(65:1.762283)4:0.000000;(72:1.885165)5;((((((76:0.534765,80:0.000000)46:0.253820,67:0.618309)36:0.375403,73:1.107300)27:0.060102,79:0.000000)24:0.652522,26:0.706297)10:0.095550)7;",
        "((((((((([&&PhyloPOMP type=sample deme=I]33:0.121768)[&&PhyloPOMP type=node deme=I]30:0.256093)[&&PhyloPOMP type=sample deme=I]23:0.048347)[&&PhyloPOMP type=node deme=I]21:0.000896)[&&PhyloPOMP type=node deme=I]20:0.272773)[&&PhyloPOMP type=node deme=I]14:0.039280)[&&PhyloPOMP type=node deme=I]13:0.100483)[&&PhyloPOMP type=node deme=I]12:0.001259)[&&PhyloPOMP type=node deme=I]11:0.058041)[&&PhyloPOMP type=root deme=0]9:0.000000;(([&&PhyloPOMP type=sample deme=I]22:0.218679)[&&PhyloPOMP type=node deme=I]16:0.269299)[&&PhyloPOMP type=root deme=0]8:0.000000;(([&&PhyloPOMP type=sample deme=I]18:0.291933)[&&PhyloPOMP type=sample deme=I]10:0.039330)[&&PhyloPOMP type=root deme=0]5:0.000000;(((([&&PhyloPOMP type=sample deme=I]39:0.011566)[&&PhyloPOMP type=sample deme=I]38:0.111957)[&&PhyloPOMP type=sample deme=I]32:0.089940)[&&PhyloPOMP type=node deme=E]28:0.759231)[&&PhyloPOMP type=root deme=0]4:0.000000;",
    ];

    @enum SEIR E I
    @test_throws "invalid base 10 digit" parse_newick(x3,Val(SEIR),5)
    @time g3 = parse_newick(x3,Val(SEIR),5.0,7.0);
    @test isa(g3,Genealogy)
    @test ismissing(g3.nodes[3].deme)
    @test sum(map(x->ismissing(x.deme),g3.nodes))==23
    @test sum(map(x->!ismissing(x.deme),g3.nodes))==17
    @test sum(map(x->x.deme===I,g3.nodes))==16
    @test sum(map(x->x.deme===E,g3.nodes))==1
    @test length(g3.nodes)==40
    @test sum(map(x->x.type,g3.nodes).==PhyloPOMP.Sample)==19
    @test sum(map(x->x.type,g3.nodes).==PhyloPOMP.Node)==21

    g4=parse_newick("():0.1;",Val(Strains),0.0);
    @test length(g4.nodes)==2

    @test_throws "unbalanced parentheses" parse_newick("(:0.1;",Val(Strains),0.0)
    @test_throws "unbalanced square brackets" parse_newick("[bob=3 tom=[]:0.1;",Val(Strains),0.0)

end
