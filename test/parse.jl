using PhyloPOMP
using Test

@info "testing Newick parser"
@testset verbose=true "Newick parser" begin

    x = [
        "((50:0.195305,(64:0.155799,81:0.000000)55:0.308953)41:1.296678)1:0.000000;(65:1.762283)4:0.000000;(72:1.885165)5;((((((76:0.534765,80:0.000000)46:0.253820,67:0.618309)36:0.375403,73:1.107300)27:0.060102,79:0.000000)24:0.652522,26:0.706297)10:0.095550)7;",
        "(((((((((((((((([&&PhyloPOMP type=sample deme=I]101:0.348292)[&&PhyloPOMP type=sample deme=I]71:0.043151)[&&PhyloPOMP type=node deme=I]64:0.248533)[&&PhyloPOMP type=node deme=I]49:0.009989)[&&PhyloPOMP type=node deme=E]46:0.005282,((([&&PhyloPOMP type=sample deme=I]75:0.053987)[&&PhyloPOMP type=node deme=I]68:0.070730)[&&PhyloPOMP type=node deme=I]59:0.191481)[&&PhyloPOMP type=node deme=I]51:0.029079)[&&PhyloPOMP type=node deme=I]42:0.186568,((([&&PhyloPOMP type=sample deme=I]74:0.254742)[&&PhyloPOMP type=node deme=I]56:0.059999)[&&PhyloPOMP type=node deme=I]48:0.000421)[&&PhyloPOMP type=node deme=E]47:0.201246)[&&PhyloPOMP type=node deme=I]25:0.166445)[&&PhyloPOMP type=node deme=I]21:0.135240)[&&PhyloPOMP type=node deme=I]15:0.175788)[&&PhyloPOMP type=node deme=I]9:0.134695,(((((((((([&&PhyloPOMP type=sample deme=I]102:0.027845)[&&PhyloPOMP type=sample deme=I]100:0.629270)[&&PhyloPOMP type=node deme=I]50:0.016357)[&&PhyloPOMP type=node deme=I]43:0.025963)[&&PhyloPOMP type=node deme=I]41:0.002351,((([&&PhyloPOMP type=sample deme=I]107:0.213945)[&&PhyloPOMP type=node deme=I]91:0.038432)[&&PhyloPOMP type=node deme=I]88:0.172176)[&&PhyloPOMP type=node deme=E]70:0.325434)[&&PhyloPOMP type=node deme=I]40:0.025929)[&&PhyloPOMP type=node deme=I]37:0.012487)[&&PhyloPOMP type=node deme=I]35:0.004896)[&&PhyloPOMP type=node deme=E]34:0.485078)[&&PhyloPOMP type=node deme=I]12:0.053967)[&&PhyloPOMP type=node deme=E]10:0.189925)[&&PhyloPOMP type=node deme=I]6:0.095254)[&&PhyloPOMP type=node deme=I]5:0.000298)[&&PhyloPOMP type=sample deme=I]4:0.044483)[&&PhyloPOMP type=node deme=E]3:0.066014,((((([&&PhyloPOMP type=sample deme=I]30:0.192045,(((([&&PhyloPOMP type=sample deme=I]106:0.135813)[&&PhyloPOMP type=node deme=I]95:0.097167)[&&PhyloPOMP type=node deme=E]89:0.087857)[&&PhyloPOMP type=node deme=I]83:0.397810)[&&PhyloPOMP type=node deme=E]44:0.296515)[&&PhyloPOMP type=node deme=I]22:0.109389)[&&PhyloPOMP type=node deme=I]20:0.009961)[&&PhyloPOMP type=node deme=I]19:0.279611)[&&PhyloPOMP type=node deme=I]8:0.045107,((((([&&PhyloPOMP type=sample deme=I]92:0.202637)[&&PhyloPOMP type=node deme=E]72:0.477023)[&&PhyloPOMP type=node deme=I]26:0.270156,(((((((([&&PhyloPOMP type=sample deme=I]82:0.029926)[&&PhyloPOMP type=node deme=I]80:0.071942)[&&PhyloPOMP type=node deme=E]69:0.031919)[&&PhyloPOMP type=node deme=I]65:0.028600)[&&PhyloPOMP type=sample deme=I]61:0.148185)[&&PhyloPOMP type=node deme=E]57:0.031422)[&&PhyloPOMP type=node deme=I]54:0.136358)[&&PhyloPOMP type=node deme=I]33:0.076175)[&&PhyloPOMP type=node deme=E]27:0.280337)[&&PhyloPOMP type=node deme=I]17:0.023589)[&&PhyloPOMP type=node deme=I]16:0.079893)[&&PhyloPOMP type=node deme=E]13:0.197383)[&&PhyloPOMP type=node deme=I]7:0.267953)[&&PhyloPOMP type=node deme=I]2:0.063324)[&&PhyloPOMP type=node deme=I]1:0.077664)[&&PhyloPOMP type=root]0:0;",
    ];

    @demes SEIR E I
    @test_throws r"final time from data \(.+\) exceeds" parse_newick(x,demes=SEIR,t0=5.0,time=6);
    g = parse_newick(x,demes=SEIR,t0=5.0,time=7.0);
    @test g isa Genealogy{SEIR.T}
    @test ismissing(g.nodes[3].deme)
    @test sum(map(x->ismissing(x.deme),g.nodes))==20
    @test sum(map(x->!ismissing(x.deme),g.nodes))==65
    @test sum(map(x->x.deme===SEIR.I,g.nodes))==52
    @test sum(map(x->x.deme===SEIR.E,g.nodes))==13
    @test length(g.nodes)==85
    @test sum(map(x->x.type,g.nodes).==PhyloPOMP.Sample)==24
    @test sum(map(x->x.type,g.nodes).==PhyloPOMP.Node)==56
    @test sum(map(x->x.type,g.nodes).==PhyloPOMP.Root)==5
    @test all((x->x.type).(g[PhyloPOMP.tips(g)]).==PhyloPOMP.Sample)
    @test all(PhyloPOMP.tips(g) .∈ Ref(PhyloPOMP.samples(g)))
    @test isempty(intersect(PhyloPOMP.tips(g),PhyloPOMP.nodes(g)))
    @test isempty(intersect(PhyloPOMP.roots(g),PhyloPOMP.nodes(g)))
    @test occursin(r"^<genealogy on .*>>"s,sprint(show,g))
    @test length(collect(eachmatch(r"time",sprint(show,g))))==length(g)
    @test occursin(r"lineage=7 .* parent=10",sprint(show,g[11]))
    @test PhyloPOMP.nsample(g)==24

    g=parse_newick("():0.1;",t0=0.0);
    @test length(g.nodes)==2
    @test g.nodes[2].type==PhyloPOMP.Sample
    g = parse_newick("A():3.2;",t0=0.0);
    @test g.nodes[2].type==PhyloPOMP.Sample
    g = parse_newick("A([&&PhyloPOMP type=node]):3.2;",t0=0.0);
    @test g.nodes[2].type==PhyloPOMP.Sample
    @test_warn "dropping zero-length branch" g = parse_newick("A()[&&PhyloPOMP type=sample]:3.2;",t0=0.0);
    @test g.nodes[2].type==PhyloPOMP.Sample
    g = parse_newick("A:4;()B:3;",t0=0.0,time=5.0);
    @test length(g.nodes)==4
    g = parse_newick("(:[Bob:yes]2.5)A:3.2[&&NHX:type=node];",t0=0.0);
    @test length(g)==3
    @test g.time==5.7
    @test g.nodes[3].type==PhyloPOMP.Sample
    @test occursin(r"^<genealogy on .*>>"s,sprint(show,g))

    @test_throws "unbalanced parentheses" parse_newick("(:0.1;",demes=SEIR,t0=0.0)
    @test_throws "unbalanced square brackets" parse_newick(")3:1;(((:0.1)),[&&PhyloPOMP:deme=E]type=sample]:1.00,(((:0.3,:0.1),),):0.3)a:0.5;",demes=SEIR,t0=0.0)
    @test_throws "unbalanced square brackets" parse_newick("yloPOMP:deme=E|type=node]:1.00,(((:0.3,:0.1),),):0.3)a:0.5;",demes=SEIR,t0=0.0)
    @test_throws "unbalanced parentheses" parse_newick(")3:1;(((:0.1)),[&&PhyloPOMP:deme=E|type=sample]:1.00,(((:0.3,:0.1),),):0.3)a:0.5;",demes=SEIR,t0=0.0)
    @test_throws "unbalanced parentheses" parse_newick(")3:1;(((:0.1)),[&&PhyloPOMP:deme=E|type=sample]:1.00,(((:0.3,:0.1),),):0.3)a:0.5;",demes=SEIR,t0=0.0)
    @test_throws "misplaced comma or" parse_newick("((:0.1)),([&&PhyloPOMP:deme=E|type=sample]:1.00,(((:0.3,:0.1),),):0.3)a:0.5;",demes=SEIR,t0=0.0)
    @test_throws "unrecognized deme" parse_newick("([&&PhyloPOMP:deme=tim|type=bob]:1.0000):0.5;",demes=SEIR,t0=0.0)
    @test_throws "unrecognized type" parse_newick("([&&PhyloPOMP:deme=E|type=bob]:1.0000):0.5;",demes=SEIR,t0=0.0)
    @test_throws "unrecognized type" parse_newick("(()[&&PhyloPOMP:deme=I|type=extant]:1.0000):0.5;",demes=SEIR,t0=0.0)
    @test_throws "no final semicolon" parse_newick("([&&PhyloPOMP|deme=9|type=sample]:1.0000):0.5",demes=SEIR,t0=0.0)
    @test_throws r"cannot parse .* as Float64" parse_newick("([&&PhyloPOMP|deme=E|type=sample]:1.0000A):0.5;",demes=SEIR,t0=0.0)
    @test_throws "unbalanced parentheses" parse_newick("(:1,:1):1;((:1,:1):1;",demes=SEIR,t0=0.0)
    @test_throws "missing comma or semicolon" parse_newick("()();",demes=SEIR,t0=0.0)
    @test_throws "unbalanced parentheses" parse_newick("(A();",demes=SEIR,t0=0.0)
    @test_throws "misplaced comma or unbalanced parentheses" parse_newick(",();",demes=SEIR,t0=0.0)
    @test_throws "misplaced comma or unbalanced parentheses" parse_newick(":2 , 3 ();",demes=SEIR,t0=0.0)
    @test_throws "misplaced colon" parse_newick(":2 ,:3 ();",demes=SEIR,t0=0.0)
    @test_throws "unbalanced square brackets" parse_newick("()[bob =3 jack=[4]:0.3;",demes=SEIR,t0=0.0)
    @test_throws "negative branch length" parse_newick("():-32.5;")
    @test_warn "zero branch-length" parse_newick("():;")
    @test_throws r"cannot parse .* as Float64" parse_newick("():A;")

end
