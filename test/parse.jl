using PhyloPOMP
using Test

@info "testing Newick parser"
@testset verbose=true "Newick parser" begin

    x = [
        "((50:0.195305,(64:0.155799,81:0.000000)55:0.308953)41:1.296678)1:0.000000;(65:1.762283)4:0.000000;(72:1.885165)5;((((((76:0.534765,80:0.000000)46:0.253820,67:0.618309)36:0.375403,73:1.107300)27:0.060102,79:0.000000)24:0.652522,26:0.706297)10:0.095550)7;",
        "((((((((([&&PhyloPOMP type=sample deme=I]33:0.121768)[&&PhyloPOMP type=node deme=I]30:0.256093)[&&PhyloPOMP type=sample deme=I]23:0.048347)[&&PhyloPOMP type=node deme=I]21:0.000896)[&&PhyloPOMP type=node deme=I]20:0.272773)[&&PhyloPOMP type=node deme=I]14:0.039280)[&&PhyloPOMP type=node deme=I]13:0.100483)[&&PhyloPOMP type=node deme=I]12:0.001259)[&&PhyloPOMP type=node deme=I]11:0.058041)[&&PhyloPOMP type=root]9:0.000000;(([&&PhyloPOMP type=sample deme=I]22:0.218679)[&&PhyloPOMP type=node deme=I]16:0.269299)[&&PhyloPOMP type=root]8:0.000000;(([&&PhyloPOMP type=sample deme=I]18:0.291933)[&&PhyloPOMP type=sample deme=I]10:0.039330)[&&PhyloPOMP type=root]5:0.000000;(((([&&PhyloPOMP type=sample deme=I]39:0.011566)[&&PhyloPOMP type=sample deme=I]38:0.111957)[&&PhyloPOMP type=sample deme=I]32:0.089940)[&&PhyloPOMP type=node deme=E]28:0.759231)[&&PhyloPOMP type=root]4:0.000000;",
    ];

    @enum SEIR E I
    @test_throws r"final time from data \(.+\) exceeds" parse_newick(x,demes=SEIR,t0=5.0,time=6);
    @time g = parse_newick(x,demes=SEIR,t0=5.0,time=7.0);
    @test isa(g,Genealogy)
    @test isa(g,Genealogy{SEIR})
    @test ismissing(g.nodes[3].deme)
    @test sum(map(x->ismissing(x.deme),g.nodes))==23
    @test sum(map(x->!ismissing(x.deme),g.nodes))==17
    @test sum(map(x->x.deme===I,g.nodes))==16
    @test sum(map(x->x.deme===E,g.nodes))==1
    @test length(g.nodes)==40
    @test sum(map(x->x.type,g.nodes).==PhyloPOMP.Sample)==19
    @test sum(map(x->x.type,g.nodes).==PhyloPOMP.Node)==13
    @test sum(map(x->x.type,g.nodes).==PhyloPOMP.Root)==8
    @test all((x->x.type).(g[PhyloPOMP.tips(g)]).==PhyloPOMP.Sample)
    @test all(PhyloPOMP.tips(g) .∈ Ref(PhyloPOMP.samples(g)))
    @test isempty(intersect(PhyloPOMP.tips(g),PhyloPOMP.nodes(g)))
    @test isempty(intersect(PhyloPOMP.roots(g),PhyloPOMP.nodes(g)))
    @test occursin(r"^<genealogy on .*>>"s,sprint(show,g))
    @test length(collect(eachmatch(r"time",sprint(show,g))))==length(g)
    @test occursin(r"lineage=4 .* parent=10",sprint(show,g[11]))

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

@info "testing Newick formatter"
@testset verbose=true "Newick formatter" begin

    x = [
        "((50:0.195305,(64:0.155799,81:0.000000)55:0.308953)41:1.296678)1:0.000000;"
        "(65:1.762283)4:0.000000;"
        "(72:1.885165)5;"
        "((((((76:0.534765,80:0.000000)46:0.253820,67:0.618309)36:0.375403,73:1.107300)27:0.060102,79:0.000000)24:0.652522,26:0.706297)10:0.095550)7;"
        "((((((((([&&PhyloPOMP type=sample deme=I]33:0.121768)[&&PhyloPOMP type=node deme=I]30:0.256093)[&&PhyloPOMP type=sample deme=I]23:0.048347)[&&PhyloPOMP type=node deme=I]21:0.000896)[&&PhyloPOMP type=node deme=I]20:0.272773)[&&PhyloPOMP type=node deme=I]14:0.039280)[&&PhyloPOMP type=node deme=I]13:0.100483)[&&PhyloPOMP type=node deme=I]12:0.001259)[&&PhyloPOMP type=node deme=I]11:0.058041)[&&PhyloPOMP type=root]9:0.000000;"
        "(([&&PhyloPOMP type=sample deme=I]22:0.218679)[&&PhyloPOMP type=node deme=I]16:0.269299)[&&PhyloPOMP type=root]8:0.000000;"
        "(([&&PhyloPOMP type=sample deme=I]18:0.291933)[&&PhyloPOMP type=sample deme=I]10:0.039330)[&&PhyloPOMP type=root]5:0.000000;"
        "(((([&&PhyloPOMP type=sample deme=I]39:0.011566)[&&PhyloPOMP type=sample deme=I]38:0.111957)[&&PhyloPOMP type=sample deme=I]32:0.089940)[&&PhyloPOMP type=node deme=E]28:0.759231)[&&PhyloPOMP type=root]4:0.000000;"
    ];

    @time g = parse_newick(x,demes=SEIR,t0=5.0);
    @time n = newick(g);
    @test length(n)==length(x)
    @test all(isa.(n,Ref(String)))
    nn = (x->count(")",x)).(n);
    nx = (x->count(")",x)).(x);
    @test nn==reverse(nx)
    @time n = newick(g,extended=false,sigdigits=4);
    @test all(isa.(n,Ref(String)))

end
