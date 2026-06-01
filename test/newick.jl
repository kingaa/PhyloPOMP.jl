module NewickTest

import ..Main: h1, h2

@info h1("testing Newick formatter")

using PhyloPOMP
using Test

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

    @demes SEIR E I
    g = parse_newick(x,demes=SEIR,t0=0.0);
    n = newick(g);
    @test length(n)==length(x)
    @test all(isa.(n,Ref(String)))
    nn = (x->count(")",x)).(n);
    nx = (x->count(")",x)).(x);
    @test nn==reverse(nx)
    n = newick(g,extended=false,sigdigits=4);
    @test all(isa.(n,Ref(String)))

end

end
