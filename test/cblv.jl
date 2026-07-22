module ParseTest

import ..Main: h1, h2

@info h1("testing CBLV representation")

using PhyloPOMP
using Test

@testset verbose=true "CBLV representation" begin

    n1 = parse_newick("(a:4,(:1,(b:2,c:1):1):2):1;(:8,:3):1;");
    x1 = cblv(n1)
    n2 = parse_newick("(:8,:3):1;(a:4,(:1,(b:2,c:1):1):2):1;");
    x2 = cblv(n2)
    n3 = parse_newick("(:8,:3):1;((:1,(:2,:1):1):2,:4):1;");
    x3 = cblv(n3)
    n4 = parse_newick("(:8,:3):1;(((:2,:1):1,:1):2,:4):1;");
    x4 = cblv(n4)
    x5 = cblv(parse_cblv(x4...));
    x6 = cblv(parse_newick(newick(parse_cblv(x5...))));
    @test x1 == x2 == x3 == x4 == x5 == x6

    n = parse_newick("((a:4,((b:2,c:1):1)):0.5,d:2.5):0.5;");
    x,y = cblv(n);
    [x y]
    @test x==[5.0,3.0,1.0,2.5]
    @test y==[1.0,2.0,0.5,0.0]

    n = parse_newick("((a:4,((b:2,c:1):1)):0.5,d:2.5):0;");
    x,y = cblv(n);
    [x y]
    @test x==[4.5,3.0,1.0,2.5]
    @test y==[0.5,1.5,0.0,0.0]

    ## one ternary node, one zero-length branch
    n1 = parse_newick("((a:4,((b:2,c:1,e:3):1)):1,d:3):1;");
    x1,y1 = cblv(n1);
    [x1 y1]
    @test x1==[6.0,2.0,1.0,4.0,3.0]
    @test y1==[3.0,3.0,2.0,1.0,0.0]
    x2,y2 = cblv(parse_cblv(x1,y1));
    @test x1 == x2 && y1 == y2

    ## one ternary node, without zero-length branch
    n2 = parse_newick("((a:4,(b:2,c:1,e:3):1):1,d:3):1;");
    x2,y2 = cblv(n2);
    [x2 y2]
    @test x1 == x2 && y1 == y2

    n = parse_newick("(d:3,(a:4,((c:1,b:2,e:3):1)):1):1;f:4;");
    x1,y1 = cblv(n);
    [x1 y1]
    @test x1==[6.0,2.0,1.0,4.0,3.0,4.0]
    @test y1==[3.0,3.0,2.0,1.0,0.0,0.0]
    x2,y2 = cblv(parse_cblv(x1,y1));
    @test x1 == x2 && y1 == y2

    n1 = parse_newick("(:4,(:2,:1):1):3;");
    n2 = parse_newick("((:4,:3):2,:1):1;");
    @test cblv(n1) != cblv(n2)

    n = [
        "((((:0.041,(:0.044,:0.32):0.62):0.2,((:0.35,:0.71):0.058,:0.21):0.54):0.064,(:0.54,:0.99):0.37):0.091);",
        "((((:0.33,:0.93):0.076,(:0.11,:0.27):0.57):0.039,(:0.2,:0.6):0.46):0.34);",
        "((((:0.091,:0.14):0.089,:0.1):0.037,:0.36):0.86);",
    ];
    x,y = cblv(parse_newick(n));
    x1,y1 = cblv(parse_newick(n[1]));
    x2,y2 = cblv(parse_newick(n[2]));
    x3,y3 = cblv(parse_newick(n[3]));
    @test x == [x1;x2;x3] && y == [y1;y2;y3]

    @test parse_cblv([3,2,1],[1,2,0]) isa Genealogy
    @test_throws "mismatch" parse_cblv([3,2,1],[1,0])
    @test_throws "empty" parse_cblv(Int[],Int[])
    @test_throws "improper" parse_cblv([1,2,3],[1,2,0])
    @test_throws "improper" parse_cblv([3,2,1],[1,5,0])
    @test_throws "improper" parse_cblv([3,-2,1],[1,2,0])
    @test_throws "improper" parse_cblv([3,2,1],[1,2,-1])
    @test_throws AssertionError parse_cblv([3,2,1],[1,2,0],demes=PhyloPOMP)
    @test_throws "exceeds" parse_cblv([3,2,1],[1,2,0],t0=5,time=6)
    @test_throws "invalid" parse_cblv([4.0,1,2],[1.0,0,3.5])
    @test_throws "invalid" parse_cblv([4.0,1,2],[1.0,0,0.5])

end

end
