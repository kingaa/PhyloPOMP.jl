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
    n3 = parse_newick("(:8,:3):1;((:1,(b:2,c:1):1):2,a:4):1;");
    x3 = cblv(n3)
    n4 = parse_newick("(:8,:3):1;(((b:2,c:1):1,:1):2,a:4):1;");
    x4 = cblv(n4)
    @test x1 == x2 == x3 == x4

    n = parse_newick("((a:4,((b:2,c:1):1)):0.5,d:2.5):0.5;");
    x,y = cblv(n);
    [x y]
    @test x==[5.0,3.0,1.0,2.5]
    @test y==[1.0,1.0,0.5,0.0]

    n = parse_newick("((a:4,((b:2,c:1):1)):0.5,d:2.5):0;");
    x,y = cblv(n);
    [x y]
    @test x==[4.5,3.0,1.0,2.5]
    @test y==[0.5,1.0,0.0,0.0]

    ## one ternary node, one zero-length branch
    n1 = parse_newick("((a:4,((b:2,c:1,e:3):1)):1,d:3):1;");
    x1,y1 = cblv(n1);
    [x1 y1]
    @test x1==[6.0,2.0,1.0,4.0,3.0]
    @test y1==[3.0,3.0,2.0,1.0,0.0]

    ## one ternary node, without zero-length branch
    n2 = parse_newick("((a:4,(b:2,c:1,e:3):1):1,d:3):1;");
    x2,y2 = cblv(n2);
    [x2 y2]
    @test x1 == x2 && y1 == y2

    n = parse_newick("(d:3,(a:4,((c:1,b:2,e:3):1)):1):1;f:4;");
    x,y = cblv(n);
    [x y]
    @test x==[6.0,2.0,1.0,4.0,3.0,4.0]
    @test y==[3.0,3.0,2.0,1.0,0.0,0.0]

end

end
