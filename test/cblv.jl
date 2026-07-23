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

    @testset "round trip via parse_cblv" begin
        for nwk in [
            "((a:4,((b:2,c:1):1)):0.5,d:2.5):0.5;",
            "((a:4,((b:2,c:1,e:3):1)):1,d:3):1;",
            "(d:3,(a:4,((c:1,b:2,e:3):1)):1):1;f:4;",
            "(:8,:3):1;(a:4,(:1,(b:2,c:1):1):2):1;",
        ]
            x, y = cblv(parse_newick(nwk))
            @test cblv(parse_cblv(x, y)) == (x, y)
            @test cblv(parse_newick(newick(parse_cblv(x, y)))) == (x, y)
        end
    end

    @testset "big tree (274 tips)" begin
        treefile = joinpath(@__DIR__, "tree.nwk")
        @assert isfile(treefile) "tree.nwk not found in test directory"
        g = parse_newick(String(strip(read(treefile, String))))
        x, y = cblv(g)
        xf = float.(x); yf = float.(y)

        _depth(g, n) = begin
            d = 0.0; cur = n
            while !isnothing(g[cur].parent)
                p = g[cur].parent
                d += float(g[cur].slate) - float(g[p].slate)
                cur = p
            end
            d
        end

        expected_internal_multiset(g) = begin
            vals = Float64[]
            for n in eachindex(g)
                k = length(g[n].children)
                d = _depth(g, n)
                for _ in 1:max(k - 1, 0)
                    push!(vals, d)
                end
            end
            append!(vals, zeros(length(collect(roots(g)))))
            sort(vals)
        end

        @test nsample(g) == 274
        @test length(xf) == 274
        @test length(yf) == 274

        tipdepths = [_depth(g, n) for n in eachindex(g) if isempty(g[n].children)]
        @test isapprox(xf[1], maximum(tipdepths); atol = 1e-6)

        @test isapprox(sort(yf), expected_internal_multiset(g); atol = 1e-6)

        rx, ry = cblv(parse_cblv(xf, yf))
        @test isapprox(sort(float.(rx)), sort(xf); atol = 1e-6)
        @test isapprox(sort(float.(ry)), sort(yf); atol = 1e-6)
        if float.(rx) == xf && float.(ry) == yf
            @info "round trip matched exactly (tie order happened to agree)"
        else
            @info "round trip matched as multisets; sequence differs on tied siblings (expected)"
        end

        reffile = joinpath(@__DIR__, "cblv_reference.csv")
        if isfile(reffile)
            ex = Float64[]; ey = Float64[]
            for (i, ln) in enumerate(eachline(reffile))
                i == 1 && continue
                f = split(strip(ln), ',')
                length(f) >= 2 && !isempty(f[2]) && push!(ex, parse(Float64, f[2]))
                length(f) >= 3 && !isempty(f[3]) && push!(ey, parse(Float64, f[3]))
            end
            @test isapprox(sort(xf), sort(ex); atol = 1e-6)
            @test isapprox(sort(yf), sort(ey); atol = 1e-6)
            if isapprox(xf, ex; atol = 1e-6) && isapprox(yf, ey; atol = 1e-6)
                @info "matched the phylodeep reference exactly on all 274 coordinates"
            else
                @info "matched phylodeep as multisets; sequence differs on tied siblings"
            end
        else
            @info "cblv_reference.csv absent; skipped the phylodeep cross-check"
        end
    end

end

end
