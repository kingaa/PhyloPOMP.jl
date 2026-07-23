# run_cblv_checks.jl
#
# Confirms the CBLV encoder on PhyloPOMP.jl devel (>= f5a5575) behaves correctly,
# using the 274-tip tree plus round-trip checks through parse_cblv.
#
#   julia --project run_cblv_checks.jl
#
# Requires PhyloPOMP.jl at devel >= 570c10e (parse_cblv). Files expected in the
# working directory:
#   tree.nwk            the 274-tip newick (shipped alongside)
#   cblv_reference.csv  optional: phylodeep ground truth for an external
#                       cross-check; if absent, that one check is skipped.
#
# The internal-node fix (push!(y,node.slate), commit f5a5575) is already upstream.
# This script verifies it rather than assuming it; check 3 below is the one that
# fails if you are on an older commit.

using PhyloPOMP
using Test

# --------------------------------------------------------------- helpers

# distance from n to the root of its tree, computed without touching cblv
_depth(g, n) = begin
    d = 0.0
    cur = n
    while !isnothing(g[cur].parent)
        p = g[cur].parent
        d += float(g[cur].slate) - float(g[p].slate)
        cur = p
    end
    d
end

# Expected internal channel as a multiset: a k-ary node contributes its depth
# (k-1) times, and each root adds one trailing zero separator. Order-free, so it
# is immune to ladderization tie-breaking.
function expected_internal_multiset(g)
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

read_ref(path) = begin
    ex = Float64[]; ey = Float64[]
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue
        f = split(strip(ln), ',')
        length(f) >= 2 && !isempty(f[2]) && push!(ex, parse(Float64, f[2]))
        length(f) >= 3 && !isempty(f[3]) && push!(ey, parse(Float64, f[3]))
    end
    ex, ey
end

# --------------------------------------------------------------- tests

@testset verbose = true "CBLV" begin

    @testset "golden cases" begin
        # permutation invariance across tree order and child order
        n1 = parse_newick("(a:4,(:1,(b:2,c:1):1):2):1;(:8,:3):1;")
        n2 = parse_newick("(:8,:3):1;(a:4,(:1,(b:2,c:1):1):2):1;")
        n3 = parse_newick("(:8,:3):1;((:1,(:2,:1):1):2,:4):1;")
        n4 = parse_newick("(:8,:3):1;(((:2,:1):1,:1):2,:4):1;")
        @test cblv(n1) == cblv(n2) == cblv(n3) == cblv(n4)

        # stemmed: the MRCA emits the 0.5 stem, not 0
        x, y = cblv(parse_newick("((a:4,((b:2,c:1):1)):0.5,d:2.5):0.5;"))
        @test x == [5.0, 3.0, 1.0, 2.5]
        @test y == [1.0, 2.0, 0.5, 0.0]

        # stem dropped to 0: every internal value falls by 0.5
        x, y = cblv(parse_newick("((a:4,((b:2,c:1):1)):0.5,d:2.5):0;"))
        @test x == [4.5, 3.0, 1.0, 2.5]
        @test y == [0.5, 1.5, 0.0, 0.0]

        # ternary node with a zero-length branch
        x1, y1 = cblv(parse_newick("((a:4,((b:2,c:1,e:3):1)):1,d:3):1;"))
        @test x1 == [6.0, 2.0, 1.0, 4.0, 3.0]
        @test y1 == [3.0, 3.0, 2.0, 1.0, 0.0]

        # same shape without the zero-length branch
        x2, y2 = cblv(parse_newick("((a:4,(b:2,c:1,e:3):1):1,d:3):1;"))
        @test x1 == x2 && y1 == y2

        # forest with a ternary node
        x, y = cblv(parse_newick("(d:3,(a:4,((c:1,b:2,e:3):1)):1):1;f:4;"))
        @test x == [6.0, 2.0, 1.0, 4.0, 3.0, 4.0]
        @test y == [3.0, 3.0, 2.0, 1.0, 0.0, 0.0]

        # distinct trees must not collide
        @test cblv(parse_newick("(:4,(:2,:1):1):3;")) != cblv(parse_newick("((:4,:3):2,:1):1;"))
    end

    @testset "round trip via parse_cblv" begin
        for nwk in [
            "((a:4,((b:2,c:1):1)):0.5,d:2.5):0.5;",
            "((a:4,((b:2,c:1,e:3):1)):1,d:3):1;",
            "(d:3,(a:4,((c:1,b:2,e:3):1)):1):1;f:4;",
            "(:8,:3):1;(a:4,(:1,(b:2,c:1):1):2):1;",
        ]
            x, y = cblv(parse_newick(nwk))
            # encode -> decode -> encode must be a fixed point
            @test cblv(parse_cblv(x, y)) == (x, y)
            # and survive a trip through newick as well
            @test cblv(parse_newick(newick(parse_cblv(x, y)))) == (x, y)
        end
    end

    @testset "forest is concatenation" begin
        # Trees parsed from separate newick strings each start at their own t0,
        # so the forest CBLV is the plain concatenation of the per-tree CBLVs.
        # This is the upstream-intended semantics; it also means no
        # between-introduction offset is carried on this path.
        n = [
            "((((:0.041,(:0.044,:0.32):0.62):0.2,((:0.35,:0.71):0.058,:0.21):0.54):0.064,(:0.54,:0.99):0.37):0.091);",
            "((((:0.33,:0.93):0.076,(:0.11,:0.27):0.57):0.039,(:0.2,:0.6):0.46):0.34);",
            "((((:0.091,:0.14):0.089,:0.1):0.037,:0.36):0.86);",
        ]
        x, y = cblv(parse_newick(n))
        x1, y1 = cblv(parse_newick(n[1]))
        x2, y2 = cblv(parse_newick(n[2]))
        x3, y3 = cblv(parse_newick(n[3]))
        @test x == [x1; x2; x3] && y == [y1; y2; y3]
    end

    @testset "big tree (274 tips)" begin
        @assert isfile("tree.nwk") "tree.nwk not found in the working directory"
        g = parse_newick(String(strip(read("tree.nwk", String))))
        x, y = cblv(g)
        xf = float.(x); yf = float.(y)

        # 1. shape
        @test nsample(g) == 274
        @test length(xf) == 274
        @test length(yf) == 274

        # 2. the first-visited tip is the deepest, so x[1] is the tree height
        tipdepths = [_depth(g, n) for n in eachindex(g) if isempty(g[n].children)]
        @test isapprox(xf[1], maximum(tipdepths); atol = 1e-6)

        # 3. THE REGRESSION CHECK. Order-free: the internal channel is exactly
        #    the multiset of internal-node depths. Passes with node.slate,
        #    fails with added_branch_length (244/273 coordinates wrong).
        @test isapprox(sort(yf), expected_internal_multiset(g); atol = 1e-6)

        # 4. round trip through the decoder.
        #    NOTE: compared as multisets, not sequences. This tree has 37 pairs
        #    of sibling subtrees with EXACTLY equal height (simultaneous
        #    sampling produces identical branch lengths), so ladderization has
        #    genuine ties. Tie order depends on node insertion order, which
        #    differs between the parsed tree and the decoded one, so exact
        #    sequence equality is not guaranteed even when both encodings are
        #    correct. The multiset comparison is tie-invariant.
        rx, ry = cblv(parse_cblv(xf, yf))
        @test isapprox(sort(float.(rx)), sort(xf); atol = 1e-6)
        @test isapprox(sort(float.(ry)), sort(yf); atol = 1e-6)
        if float.(rx) == xf && float.(ry) == yf
            @info "round trip matched exactly (tie order happened to agree)"
        else
            @info "round trip matched as multisets; sequence differs on tied siblings (expected)"
        end

        # 5. optional external cross-check against phylodeep.
        #    Also tie-sensitive: phylodeep breaks ties by first-max-found while
        #    scanning leaves, Julia's sortperm is stable on the children vector.
        #    Compared as multisets for the same reason as check 4; the exact
        #    sequence agreed when this reference was generated, so a sequence
        #    match is reported when it holds.
        if isfile("cblv_reference.csv")
            ex, ey = read_ref("cblv_reference.csv")
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
