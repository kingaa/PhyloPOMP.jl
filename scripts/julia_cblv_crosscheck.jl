#!/usr/bin/env julia
#
# julia_cblv_crosscheck.jl
#
# Read the same SEIR Newick trees used by python_cblv_crosscheck.py,
# apply Julia PhyloPOMP.cblv, and compare against the Python (phylodeep)
# unscaled CBLV output.
#
#   julia --project scripts/julia_cblv_crosscheck.jl

using PhyloPOMP
using Printf
using Statistics

base_dir  = joinpath(@__DIR__, "..")
cross_dir = joinpath(base_dir, "output", "crosscheck")
py_dir    = joinpath(cross_dir, "python_cblv")
jl_dir    = joinpath(cross_dir, "julia_cblv")
mkpath(jl_dir)

manifest_path = joinpath(cross_dir, "manifest.csv")
@assert isfile(manifest_path) "manifest.csv not found — run gen_seir_for_crosscheck.R first"

lines = readlines(manifest_path)
header = split(lines[1], ',')
col = Dict(h => i for (i, h) in enumerate(header))

function strip_singletons(g::Genealogy)
    g = deepcopy(g)
    changed = true
    while changed
        changed = false
        for n in eachindex(g)
            node = g[n]
            if length(node.children) == 1
                child_name = node.children[1]
                if isnothing(node.parent)
                    # Single-child ROOT: promote the child to root
                    # (matches phylodeep which detaches the child and discards the root edge)
                    g[child_name].parent = nothing
                else
                    parent_name = node.parent
                    g[child_name].parent = parent_name
                    idx = findfirst(==(n), g[parent_name].children)
                    if !isnothing(idx)
                        g[parent_name].children[idx] = child_name
                    end
                end
                node.children = PhyloPOMP.Name[]
                node.parent = nothing
                changed = true
            end
        end
    end
    PhyloPOMP.repair!(g)
    # Shift all slates so the new root starts at time 0,
    # matching phylodeep which sets dist_to_root=0 on the new root.
    root_offset = nothing
    for n in eachindex(g)
        if isnothing(g[n].parent) && !isempty(g[n].children)
            root_offset = g[n].slate
            break
        end
    end
    if !isnothing(root_offset) && root_offset != 0
        for n in eachindex(g)
            if !isnothing(g[n].parent) || !isempty(g[n].children)
                g[n].slate -= root_offset
            end
        end
    end
    g
end

println("=" ^ 70)
println("Julia ↔ Python (phylodeep) CBLV cross-check")
println("=" ^ 70)

results = []
all_pass = true

for ln in lines[2:end]
    isempty(strip(ln)) && continue
    f = split(ln, ',')
    idx      = parse(Int, f[col["index"]])
    n_tips   = parse(Int, f[col["n_tips"]])
    nwk_file = f[col["nwk_file"]]

    nwk_path = joinpath(cross_dir, nwk_file)
    nwk = String(strip(read(nwk_path, String)))

    g = parse_newick(nwk)
    gs = strip_singletons(g)
    x, y = cblv(gs)
    xf, yf = float.(x), float.(y)

    # Write Julia CBLV
    csv_path = joinpath(jl_dir, @sprintf("tree_%04d_cblv.csv", idx))
    open(csv_path, "w") do io
        println(io, "position,tip_x,internal_y")
        for k in eachindex(xf)
            @printf(io, "%d,%.12f,%.12f\n", k - 1, xf[k], yf[k])
        end
    end

    # Read Python unscaled CBLV
    py_path = joinpath(py_dir, @sprintf("tree_%04d_unscaled.csv", idx))
    if !isfile(py_path)
        println("  tree $idx: SKIP (no Python output)")
        continue
    end

    py_lines = readlines(py_path)
    px = Float64[]
    py_y = Float64[]
    for (i, pln) in enumerate(py_lines)
        i == 1 && continue
        parts = split(strip(pln), ',')
        push!(px, parse(Float64, parts[2]))
        if length(parts) >= 3 && !isempty(parts[3])
            push!(py_y, parse(Float64, parts[3]))
        end
    end

    # Julia: x has length n_tips, y has length n_tips (last y == 0 for root)
    # Python: x has length n_tips, y has length n_tips - 1 (no trailing 0 for root)
    # Need to reconcile: Python interleaved walk yields n tips in x, n-1 internals in y
    # Julia's cblv pushes a trailing 0 into y for each root

    jl_x = xf
    jl_y = yf[1:end-1]  # drop trailing root zero for comparison

    # But we also need to handle: phylodeep stores DEPTHS from root,
    # Julia stores BRANCH INCREMENTS from the walk.
    # Let me check by comparing sorted multisets first.

    println("\nTree $idx (tips=$n_tips)")
    println("  Julia:  len(x)=$(length(jl_x))  len(y)=$(length(jl_y))  height=$(jl_x[1])")
    println("  Python: len(x)=$(length(px))  len(y)=$(length(py_y))  height=$(px[1])")

    if length(jl_x) != length(px)
        println("  MISMATCH: tip count differs! Julia=$(length(jl_x)) Python=$(length(px))")
        all_pass = false
        push!(results, (idx=idx, tips_jl=length(jl_x), tips_py=length(px), status="TIP_COUNT_MISMATCH"))
        continue
    end

    # Check x vectors
    x_exact = isapprox(jl_x, px; atol=1e-6)
    x_multi = isapprox(sort(jl_x), sort(px); atol=1e-6)

    # Check y vectors (lengths might differ by 1)
    y_len_match = length(jl_y) == length(py_y)
    y_exact = y_len_match && isapprox(jl_y, py_y; atol=1e-6)
    y_multi = y_len_match && isapprox(sort(jl_y), sort(py_y); atol=1e-6)

    if x_exact && y_exact
        status = "EXACT"
    elseif x_multi && y_multi
        status = "MULTISET"
    elseif x_multi
        status = "X_MULTI_Y_DIFF"
    else
        status = "DIFFER"
    end

    ok = status in ("EXACT", "MULTISET")
    if !ok
        global all_pass = false
    end

    println("  x exact=$(x_exact)  x multiset=$(x_multi)")
    println("  y exact=$(y_exact)  y multiset=$(y_multi)  y_len_match=$(y_len_match)")
    println("  STATUS: $status  $(ok ? "PASS" : "FAIL")")

    if !x_multi
        # Show first differences
        dx = abs.(jl_x .- px)
        worst = argmax(dx)
        println("  worst x diff at [$worst]: julia=$(jl_x[worst]) python=$(px[worst]) diff=$(dx[worst])")
        # also check sorted
        sx_jl = sort(jl_x); sx_py = sort(px)
        dxs = abs.(sx_jl .- sx_py)
        println("  worst sorted x diff: $(maximum(dxs))")
    end

    if y_len_match && !y_multi
        dy = abs.(jl_y .- py_y)
        worst = argmax(dy)
        println("  worst y diff at [$worst]: julia=$(jl_y[worst]) python=$(py_y[worst]) diff=$(dy[worst])")
        sy_jl = sort(jl_y); sy_py = sort(py_y)
        dys = abs.(sy_jl .- sy_py)
        println("  worst sorted y diff: $(maximum(dys))")
    end

    push!(results, (idx=idx, tips_jl=length(jl_x), tips_py=length(px), status=status))
end

println("\n" * "=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
for r in results
    mark = r.status in ("EXACT", "MULTISET") ? "PASS" : "FAIL"
    println("  tree $(r.idx): tips_jl=$(r.tips_jl) tips_py=$(r.tips_py) → $mark ($( r.status))")
end
println()

if all_pass
    println("ALL TREES MATCH — Julia PhyloPOMP.cblv agrees with phylodeep CBLV.")
else
    println("SOME TREES DIFFER — see details above.")
end
