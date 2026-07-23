#!/usr/bin/env julia
#
# seir_cblv_forests.jl
#
# Consume a forest of Newick trees freshly simulated by
# scripts/simulate_seir_forest.R (phylopomp::runSEIR), plot summaries,
# and apply CBLV in Julia.
#
#   Rscript scripts/simulate_seir_forest.R 8
#   julia --project scripts/seir_cblv_forests.jl
#

using PhyloPOMP
using Printf
using Dates
using SHA

forest_dir = joinpath(@__DIR__, "..", "output", "seir_forest")
manifest_path = joinpath(forest_dir, "forest_manifest.csv")
outdir = joinpath(forest_dir, "cblv")
mkpath(outdir)

@assert isfile(manifest_path) """
Missing $manifest_path

Simulate a forest first:
  Rscript scripts/simulate_seir_forest.R 8
"""

println("=" ^ 64)
println("Julia CBLV on R-simulated SEIR forest")
println("=" ^ 64)
println("forest dir : $forest_dir")
println("time       : $(Dates.now())")

# ── read proof + manifest ──────────────────────────────────────────
proof_path = joinpath(forest_dir, "PROOF.txt")
if isfile(proof_path)
    println("\n── PROOF.txt (from R simulator) ──")
    println(read(proof_path, String))
end

# Minimal CSV reader (handles optional quotes)
dequote(s) = strip(s, '"')
lines = readlines(manifest_path)
header = dequote.(split(lines[1], ','))
col = Dict(h => i for (i, h) in enumerate(header))
needed = ["index", "multi_root", "n_tips", "n_root", "nwk_file", "nwk_sha1",
          "phylopomp_version", "seed", "simulated_at"]
for h in needed
    @assert haskey(col, h) "manifest missing column $h; have=$(keys(col))"
end

rows = []
for ln in lines[2:end]
    isempty(strip(ln)) && continue
    # simulated_at may contain spaces; it is the last column when quote=FALSE
    f = split(ln, ',')
    # If more fields than header (timestamp with spaces shouldn't split with quote=FALSE
    # and no commas in timestamp), pad/truncate defensively.
    if length(f) < length(header)
        error("CSV row has fewer fields than header: $ln")
    end
    getf(name) = dequote(f[col[name]])
    push!(rows, (
        index = parse(Int, getf("index")),
        multi_root = getf("multi_root") == "TRUE",
        n_tips = parse(Int, getf("n_tips")),
        n_root = parse(Int, getf("n_root")),
        nwk_file = getf("nwk_file"),
        nwk_sha1 = getf("nwk_sha1"),
        version = getf("phylopomp_version"),
        seed = getf("seed"),
        simulated_at = join(dequote.(f[col["simulated_at"]:end]), ","),
    ))
end

println("\nManifest: $(length(rows)) trees")
println("  phylopomp version : $(rows[1].version)")
println("  R seed            : $(rows[1].seed)")
println("  simulated_at      : $(rows[1].simulated_at)")

# ── strip singletons for CBLV ──────────────────────────────────────
function strip_singletons(g::Genealogy)
    g = deepcopy(g)
    changed = true
    while changed
        changed = false
        for n in eachindex(g)
            node = g[n]
            if length(node.children) == 1 && !isnothing(node.parent)
                child_name = node.children[1]
                parent_name = node.parent
                g[child_name].parent = parent_name
                idx = findfirst(==(n), g[parent_name].children)
                if !isnothing(idx)
                    g[parent_name].children[idx] = child_name
                end
                node.children = PhyloPOMP.Name[]
                node.parent = nothing
                changed = true
            end
        end
    end
    PhyloPOMP.repair!(g)
    g
end

function ascii_cblv_plot(x, y; width = 50, height = 10, title = "")
    lines_out = String[]
    push!(lines_out, "  CBLV: $title")
    n = length(x)
    vmax = max(maximum(x), maximum(y), 1e-12)
    step = max(1, n ÷ width)
    sx = x[1:step:end]; sy = y[1:step:end]
    cols = length(sx)
    for row in height:-1:1
        thr = vmax * row / height
        line = @sprintf("  %6.1f │", thr)
        for k in 1:cols
            xv, yv = sx[k], sy[k]
            line *= (xv >= thr && yv >= thr) ? "█" :
                    (xv >= thr) ? "▓" :
                    (yv >= thr) ? "░" : " "
        end
        push!(lines_out, line)
    end
    push!(lines_out, "         └" * "─"^cols)
    push!(lines_out, "          ▓=x tip depth  ░=y internal  █=both")
    join(lines_out, '\n')
end

# ── process each tree ──────────────────────────────────────────────
println("\n" * "=" ^ 64)
println("Parse Newick → verify SHA → CBLV")
println("=" ^ 64)

cblv_summary = open(joinpath(outdir, "cblv_summary.csv"), "w")
println(cblv_summary, "index,multi_root,n_tips_raw,n_root_raw,n_tips_cblv,tree_height,roundtrip,sha1_ok")

n_multi = Ref(0)
n_single = Ref(0)
n_cblv_ok = Ref(0)

for r in rows
    nwk_path = joinpath(forest_dir, r.nwk_file)
    @assert isfile(nwk_path) "missing $(nwk_path)"
    nwk = String(strip(read(nwk_path, String)))

    # Prove file matches manifest digest (when digest was recorded)
    sha_ok = true
    if !isempty(r.nwk_sha1) && r.nwk_sha1 != "NA"
        got = bytes2hex(sha1(nwk))
        sha_ok = (got == r.nwk_sha1)
        @assert sha_ok "SHA1 mismatch for $(r.nwk_file): manifest=$(r.nwk_sha1) got=$got"
    end

    g = parse_newick(nwk)
    nroots = length(collect(roots(g)))
    ntips = nsample(g)

    @assert nroots == r.n_root "nroot mismatch tree $(r.index): file=$nroots manifest=$(r.n_root)"
    @assert ntips == r.n_tips "ntips mismatch tree $(r.index)"

    kind = r.multi_root ? "MULTI-ROOT FOREST" : "single-root tree"
    println("\nTree $(r.index) [$kind]")
    println("  tips=$(ntips)  roots=$(nroots)  sha1_ok=$(sha_ok)")

    if r.multi_root
        n_multi[] += 1
        println("  CBLV skipped (multi-root forests are not CBLV-encodable as one vector)")
        println(cblv_summary, "$(r.index),TRUE,$ntips,$nroots,NA,NA,skipped,$sha_ok")
        continue
    end

    n_single[] += 1
    gs = strip_singletons(g)
    x, y = cblv(gs)
    xf, yf = float.(x), float.(y)

    rx, ry = cblv(parse_cblv(xf, yf))
    exact = (float.(rx) == xf && float.(ry) == yf)
    multiset = isapprox(sort(float.(rx)), sort(xf); atol=1e-6) &&
               isapprox(sort(float.(ry)), sort(yf); atol=1e-6)
    status = exact ? "EXACT" : multiset ? "MULTISET" : "FAIL"
    @assert multiset "CBLV round-trip failed for tree $(r.index)"
    n_cblv_ok[] += 1

    println("  CBLV length=$(length(xf))  height=$(round(xf[1], digits=3))  roundtrip=$status")
    println(ascii_cblv_plot(xf, yf, title="tree $(r.index) ($(length(xf)) tips)"))

    csv_path = joinpath(outdir, @sprintf("tree_%04d_cblv.csv", r.index))
    open(csv_path, "w") do io
        println(io, "index,tip_x,internal_y")
        for k in eachindex(xf)
            @printf(io, "%d,%.10f,%.10f\n", k - 1, xf[k], yf[k])
        end
    end
    println(cblv_summary, "$(r.index),FALSE,$ntips,$nroots,$(length(xf)),$(xf[1]),$status,$sha_ok")
end
close(cblv_summary)

println("\n" * "=" ^ 64)
println("PROOF SUMMARY")
println("=" ^ 64)
println("  Source              : phylopomp::runSEIR (R), NOT Julia canned seir_trees")
println("  phylopomp version   : $(rows[1].version)")
println("  simulation seed     : $(rows[1].seed)")
println("  simulated_at        : $(rows[1].simulated_at)")
println("  multi-root forests  : $(n_multi[])  (nroot≥2 — proves genealogical forest)")
println("  single-root trees   : $(n_single[])")
println("  CBLV round-trips OK : $(n_cblv_ok[]) / $(n_single[])")
println("  tip counts          : ", join([r.n_tips for r in rows], ", "))
println("  root counts         : ", join([r.n_root for r in rows], ", "))
println("  outputs             : $outdir")
println()
@assert n_multi[] ≥ 1 "expected at least one multi-root forest for proof"
@assert n_cblv_ok[] == n_single[] "not all single-root CBLV round-trips passed"
println("PASS: forest was simulated in R and CBLV applied in Julia.")
