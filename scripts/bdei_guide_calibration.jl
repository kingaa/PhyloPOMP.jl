# Guide calibration for the generic BDEI filter (`mgp_filter_pomp`), judged against the exact likelihood (`mtbd_exact`).
#
#   julia -t 12 --project=pipeline scripts/bdei_guide_calibration.jl pilot [seed = 11] [ntips = 20] [nrep = 24]
#   julia -t 12 --project=pipeline scripts/bdei_guide_calibration.jl grid [seed = 11] [ntips = 20] [nrep = 24]
#   julia -t 12 --project=pipeline scripts/bdei_guide_calibration.jl scale <out.tsv> [sizes = 20,50,80,150]
#   julia -t 12 --project=pipeline scripts/bdei_guide_calibration.jl moderate <out.tsv>
#
# A guide over (E, I) with backward-time switching rates a (E → I) and b (I → E) is
# fsmarkov(E => b, I => a, (E, I) => a + b): make_generator gives Q = Diagonal(π)·C, so the rate I → E is π_E·(a + b) = b
# and E → I is π_I·(a + b) = a. Candidate from θ: in an epidemic growing at rate r, an E host's age is Exp(σ + r) and an I
# host's age Exp(γ + r), so a lineage traced back leaves E at a = σ + r (to its infector, in I) and I at b = γ + r (to E).
# The default guide of the earlier experiments, fsmarkov(E => 0.1, I => 1, (E, I) => 1), is a = 1/1.1, b = 0.1/1.1.
#
# Trees: the BDEI design of generate_trees.jl (stop at the t-th sample), p ≥ 0.1 to keep the hidden population small.
# Per filter: nrep runs at Np particles; mean exp(ll − L) ± SE (1 = unbiased), sd(ll), median(ll) − L (the typical
# shortfall), runs with a finite ll. The filters are not seeded, so reruns agree in distribution, not digit for digit.
#
# Results (2026-10-09, Np 1000 unless noted; median(ll) − L):
#   20 tips: naive −3.3 to −18; θ guide −0.03 to −1.8 (one tree with R0 4.7: −12, and −1.6 at Np 4000); sd ≈ 3× smaller.
#   grid (20 tips): the θ point lies in the region of smallest sd (0.8–1.3); a × 2 or × 4 is clearly worse.
#   50–80 tips, R0 1.6–4.6 (9 trees): every filter fails; θ guide best, −14 to −160 (Np 4000: −4 to −83); naive −Inf or worse.
#   For comparison, the naive LBDP filter at R0 4.5 has sd(ll) 0.3–0.9 on 50 and 120 tips: the cost is BDEI's E/I
#   lineage structure, not the size of the hidden population.

using PhyloPOMP, Printf, Random, Statistics
import PartiallyObservedMarkovProcesses as POMP
include(joinpath(@__DIR__, "..", "pipeline", "generate_trees.jl"))      # sim_tips, draw
const D = PhyloPOMP.MGPDemes2

growth(θ) = (σ = θ.σ; λ = θ.λ; γ = θ.μ + θ.χ; (-(σ + γ) + sqrt((σ - γ)^2 + 4σ * λ)) / 2)
guide_ab(a, b) = fsmarkov(D.d1 => b, D.d2 => a, (D.d1, D.d2) => a + b)
theory_ab(θ) = (r = growth(θ); (θ.σ + r, θ.μ + θ.χ + r))
const DEFAULT = (1 / 1.1, 0.1 / 1.1)

function runs(g, θ, L; guide = nothing, proposal = :guided, Np = 1000, nrep = 16)
    P = mgp_filter_pomp(g, PhyloPOMP.BDEI; θ, x0 = (E = 0, I = 1), guide, proposal, maxpop = 100 * nsample(g) + 1000)
    ll = zeros(nrep)
    t = @elapsed Threads.@threads :dynamic for k in 1:nrep
        ll[k] = POMP.logLik(POMP.pfilter(P, Np = Np))
    end
    fin = filter(isfinite, ll)
    w = exp.(ll .- L)
    (m = mean(w), se = std(w) / sqrt(nrep), sd = length(fin) > 1 ? std(fin) : NaN,
     med = median(ll) - L, nfin = length(fin), t = t * min(Threads.nthreads(), nrep) / nrep)
end
show_row(name, s) = @printf("  %-34s mean exp(ll-L) %6.3f ± %5.3f  sd(ll) %6.2f  median(ll)-L %+8.2f  finite %2d  %5.1f s/run\n",
                            name, s.m, s.se, s.sd, s.med, s.nfin, s.t)

function tree_of(ntips, rng; pmin = 0.1)
    while true
        targets, p, θ, x0, graft = draw(Val(:BDEI), rng)
        p < pmin && continue
        g = sim_tips(θ, ntips, rng; M = PhyloPOMP.BDEI, x0, graft)
        g === nothing || return g, θ, targets, p
    end
end

describe(g, θ, targets, p, L) =
    @printf("tree %d tips: R0 %.2f IP %.2f INC %.2f p %.2f; σ %.3f γ %.3f r %.3f; θ guide a %.3f b %.3f; L = %.3f\n",
            nsample(g), targets..., p, θ.σ, θ.μ + θ.χ, growth(θ), theory_ab(θ)..., L)

function pilot(seed, ntips, nrep)
    g, θ, targets, p = tree_of(ntips, Xoshiro(seed))
    L = mtbd_exact(g, PhyloPOMP.BDEI; θ, founder = :I)
    describe(g, θ, targets, p, L)
    show_row("naive", runs(g, θ, L; nrep))
    for prop in (:guided, :soft, :hard)
        show_row("default $prop", runs(g, θ, L; guide = guide_ab(DEFAULT...), proposal = prop, nrep))
        show_row("θ guide $prop", runs(g, θ, L; guide = guide_ab(theory_ab(θ)...), proposal = prop, nrep))
        flush(stdout)
    end
end

function grid(seed, ntips, nrep)
    g, θ, targets, p = tree_of(ntips, Xoshiro(seed))
    L = mtbd_exact(g, PhyloPOMP.BDEI; θ, founder = :I)
    describe(g, θ, targets, p, L)
    a0, b0 = theory_ab(θ)
    mults = (0.25, 0.5, 1.0, 2.0, 4.0)
    for prop in (:guided, :hard)
        println("proposal $prop: sd(ll) [median(ll) − L]; rows a × multiplier, columns b × multiplier ", mults)
        for ma in mults
            cells = [(s = runs(g, θ, L; guide = guide_ab(ma * a0, mb * b0), proposal = prop, nrep);
                      @sprintf("%5.2f [%+6.2f]%s", s.sd, s.med, s.nfin < nrep ? "*" : " ")) for mb in mults]
            @printf("  a×%-4s %s\n", ma, join(cells, "  ")); flush(stdout)
        end
    end
end

## `trees`: (ntips, label, rng) for each tree; each tree is run with the filters of `filters`
function scale(out, trees, filters; nrep = 16)
    io = open(out, "w")
    println(io, "ntips\tlabel\tR0\tIP\tINC\tp\tL\tfilter\tNp\tmean_w\tse_w\tsd_ll\tmed_short\tnfinite\tnrep\tsec_per_run")
    for (ntips, label, pick) in trees
        g, θ, targets, p = pick()
        L = mtbd_exact(g, PhyloPOMP.BDEI; θ, founder = :I)
        describe(g, θ, targets, p, L); flush(stdout)
        a, b = theory_ab(θ)
        for (name, guide, prop, Np) in filters(a, b)
            s = runs(g, θ, L; guide, proposal = prop, Np, nrep)
            show_row("$name Np $Np", s); flush(stdout)
            println(io, join((nsample(g), label, targets..., p, L, name, Np, s.m, s.se, s.sd, s.med, s.nfin, nrep, s.t), '\t'))
            flush(io)
        end
    end
    close(io)
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = ARGS[1]
    num(i, d) = length(ARGS) >= i ? parse(Int, ARGS[i]) : d
    if mode == "pilot"
        pilot(num(2, 11), num(3, 20), num(4, 24))
    elseif mode == "grid"
        grid(num(2, 11), num(3, 20), num(4, 24))
    elseif mode == "scale"
        sizes = length(ARGS) >= 3 ? parse.(Int, split(ARGS[3], ',')) : [20, 50, 80, 150]
        trees = [(n, s, () -> tree_of(n, Xoshiro(1000n + s))) for n in sizes for s in 1:3]
        filters(a, b) = (("naive", nothing, :guided, 1000), ("default guided", guide_ab(DEFAULT...), :guided, 1000),
                         ("theory guided", guide_ab(a, b), :guided, 1000), ("theory hard", guide_ab(a, b), :hard, 1000),
                         ("theory guided", guide_ab(a, b), :guided, 4000))
        scale(ARGS[2], trees, filters)
    elseif mode == "moderate"
        ## R0 in [1.5, 3.5]: the first tree from successive seeds that has it
        pickmod(n, k) = () -> begin
            s = 0
            while true
                s += 1
                tr = tree_of(n, Xoshiro(9_000_000 + 1000n + 100k + s))
                1.5 <= tr[3][1] <= 3.5 && return tr
            end
        end
        trees = [(n, k, pickmod(n, k)) for n in (50, 80) for k in 1:2]
        filters(a, b) = (("naive", nothing, :guided, 1000), ("theory guided", guide_ab(a, b), :guided, 1000),
                         ("theory hard", guide_ab(a, b), :hard, 1000), ("theory guided", guide_ab(a, b), :guided, 4000))
        scale(ARGS[2], trees, filters)
    else
        error("mode must be pilot, grid, scale or moderate")
    end
end
