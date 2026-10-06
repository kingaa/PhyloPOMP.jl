module MTBDTest

import ..Main: h1, h2

@info h1("MTBD model: naïve and guided filters against the exact likelihood")

using Test
using Random: seed!
using PhyloPOMP
using PhyloPOMP.NaiveMTBD, PhyloPOMP.GuidedMTBD
import PartiallyObservedMarkovProcesses as POMP

heavy = occursin(r"y|yes|t|true", get(ENV,"RUN_HEAVY_TESTS","no"))

## Exact log likelihoods ("ranked" convention, which is phylopomp's and
## this package's) from the ODE recursion of Kuehnert et al. (2016), as
## implemented in mtbd2_exact.py (mers_private/mtbd-analysis); that
## implementation reproduces phylopomp's lbdp_exact in the single-type
## limit.  All cases: r = 1, one camel founder (I1_0 = 1, I2_0 = 0).
const T1 = "(((t1[&&PhyloPOMP deme=1]:0.5,t2[&&PhyloPOMP deme=1]:0.8):0.4,t3[&&PhyloPOMP deme=1]:1.1):0.3,(t4[&&PhyloPOMP deme=1]:0.9,t5[&&PhyloPOMP deme=1]:0.2):0.6):0.5;"
const T2 = "((t1[&&PhyloPOMP deme=1]:0.3,t2[&&PhyloPOMP deme=2]:1.0)B:0.6,t3[&&PhyloPOMP deme=2]:1.5)A:0.4;"
const T4 = "((t1[&&PhyloPOMP deme=1]:0.3,t2[&&PhyloPOMP deme=2]:1.0)B:0.6,t3[&&PhyloPOMP deme=2]:1.5)A:0;"
const single = (lambda11=1.5, lambda12=0.0, lambda21=0.0, lambda22=0.0,
                mu1=0.6, mu2=0.0, psi1=0.5, psi2=0.0, m12=0.0, m21=0.0)
const two = (lambda11=1.4, lambda12=0.5, lambda21=0.3, lambda22=1.9,
             mu1=0.5, mu2=0.7, psi1=0.4, psi2=0.6, m12=0.0, m21=0.0)
const cases = [
    ## name,                       tree,    parameters,                       exact
    ("single type, 5 tips",        T1,      single,                            -2.791776834),
    ("two types, 3 tips",          T2,      two,                               -4.929549691),
    ("two types with migration",   T2,      merge(two, (m12=0.2, m21=0.35)),   -4.822551362),
    ("two types, zero stem",       T4,      two,                               -4.199766476),
    ## One tip each.  In "camel tip" the only camel is tracked at the
    ## start, so a camel-to-human birth that leaves the lineage in place
    ## is possible although no untracked camel exists: a proposal that
    ## ignores this underestimates the likelihood (see event_rates! in
    ## mtbd_naive.jl).
    ("camel tip",                  "t1[&&PhyloPOMP deme=1]:1.0;",  two,     -1.243148680),
    ("human tip",                  "t1[&&PhyloPOMP deme=2]:1.0;",  two,     -1.919820527),
    ("camel/human cherry",         "(t1[&&PhyloPOMP deme=1]:0.5,t2[&&PhyloPOMP deme=2]:0.5):0.5;", two, -1.547447482),
    ("human/human cherry",         "(t1[&&PhyloPOMP deme=2]:0.5,t2[&&PhyloPOMP deme=2]:0.5):0.5;", two, -1.524962407),
]

## Particle-filter likelihood estimates are unbiased on the natural scale,
## so average replicates with logmeanexp and compare to the exact value
## in units of the resulting standard error.
zscore(p, exact; reps=200, Np=1000) = begin
    ll = [POMP.logLik(pfilter(p, Np=Np)) for _ ∈ 1:reps]
    @test all(isfinite, ll)
    est, se = POMP.logmeanexp(ll, se=true)
    (est - exact)/se
end

## Any guide gives an unbiased filter; it only changes variance. Conductance 3
## keeps the guide informative and the replicate spread small.
guide_proc(M) = fsmarkov(M.Camel=>0.5, M.Human=>0.5, (M.Camel,M.Human)=>3.0)

@testset verbose=true "MTBD model" begin

    seed!(5048671)

    @testset "parameters and genealogy checks" begin
        e = NaiveMTBD.epi_params(R11=1.0, R12=0.5, R21=0.2, R22=0.9,
                       delta1=10.0, delta2=20.0, s1=0.1, s2=0.05)
        @test e.lambda12 ≈ 5.0 && e.lambda21 ≈ 4.0
        @test e.mu1 + e.psi1 ≈ 10.0 && e.psi2 ≈ 1.0
        g = parse_newick(T2, t0=0, demes=NaiveMTBD.Demes)
        @test NaiveMTBD.check(g) === nothing
        ## a sampled ancestor would need r < 1, which is not implemented
        sa = parse_newick("((t1[&&PhyloPOMP deme=1]:0.5)s[&&PhyloPOMP deme=1]:0.5);",
                          t0=0, demes=NaiveMTBD.Demes)
        @test_throws AssertionError NaiveMTBD.check(sa)
        ## no founder: the tree is impossible
        p = NaiveMTBD.filter_pomp(g; I1_0=0, I2_0=0, two...)
        @test logLik(pfilter(p, Np=100)) == -Inf
    end

    for (name, nwk, th, exact) ∈ cases
        @testset "$name" begin
            gn = parse_newick(nwk, t0=0, demes=NaiveMTBD.Demes)
            gg = parse_newick(nwk, t0=0, demes=GuidedMTBD.Demes)
            zn = zscore(NaiveMTBD.filter_pomp(gn; th...), exact)
            zg = zscore(GuidedMTBD.filter_pomp(gg, guide_proc(GuidedMTBD); th...), exact)
            @info h2("$name: z = $(round(zn,digits=2)) (naïve), $(round(zg,digits=2)) (guided)")
            @test abs(zn) < 4
            @test abs(zg) < 4
        end
    end

    if heavy
        @info h2("274-tip MERS tree at the maximum-likelihood estimate (known failure)")
        ## Both filters return -Inf on the 274-tip MERS tree (known failure).
        ## Exact log likelihood at mers_mle with a 1-year stem: +381.8742.
        nwk = replace(PhyloPOMP.GuidedMERS.mers_newick, r":0;$" => ":1;")
        g = parse_newick(nwk, t0=0, demes=GuidedMTBD.Demes)
        p = GuidedMTBD.filter_pomp(g, fsmarkov(GuidedMTBD.Camel=>0.99, GuidedMTBD.Human=>0.01,
                                               (GuidedMTBD.Camel,GuidedMTBD.Human)=>0.01))
        @test_broken isfinite(logLik(pfilter(p, Np=2000)))
    end

end

end
