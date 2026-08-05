module HardMERSTest

import ..Main: h1, h2

@info h1("MERS model with hard-guided proposals")

using Test
using BenchmarkTools
using Random: seed!
using PhyloPOMP
using PhyloPOMP.HardMERS
using PhyloPOMP.HardMERS.Demes: Camel, Human
import PartiallyObservedMarkovProcesses as POMP

heavy = occursin(r"y|yes|t|true", get(ENV,"RUN_HEAVY_TESTS","yes"))

## See test/mers_soft.jl for why a small reduced fixture (rather than the
## full empirical `mers_tree`) is used for the simulate/pfilter checks.
const small_tree = "(([&&PhyloPOMP deme=camel]1:1.0,[&&PhyloPOMP deme=camel]2:1.0):0.5,[&&PhyloPOMP deme=human]3:1.5);"

@testset verbose=true "MERS model with hard-guided proposals" begin

    seed!(2121916527)

    m0 = fsmarkov(Camel=>0.3, Human=>0.7, (Camel,Human)=>1)
    p = HardMERS.filter_pomp(HardMERS.mers_tree, m0, I_c0=0, I_h0=0)
    @test p isa POMP.PompObject
    @test logLik(pfilter(p,Np=50))==-Inf

    @info h2("constructor validation")
    ## proposal_floor=0 is legal (disables the floor -- "no floor, let it
    ## fail"); only out-of-[0,1] values are rejected.
    p0 = HardMERS.filter_pomp(HardMERS.mers_tree, m0, proposal_floor=0.0)
    @test p0 isa POMP.PompObject
    @test_throws ArgumentError HardMERS.filter_pomp(HardMERS.mers_tree, m0, proposal_floor=-0.1)
    @test_throws ArgumentError HardMERS.filter_pomp(HardMERS.mers_tree, m0, proposal_floor=1.5)

    g = parse_newick(small_tree, demes=HardMERS.Demes)
    m = fsmarkov(Camel=>0.3, Human=>0.7, (Camel,Human)=>1)
    common = (
        β_hc=1.0, β_ch=1.0, χ_c=1.0, χ_h=1.0,
        S_c0=1.0, S_h0=1.0, I_c0=0.5, I_h0=0.5, N_c=100.0, N_h=100.0,
    )
    p = HardMERS.filter_pomp(g, m; common...)
    @test p isa POMP.PompObject

    @info h2("simulate test")
    sm = simulate(p, nsim = 3)
    @time sm = simulate(p, nsim = 3)
    @test sm isa Matrix{<:POMP.PompObject}

    @info h2("pfilter test")
    ## Np is large enough that the per-particle probability of matching the
    ## required cross-deme transmission makes total collapse to -Inf
    ## negligible regardless of seed (see test/mers_soft.jl).
    pf = pfilter(p, Np = 5000)
    @time pf = pfilter(p, Np = 5000)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(logLik(pf))

    @info h2("deterministic kernel tests")

    ## ell=0: no branches to boost, so the identity outcome is certain.
    kappa0_e, _ = HardMERS.cross_proposal(9, 0, 0.05)
    @test kappa0_e == 1
    @test isempty(HardMERS.hard_branch_law(Float64[], 0, 9, 0.05))

    ## I=ell boundary and vanishing/non-finite hazards: every outcome keeps
    ## strictly positive auxiliary intensity.
    I, ell, eps = 4, 4, 0.05
    kappa0, _ = HardMERS.cross_proposal(I, ell, eps)
    @test kappa0 > 0
    kappaC = HardMERS.hard_branch_law([0.0, NaN, Inf, 2.0], ell, I, eps)
    @test all(isfinite.(kappaC)) && all(kappaC .> 0)

    @info h2("hard auxiliary-intensity / compensating-decay identity")
    ## beta_j = alpha_population * kappa_epsilon(j); since the raw relative
    ## hazards are not renormalized to ell/I (contrast SoftMERS), the
    ## kappa's need not sum to 1, and the shortfall/excess relative to the
    ## population rate must appear entirely in the compensating decay.
    I2, ell2, eps2 = 10, 3, 0.1
    rh = [0.2, 5.0, 1.0]
    kappa0_2, _ = HardMERS.cross_proposal(I2, ell2, eps2)
    kappaC_2 = HardMERS.hard_branch_law(rh, ell2, I2, eps2)
    alpha_population = 4.0
    beta = alpha_population .* vcat(kappa0_2, kappaC_2)
    decay_compensation = alpha_population - sum(beta)
    @test alpha_population*(kappa0_2 + sum(kappaC_2)) + decay_compensation ≈ alpha_population
    ## realized correction for branch b: log(alpha_population*Phi_b/beta_b)
    ## must equal log(Phi_b/kappa_epsilon(b)) -- i.e. the alpha_population
    ## factor cancels exactly.
    Phi = [0.7, 0.05, 0.9]
    for b ∈ eachindex(Phi)
        @test isapprox(
            log(alpha_population*Phi[b]/beta[b+1]),
            log(Phi[b]/kappaC_2[b]),
        )
    end

    ## Aggregate-plus-conditional equivalence: drawing the aggregate
    ## "boosted swap" outcome (probability sum(kappaC)) and then the
    ## specific branch conditionally on kappaC (probability
    ## kappaC[b]/sum(kappaC)) combines to the full per-outcome
    ## kappa_epsilon(b), exactly as used in `HardMERS.regular_part!`.
    aggregate = sum(kappaC_2)
    conditional = kappaC_2 ./ aggregate
    @test all((aggregate .* conditional) .≈ kappaC_2)

    @info h2("proposal_floor=0: no floor, let it fail")
    ## epsilon=0 disables the floor: cross_proposal reduces to the plain
    ## naive split, and a branch with genuinely zero relative hazard gets
    ## exactly zero auxiliary intensity -- unlike eps>0 above, where every
    ## entry (including the [0.0, NaN, Inf] cases) stayed strictly positive.
    kappa0_3, _ = HardMERS.cross_proposal(I2, ell2, 0.0)
    @test kappa0_3 ≈ (I2-ell2)/I2
    kappaC_3 = HardMERS.hard_branch_law([0.0, 5.0, 1.0], ell2, I2, 0.0)
    @test kappaC_3[1] == 0
    @test all(kappaC_3[2:end] .> 0)

    if heavy
        @info h2("pfilter benchmark")
        @btime pfilter($p, Np = 1000)

        @time ll = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
        llest,llse = logmeanexp(ll,se=true)
        @info "logLik = $(round(llest,digits=2)) ± $(round(llse,sigdigits=3))"
    end

end

end
