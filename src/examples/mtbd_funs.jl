## Definitions shared by the two-type multitype birth-death-sampling
## (MTBD) filters, NaiveMTBD (mtbd_naive.jl) and GuidedMTBD
## (mtbd_guided.jl).  Each module `include`s this file after declaring
## its demes, so `Camel` and `Human` refer to that module's own `Demes`.
##
## The model is the linear two-type birth-death-sampling process of
## Kuehnert et al. (2016) and Vaughan & Stadler (2025): a type-i host
## infects a new type-j host at rate lambda_ij, is removed unsampled at
## rate mu_i, is sampled (and removed: r = 1) at rate psi_i, and
## changes type at rate m_ij.  Type 1 is Camel and type 2 is Human, so
## trees may label tips `deme=camel`/`deme=human` or `deme=1`/`deme=2`.
##
## There are no susceptibles: every rate is per infected host.

"""
    epi_params(; R11, R12, R21, R22, delta1, delta2, s1, s2)

Convert the epidemiological parameterization of Vaughan & Stadler
(2025) to the event rates used by `filter_pomp`.  `Rij` is the expected
number of type-`j` hosts infected by one type-`i` host over its
infectious period, `deltai` the rate of becoming uninfectious, and `si`
the sampling proportion, so that

    lambdaij = Rij*deltai,   mui = deltai*(1-si),   psii = si*deltai.

Migration rates are zero.
"""
epi_params(; R11, R12, R21, R22, delta1, delta2, s1, s2) = (
    lambda11 = R11*delta1, lambda12 = R12*delta1,
    lambda21 = R21*delta2, lambda22 = R22*delta2,
    mu1 = delta1*(1-s1), mu2 = delta2*(1-s2),
    psi1 = s1*delta1, psi2 = s2*delta2,
    m12 = 0.0, m21 = 0.0,
)

"""
    check(gen)

Assert that the genealogy is one these filters can handle: one child at
each root, two at each branch point, and sample tips that are leaves
with a known deme.  Sampled ancestors would need removal probability
`r < 1`, which is not implemented.
"""
check(gen::Genealogy) = begin
    for node ∈ eachindex(gen)
        n = gen[node]
        if n.type == Root
            @assert length(n.children)==1 "wrong number of children ($(length(n.children)) ≠ 1) at root $(n.name)"
        elseif n.type == Sample
            @assert isempty(n.children) "sample $(n.name) has children: sampled ancestors need r < 1, which is not implemented"
            @assert n.deme ∈ (Camel, Human) "sample $(n.name) has unknown deme"
        elseif n.type == Node
            @assert length(n.children)==2 "wrong number of children ($(length(n.children)) ≠ 2) at node $(n.name)"
        end
    end
    nothing
end

## Initial numbers of infected hosts of each type.  The MTBD likelihood
## of Vaughan & Stadler starts from a single host, so the usual choice is
## I1_0 = 1, I2_0 = 0 (one camel).
mtbd_rinit(; I1_0, I2_0, _...) = (
    I1 = round(Int64, I1_0),
    I2 = round(Int64, I2_0),
)

## Default parameters: MLE for the 274-tip MERS tree (stem 1.0, sampling proportions 0.05).
const mers_mle = epi_params(
    R11 = 1.0385, R12 = 0.1393, R21 = 0.0, R22 = 0.9226,
    delta1 = 7.4170, delta2 = 60.9657, s1 = 0.05, s2 = 0.05,
)
