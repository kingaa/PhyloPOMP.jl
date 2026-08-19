# SEIR forward-simulator cross-validation: Julia vs R phylopomp

N = 1000 draws per side. Julia seed: `MersenneTwister(20260813)` (single
continuing stream). R seed: `set.seed(20260813)` (see
`scripts/seir_crossvalidate.R`). Parameters matched exactly on both sides
(β=4, σ=1, γ=1, ψ=0.30, χ=0, ω=1, N=100, S0=0.99, E0=0, I0=0.01, R0=0,
t0=0, tmax=20); see that script's header comment for why E0=0 (single
founding lineage) was chosen deliberately.

## Extinction fraction (reported separately from KS -- see script header)

| side  | fraction | 95% CI |
|-------|----------|--------|
| Julia | 0.243 | (0.216, 0.27) |
| R     | 0.237 | (0.211, 0.263) |

Two-proportion z-statistic (difference): z = 0.314
(not significant at α=0.05)

## nsample distribution (primary statistic)

|       | n   | mean | Q25/median/Q75 |
|-------|-----|------|-----------------|
| Julia | 1000 | 89.8 | (1, 115, 128) |
| R     | 1000 | 89.5 | (1, 113, 129) |

Two-sample KS: D = 0.036, p = 0.5361
(no evidence of a distributional difference at α=0.05)

## First-sample-time distribution (secondary statistic, conditional on nsample≥1)

|       | n (nonextinct) | mean | Q25/median/Q75 |
|-------|-----|------|-----------------|
| Julia | 757 | 1.794 | (0.7874397499327034, 1.5255642771171054, 2.394776420443183) |
| R     | 763 | 1.84 | (0.770735469576833, 1.6436067756519763, 2.5674056330046153) |

Two-sample KS: D = 0.0516, p = 0.2641
(no evidence of a distributional difference at α=0.05)
