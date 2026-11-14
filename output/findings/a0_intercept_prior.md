# A0 under three intercept priors (D25)

Produced by `analyses/probe_latent_groups/a0_intercept_prior.R` on 2026-10-18, real basin data (28 assemblages, counts sha256 `41bd753fe1c75396b09df992947a760de16739a9dcae0c7f5827f7adf351c7a7`). Exploratory design validation; the reported A0 fit is decided from this table.

Pooled log-odds against the reference class, m_c: -0.54, -1.84, -3.94, -3.58, -4.02, -5.06, -4.74, -5.42, -6.06.

| quantity | P2: Normal(0, 2), inherited | P10: Normal(0, 10) | PC: Normal(pooled log-odds, 2) |
|---|---|---|---|
| sigma | 4.134 [3.165, 5.286] | 3.015 [2.041, 4.444] | 2.269 [1.723, 3.236] |
| tau | 0.466 [0.089, 0.734] | 0.453 [0.075, 0.744] | 0.404 [0.057, 0.720] |
| rho (km) | 204.142 [109.260, 287.255] | 98.139 [37.574, 236.734] | 50.596 [24.522, 126.029] |
| sigma^2/rho | 0.088 [0.056, 0.131] | 0.095 [0.058, 0.146] | 0.103 [0.061, 0.159] |
| spatial share sigma^2/(sigma^2+tau^2) | 0.988 [0.961, 1.000] | 0.979 [0.925, 0.999] | 0.971 [0.896, 0.999] |
| field variance after removing its basin mean, sigma^2 c(rho) | 2.879 [1.866, 4.172] | 2.693 [1.767, 3.872] | 2.445 [1.632, 3.482] |
| contrast-based share sigma^2 c(rho) / (sigma^2 c(rho) + tau^2) | 0.929 [0.794, 0.998] | 0.929 [0.783, 0.998] | 0.938 [0.784, 0.999] |
| field contrast SD at 4.3 km | 0.867 [0.692, 1.053] | 0.892 [0.701, 1.098] | 0.922 [0.714, 1.132] |
| field contrast SD at 36.6 km | 2.423 [1.946, 2.920] | 2.368 [1.913, 2.849] | 2.286 [1.858, 2.721] |
| field contrast SD at 79.1 km | 3.371 [2.738, 4.040] | 3.119 [2.505, 3.792] | 2.812 [2.290, 3.433] |

## How the class baselines are carried

| class | baseline alpha + mean f (P2) | field mean, P2 | field mean, P10 | field mean, PC | corr(alpha, field mean) P2 / P10 / PC |
|---|---|---|---|---|---|
| 2 | -0.45 (sd 0.03) | -0.30 | -0.27 | -0.06 | -1.00 / -1.00 / -1.00 |
| 3 | -1.70 (sd 0.05) | -1.43 | -0.83 | -0.25 | -1.00 / -1.00 / -1.00 |
| 4 | -4.80 (sd 0.26) | -3.47 | 0.43 | -0.11 | -0.99 / -0.99 / -0.98 |
| 5 | -3.81 (sd 0.13) | -2.94 | -0.28 | -0.12 | -1.00 / -1.00 / -0.99 |
| 6 | -5.39 (sd 0.39) | -4.21 | -0.68 | -0.66 | -0.98 / -0.99 / -0.96 |
| 7 | -6.44 (sd 0.38) | -4.68 | 0.63 | -0.10 | -0.98 / -0.99 / -0.96 |
| 8 | -5.64 (sd 0.32) | -4.34 | -0.45 | -0.43 | -0.98 / -0.99 / -0.97 |
| 9 | -5.93 (sd 0.29) | -4.50 | -0.12 | -0.13 | -0.99 / -0.99 / -0.97 |
| 10 | -7.27 (sd 0.55) | -5.58 | -0.45 | -0.45 | -0.96 / -0.98 / -0.92 |

## Diagnostics

- P2: R-hat max 1.0034 (lp__) | bulk ESS min 608 (tau[1]) | tail ESS min 1249 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.676 | stuck 0 (rung 1, 89 s, treedepth at max 0.0%)
- P10: R-hat max 1.0059 (lp__) | bulk ESS min 874 (tau[1]) | tail ESS min 1666 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.672 | stuck 0 (rung 1, 146 s, treedepth at max 0.0%)
- PC: R-hat max 1.0050 (lp__) | bulk ESS min 1053 (tau[1]) | tail ESS min 1928 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.720 | stuck 0 (rung 1, 53 s, treedepth at max 0.0%)

Visual set for P10 (Liberation Sans): `a0_alphaP10_rank.png`, `a0_alphaP10_trace.png`, `a0_alphaP10_pairs.png`, `a0_alphaP10_energy.png`.
