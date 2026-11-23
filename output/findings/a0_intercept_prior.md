# A0 under three intercept priors (D25)

Produced by `analyses/probe_latent_groups/a0_intercept_prior.R` on 2026-11-22, real basin data (28 assemblages, counts sha256 `41bd753fe1c75396b09df992947a760de16739a9dcae0c7f5827f7adf351c7a7`). Exploratory design validation; the reported A0 fit is decided from this table.

Pooled log-odds against the reference class, m_c: -0.54, -1.84, -3.94, -3.58, -4.02, -5.06, -4.74, -5.42, -6.06.

| quantity | P2: Normal(0, 2), inherited | P10: Normal(0, 10) | PC: Normal(pooled log-odds, 2) |
|---|---|---|---|
| sigma | 4.132 [3.148, 5.200] | 3.016 [2.041, 4.497] | 2.257 [1.730, 3.196] |
| tau | 0.476 [0.110, 0.736] | 0.458 [0.086, 0.736] | 0.410 [0.061, 0.710] |
| rho (km) | 202.594 [109.997, 287.749] | 99.093 [37.154, 241.148] | 49.746 [24.898, 120.715] |
| sigma^2/rho | 0.087 [0.057, 0.130] | 0.095 [0.058, 0.148] | 0.103 [0.063, 0.160] |
| spatial share sigma^2/(sigma^2+tau^2) | 0.987 [0.961, 0.999] | 0.979 [0.926, 0.999] | 0.971 [0.897, 0.999] |
| field variance after removing its basin mean, sigma^2 c(rho) | 2.860 [1.890, 4.142] | 2.699 [1.785, 3.894] | 2.449 [1.663, 3.451] |
| contrast-based share sigma^2 c(rho) / (sigma^2 c(rho) + tau^2) | 0.927 [0.796, 0.997] | 0.928 [0.787, 0.998] | 0.937 [0.789, 0.999] |
| field contrast SD at 4.3 km | 0.862 [0.697, 1.050] | 0.893 [0.700, 1.101] | 0.922 [0.725, 1.133] |
| field contrast SD at 36.6 km | 2.414 [1.959, 2.910] | 2.373 [1.922, 2.848] | 2.287 [1.879, 2.713] |
| field contrast SD at 79.1 km | 3.360 [2.742, 4.031] | 3.124 [2.530, 3.785] | 2.804 [2.312, 3.413] |

## How the class baselines are carried

| class | baseline alpha + mean f (P2) | field mean, P2 | field mean, P10 | field mean, PC | corr(alpha, field mean) P2 / P10 / PC |
|---|---|---|---|---|---|
| 2 | -0.45 (sd 0.03) | -0.38 | -0.23 | -0.01 | -1.00 / -1.00 / -1.00 |
| 3 | -1.70 (sd 0.06) | -1.46 | -0.71 | -0.26 | -1.00 / -1.00 / -1.00 |
| 4 | -4.81 (sd 0.26) | -3.50 | 0.34 | -0.08 | -0.99 / -1.00 / -0.98 |
| 5 | -3.81 (sd 0.13) | -2.93 | -0.18 | -0.08 | -1.00 / -1.00 / -0.99 |
| 6 | -5.38 (sd 0.39) | -4.19 | -0.56 | -0.67 | -0.98 / -0.99 / -0.96 |
| 7 | -6.45 (sd 0.38) | -4.61 | 0.70 | -0.13 | -0.98 / -0.99 / -0.96 |
| 8 | -5.63 (sd 0.32) | -4.36 | -0.48 | -0.43 | -0.98 / -0.99 / -0.97 |
| 9 | -5.93 (sd 0.30) | -4.48 | -0.23 | -0.16 | -0.99 / -0.99 / -0.97 |
| 10 | -7.27 (sd 0.56) | -5.60 | -0.46 | -0.47 | -0.95 / -0.98 / -0.92 |

## Diagnostics

- P2: R-hat max 1.0027 (tau[1]) | bulk ESS min 1260 (tau[1]) | tail ESS min 2241 (tau[1]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.711 | stuck 0 (rung 2, 359 s, treedepth at max 0.0%)
- P10: R-hat max 1.0025 (lp__) | bulk ESS min 898 (tau[1]) | tail ESS min 1778 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.673 | stuck 0 (rung 1, 110 s, treedepth at max 0.0%)
- PC: R-hat max 1.0040 (lp__) | bulk ESS min 878 (lp__) | tail ESS min 1894 (lp__) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.696 | stuck 0 (rung 1, 60 s, treedepth at max 0.0%)

Visual set for P10 (Arial): `a0_alphaP10_rank.png`, `a0_alphaP10_trace.png`, `a0_alphaP10_pairs.png`, `a0_alphaP10_energy.png`.
