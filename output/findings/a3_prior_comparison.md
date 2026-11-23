# A3 under the two class-intercept priors

Produced by `analyses/probe_latent_groups/a3_prior_comparison.py` from the findings tables of `validate_a3.R` and `a3_pinned_range.R`. Posterior medians of omega (SD in logits of a region-level shift), pooled by record type over gate-passing fits only. Descriptive; exploratory design validation, not evidence about the record. Prior on omega half-Normal(0, 2): median 1.35.

## free range

| records | Normal(0, 2) | Normal(0, 10) |
|---|---|---|
| gate | 39 of 39 pass (rung 1 in 25) | 39 of 39 pass (rung 1 in 25) |
| own-generator, omega 0 | mean 0.77, range 0.53 to 1.13 (n 6) | mean 0.60, range 0.45 to 0.78 (n 6) |
| own-generator, omega 0.25 | mean 0.73, range 0.45 to 0.84 (n 6) | mean 0.67, range 0.51 to 0.84 (n 6) |
| own-generator, omega 0.5 | mean 1.16, range 0.65 to 2.48 (n 6) | mean 0.98, range 0.43 to 1.66 (n 6) |
| own-generator, omega 1 | mean 1.53, range 0.80 to 2.58 (n 6) | mean 1.41, range 0.82 to 2.31 (n 6) |
| smooth slope, rho 30 km, no barrier | mean 1.82, range 1.00 to 2.38 (n 3) | mean 1.12, range 0.65 to 1.49 (n 3) |
| smooth slope, rho 60 km, no barrier | mean 0.82, range 0.76 to 0.88 (n 3) | mean 0.66, range 0.57 to 0.81 (n 3) |
| simulator S0, no barrier | mean 1.86, range 1.38 to 2.69 (n 3) | mean 1.04, range 0.81 to 1.43 (n 3) |
| simulator S1 and S2, barrier | mean 1.13, range 0.28 to 2.60 (n 6) | mean 1.05, range 0.47 to 1.38 (n 6) |

## range pinned at 6 km

| records | Normal(0, 2) | Normal(0, 10) |
|---|---|---|
| gate | 39 of 39 pass (rung 1 in 34) | 39 of 39 pass (rung 1 in 32) |
| own-generator, omega 0 | mean 3.22, range 2.87 to 3.43 (n 6) | mean 1.07, range 0.78 to 1.37 (n 6) |
| own-generator, omega 0.25 | mean 2.97, range 2.49 to 3.45 (n 6) | mean 1.05, range 0.72 to 1.31 (n 6) |
| own-generator, omega 0.5 | mean 3.50, range 3.14 to 3.81 (n 6) | mean 1.59, range 1.09 to 1.92 (n 6) |
| own-generator, omega 1 | mean 3.61, range 3.00 to 4.21 (n 6) | mean 1.75, range 1.07 to 2.52 (n 6) |
| smooth slope, rho 30 km, no barrier | mean 3.49, range 3.08 to 3.93 (n 3) | mean 1.60, range 1.08 to 1.93 (n 3) |
| smooth slope, rho 60 km, no barrier | mean 3.44, range 3.22 to 3.59 (n 3) | mean 1.48, range 1.35 to 1.75 (n 3) |
| simulator S0, no barrier | mean 3.64, range 3.35 to 3.83 (n 3) | mean 1.27, range 0.95 to 1.79 (n 3) |
| simulator S1 and S2, barrier | mean 2.36, range 0.61 to 3.96 (n 6) | mean 1.17, range 0.63 to 1.40 (n 6) |

## range pinned at 12 km

| records | Normal(0, 2) | Normal(0, 10) |
|---|---|---|
| gate | 39 of 39 pass (rung 1 in 33) | 31 of 39 pass (rung 1 in 20); fail: own_G1_om0_r3, own_L_om0.25_r1, own_L_om0.25_r2, own_L_om0_r1, own_L_om0_r2, own_L_om0_r3, own_L_om1_r3, slope_rho60_r2 |
| own-generator, omega 0 | mean 3.10, range 2.81 to 3.33 (n 6) | mean 1.03, range 0.82 to 1.24 (n 2) |
| own-generator, omega 0.25 | mean 2.86, range 2.39 to 3.40 (n 6) | mean 0.96, range 0.62 to 1.10 (n 4) |
| own-generator, omega 0.5 | mean 3.39, range 3.02 to 3.77 (n 6) | mean 1.39, range 0.72 to 1.84 (n 6) |
| own-generator, omega 1 | mean 3.53, range 2.91 to 4.21 (n 6) | mean 1.74, range 1.09 to 2.35 (n 5) |
| smooth slope, rho 30 km, no barrier | mean 3.16, range 2.58 to 3.73 (n 3) | mean 1.27, range 0.72 to 1.68 (n 3) |
| smooth slope, rho 60 km, no barrier | mean 3.13, range 3.02 to 3.31 (n 3) | mean 1.28, range 1.11 to 1.45 (n 2) |
| simulator S0, no barrier | mean 3.63, range 3.41 to 3.82 (n 3) | mean 1.24, range 0.94 to 1.74 (n 3) |
| simulator S1 and S2, barrier | mean 2.23, range 0.58 to 3.93 (n 6) | mean 1.14, range 0.59 to 1.39 (n 6) |

## Contrast scale of the generating field

| field | sigma | rho (km) | between-site contrast SD (logits) |
|---|---|---|---|
| own-generator records of the Normal(0, 10) column (the A0 medians in `a0_pinned_nuisances.json`, written from `a0_basin_46cff13716e8eebb.rds`) | 3.02 | 99.1 | 1.65 |
| own-generator records of the Normal(0, 2) column (the A0 medians in `a0_pinned_nuisances_alphaN0-2.json`, written from `a0_alpha_P2_3d93d6258d9aa823.rds`) | 4.13 | 202.6 | 1.67 |
| step-against-slope records | 3.02 | 30.0 | 2.42 |
| step-against-slope records | 3.02 | 60.0 | 1.98 |
| A0 posterior median, Normal(0, 2), inherited | 4.13 | 202.6 | 1.69 (sqrt of the median of sigma^2 c(rho)) |
| A0 posterior median, Normal(0, 10) | 3.02 | 99.1 | 1.64 (sqrt of the median of sigma^2 c(rho)) |
| A0 posterior median, Normal(pooled log-odds, 2) | 2.26 | 49.7 | 1.56 (sqrt of the median of sigma^2 c(rho)) |

