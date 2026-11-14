# A3 under the inherited and the corrected intercept prior

Produced by `analyses/probe_latent_groups/a3_prior_comparison.py` from the findings tables of `validate_a3.R` and `a3_pinned_range.R`. Posterior medians of omega (SD in logits of a region-level shift), pooled by record type over gate-passing fits only. Descriptive; exploratory design validation, not evidence about the record. Prior on omega half-Normal(0, 2): median 1.35.

## free range

| records | Normal(0, 2), inherited | Normal(0, 10), corrected |
|---|---|---|
| gate | 39 of 39 pass (rung 1 in 27) | 39 of 39 pass (rung 1 in 30) |
| own-generator, omega 0 | mean 0.77, range 0.47 to 1.16 (n 6) | mean 0.80, range 0.48 to 1.21 (n 6) |
| own-generator, omega 0.25 | mean 0.75, range 0.44 to 1.14 (n 6) | mean 0.75, range 0.43 to 1.08 (n 6) |
| own-generator, omega 0.5 | mean 1.05, range 0.64 to 2.05 (n 6) | mean 1.07, range 0.61 to 1.86 (n 6) |
| own-generator, omega 1 | mean 1.51, range 0.76 to 2.59 (n 6) | mean 1.41, range 0.70 to 2.28 (n 6) |
| smooth slope, rho 30 km, no barrier | mean 1.73, range 0.97 to 2.49 (n 3) | mean 1.44, range 0.78 to 2.22 (n 3) |
| smooth slope, rho 60 km, no barrier | mean 0.74, range 0.57 to 0.85 (n 3) | mean 0.82, range 0.58 to 1.12 (n 3) |
| simulator S0, no barrier | mean 1.85, range 1.40 to 2.65 (n 3) | mean 1.01, range 0.80 to 1.43 (n 3) |
| simulator S1 and S2, barrier | mean 1.16, range 0.27 to 2.85 (n 6) | mean 1.06, range 0.47 to 1.41 (n 6) |

## range pinned at 6 km

| records | Normal(0, 2), inherited | Normal(0, 10), corrected |
|---|---|---|
| gate | 39 of 39 pass (rung 1 in 34) | 39 of 39 pass (rung 1 in 31) |
| own-generator, omega 0 | mean 3.30, range 2.83 to 3.49 (n 6) | mean 1.30, range 0.69 to 1.58 (n 6) |
| own-generator, omega 0.25 | mean 2.98, range 2.58 to 3.35 (n 6) | mean 1.17, range 0.81 to 1.57 (n 6) |
| own-generator, omega 0.5 | mean 3.43, range 3.12 to 3.72 (n 6) | mean 1.75, range 1.42 to 2.04 (n 6) |
| own-generator, omega 1 | mean 3.66, range 3.01 to 4.20 (n 6) | mean 1.79, range 1.05 to 2.50 (n 6) |
| smooth slope, rho 30 km, no barrier | mean 3.40, range 2.92 to 3.97 (n 3) | mean 2.02, range 1.25 to 2.81 (n 3) |
| smooth slope, rho 60 km, no barrier | mean 3.34, range 2.87 to 3.58 (n 3) | mean 1.88, range 1.56 to 2.17 (n 3) |
| simulator S0, no barrier | mean 3.66, range 3.35 to 3.86 (n 3) | mean 1.27, range 0.97 to 1.79 (n 3) |
| simulator S1 and S2, barrier | mean 2.32, range 0.60 to 3.91 (n 6) | mean 1.17, range 0.62 to 1.39 (n 6) |

## range pinned at 12 km

| records | Normal(0, 2), inherited | Normal(0, 10), corrected |
|---|---|---|
| gate | 39 of 39 pass (rung 1 in 30) | 37 of 39 pass (rung 1 in 22); fail: own_G1_om0_r3, own_L_om0.25_r3 |
| own-generator, omega 0 | mean 3.15, range 2.74 to 3.40 (n 6) | mean 1.24, range 1.03 to 1.51 (n 5) |
| own-generator, omega 0.25 | mean 2.81, range 2.34 to 3.25 (n 6) | mean 1.04, range 0.57 to 1.44 (n 5) |
| own-generator, omega 0.5 | mean 3.31, range 2.91 to 3.58 (n 6) | mean 1.58, range 1.22 to 1.96 (n 6) |
| own-generator, omega 1 | mean 3.58, range 2.93 to 4.19 (n 6) | mean 1.64, range 0.81 to 2.38 (n 6) |
| smooth slope, rho 30 km, no barrier | mean 3.13, range 2.59 to 3.75 (n 3) | mean 1.62, range 0.86 to 2.52 (n 3) |
| smooth slope, rho 60 km, no barrier | mean 3.00, range 2.61 to 3.27 (n 3) | mean 1.36, range 1.01 to 1.82 (n 3) |
| simulator S0, no barrier | mean 3.65, range 3.39 to 3.87 (n 3) | mean 1.23, range 0.94 to 1.73 (n 3) |
| simulator S1 and S2, barrier | mean 2.27, range 0.57 to 3.95 (n 6) | mean 1.15, range 0.60 to 1.40 (n 6) |

## Contrast scale of the generating field

| field | sigma | rho (km) | between-site contrast SD (logits) |
|---|---|---|---|
| own-generator records (A0 point under Normal(0, 2), `a0_pinned_nuisances.json`) | 4.13 | 204.1 | 1.66 |
| step-against-slope records | 4.13 | 30.0 | 3.32 |
| step-against-slope records | 4.13 | 60.0 | 2.72 |
| A0 posterior median, Normal(0, 2), inherited | 4.13 | 204.1 | 1.70 (sqrt of the median of sigma^2 c(rho)) |
| A0 posterior median, Normal(0, 10) | 3.01 | 98.1 | 1.64 (sqrt of the median of sigma^2 c(rho)) |
| A0 posterior median, Normal(pooled log-odds, 2) | 2.27 | 50.6 | 1.56 (sqrt of the median of sigma^2 c(rho)) |

