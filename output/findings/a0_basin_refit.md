# A0 refit on the real basin data, with corrected diagnostics and the rule-20 program

Produced by `analyses/probe_latent_groups/refit_a0_basin.R` on 2026-11-22 (pilot step 4 of the latent-group probe; exploratory design validation, not evidence about the record).

## Data contract and provenance

Prior in the contract: `prior_rho_type = 2` (1 inverse-gamma, 2 proportional to 1/rho, 3 uniform), rho on [0.5, 300] km, half-Normal(0, 2.0) on sigma and tau, Normal(0, 10) on the class intercepts (`prior_alpha_mean`, `prior_alpha_sd`; D25 made the intercept prior data on 2026-10-18, default Normal(0, 10); the 2026-10-16 run of this script used the hard-coded Normal(0, 2) and is kept as `a0_basin_refit_alphaN0-2.md`), kernel 1. This settles the lineage question of inventory row 14: the contract Carl's basin fit read carries the 1/rho prior; only the two repair scripts switched to the inverse-gamma.
Software: R version 4.6.2 (2026-10-29); rstan 2.32.7; StanHeaders 2.39.1.1; posterior 1.7.0; loo 2.10.1; bayesplot 1.16.0. Stan source sha256 `35fbf97122c9ab244be792a09756b4997094defe3fa1ea181890c703dae5b902`; fit layer sha256 `79d141d356a09fd4cef38c97343b0485791af1ffcf18c65d489543f84a6ab1a1`.

## 1. The refit

Gate: PASS. Rung 1 (iter 3000, warmup 1500, adapt_delta 0.95, max_treedepth 12); 116 s wall time (4 chains in parallel); read from cache.
- a0_basin: R-hat max 1.0025 (lp__) | bulk ESS min 898 (tau[1]) | tail ESS min 1778 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.673 | stuck 0
Treedepth, measured against the configured maximum for the first time: 0 of 6000 post-warmup iterations (0.0%) at depth 12; largest depth observed 10. Gated variables: 293 (sigma, tau, rho, 9 alpha, lp__, 252 field coordinates, 28 pointwise log-likelihoods); exempt 0; stuck 0.

| quantity | median | 5% | 95% |
|---|---|---|---|
| sigma | 3.016 | 2.041 | 4.497 |
| tau[1] | 0.458 | 0.086 | 0.736 |
| rho | 99.093 | 37.154 | 241.148 |
| spatial_share | 0.979 | 0.926 | 0.999 |
| half_distance | 68.686 | 25.753 | 167.151 |
| practical_range | 297.280 | 111.463 | 723.444 |

Visual diagnostics (Arial, 300 dpi, 7 in): `a0_basin_rank.png`, `a0_basin_trace.png`, `a0_basin_pairs.png`, `a0_basin_energy.png`; plus `a0_basin_pairs_alpha.png` (pairs on the nine intercepts).
Posterior-predictive counts for 200 draws: `output/probe_latent_groups/a0_posterior_predictive_yrep.csv` (summarized by the Python script).

## 2a. Prior predictive (rule 20a)

Sampled with the model's own `likelihood = 0` switch. - a0_prior_predictive: R-hat max 1.0048 (f[4,1]) | bulk ESS min 1315 (lp__) | tail ESS min 1525 (tau[1]) | divergences 0/4000 (0.00%) | treedepth at max 10: 0/4000 (0.0%) | E-BFMI min 0.965 | stuck 0
Counts for 200 draws: `output/probe_latent_groups/a0_prior_predictive_yrep.csv`. The resemblance summaries (mean within-assemblage Gini-Simpson, mean richness, pooled Gini-Simpson, two-cluster plug-in F_ST) are computed beside the observed values by `analyses/probe_latent_groups/a0_predictive_summaries.py`, which owns those definitions (spec 3.3); its output is `output/findings/a0_predictive_summaries.md`.

## 2b. Prior sensitivity (rule 20b)

Direction of each default and what the alternative leans toward: half-Normal(0, 2) on sigma and tau admits large variance of both kinds; half-Normal(0, 0.5) leans toward little of either. The 1/rho prior puts half its mass below 12.2 km (the geometric mean of 0.5 and 300); the uniform prior on [0.5, 300] puts 90 percent of its mass above 30 km and leans toward long ranges.

| quantity | default median [5%, 95%] | half-Normal(0, 0.5) on sigma, tau | uniform rho |
|---|---|---|---|
| sigma | 3.016 [2.041, 4.497] | 1.813 [1.506, 2.206] | 3.444 [2.266, 4.820] |
| tau[1] | 0.458 [0.086, 0.736] | 0.402 [0.058, 0.691] | 0.485 [0.125, 0.742] |
| rho | 99.093 [37.154, 241.148] | 32.499 [18.737, 59.500] | 138.665 [51.004, 271.021] |
| spatial_share | 0.979 [0.926, 0.999] | 0.955 [0.857, 0.999] | 0.981 [0.937, 0.999] |
| half_distance | 68.686 [25.753, 167.151] | 22.527 [12.987, 41.242] | 96.116 [35.353, 187.858] |

- a0_sens_scale05: R-hat max 1.0063 (tau[1]) | bulk ESS min 961 (tau[1]) | tail ESS min 1813 (lp__) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.676 | stuck 0
- a0_sens_rho_uniform: R-hat max 1.0068 (tau[1]) | bulk ESS min 930 (tau[1]) | tail ESS min 1817 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.654 | stuck 0

### What is identified (D20)

The data see between-site contrasts, with variance 2 sigma^2 (1 - exp(-d/rho)). On 28 sites within 79 km the contrast scale sigma^2/rho and the contrast SD at distances inside the basin are informed; sigma, rho and the intercepts individually are not (the intercepts are anti-correlated with the field's per-class mean, correlations -0.96 to -1.00 on the default fit). Distances below are the nearest-neighbor median, the median pair and the maximum pair. Length scales in km and sigma as an amplitude are not reportable as estimates; these are.

| quantity | default prior | half-Normal(0, 0.5) on sigma, tau | uniform rho |
|---|---|---|---|
| sigma^2/rho (logit^2 per km) | 0.09 [0.06, 0.15] | 0.10 [0.06, 0.16] | 0.09 [0.05, 0.14] |
| field contrast SD at 4.3 km (logits) | 0.89 [0.70, 1.10] | 0.91 [0.70, 1.11] | 0.87 [0.68, 1.07] |
| contrast SD with nugget at 4.3 km | 1.12 [0.96, 1.33] | 1.09 [0.93, 1.29] | 1.12 [0.94, 1.32] |
| field contrast SD at 36.6 km (logits) | 2.37 [1.92, 2.85] | 2.09 [1.73, 2.44] | 2.36 [1.89, 2.86] |
| contrast SD with nugget at 36.6 km | 2.46 [2.11, 2.88] | 2.18 [1.91, 2.48] | 2.46 [2.09, 2.90] |
| field contrast SD at 79.1 km (logits) | 3.12 [2.53, 3.79] | 2.42 [2.04, 2.84] | 3.20 [2.55, 3.88] |
| contrast SD with nugget at 79.1 km | 3.20 [2.65, 3.82] | 2.50 [2.17, 2.90] | 3.27 [2.70, 3.92] |

## 2c. Recovery against the prior (rule 20c)

Data simulated from the model by `analyses/49_write_stan_data.simulate_from_model` (numpy, no Stan) at the real geometry and sherd totals; spatial share 0.2 and 0.9 crossed with rho 5 and 60 km; total latent variance 1 logit^2. Covered = truth inside the 5 to 95 percent interval.

| record | quantity | truth | median | 5% | 95% | covered |
|---|---|---|---|---|---|---|
| share0.2_rho5 | spatial_share | 0.200 | 0.369 | 0.004 | 0.991 | yes |
| share0.2_rho5 | rho | 5.000 | 1.978 | 0.609 | 150.281 | yes |
| share0.2_rho5 | sigma | 0.447 | 0.624 | 0.063 | 1.036 | yes |
| share0.2_rho5 | tau[1] | 0.894 | 0.813 | 0.095 | 1.055 | yes |
| share0.2_rho60 | spatial_share | 0.200 | 0.163 | 0.001 | 0.956 | yes |
| share0.2_rho60 | rho | 60.000 | 3.006 | 0.588 | 209.589 | yes |
| share0.2_rho60 | sigma | 0.447 | 0.443 | 0.040 | 1.064 | yes |
| share0.2_rho60 | tau[1] | 0.894 | 0.992 | 0.227 | 1.136 | yes |
| share0.9_rho5 | spatial_share | 0.900 | 0.919 | 0.486 | 0.999 | yes |
| share0.9_rho5 | rho | 5.000 | 4.271 | 2.383 | 9.794 | yes |
| share0.9_rho5 | sigma | 0.949 | 0.975 | 0.721 | 1.111 | yes |
| share0.9_rho5 | tau[1] | 0.316 | 0.295 | 0.029 | 0.777 | yes |
| share0.9_rho60 | spatial_share | 0.900 | 0.918 | 0.797 | 0.972 | yes |
| share0.9_rho60 | rho | 60.000 | 90.220 | 26.490 | 265.624 | yes |
| share0.9_rho60 | sigma | 0.949 | 1.072 | 0.686 | 1.799 | yes |
| share0.9_rho60 | tau[1] | 0.316 | 0.330 | 0.235 | 0.416 | yes |

- share0.2_rho5: R-hat max 1.0044 (alpha[2]) | bulk ESS min 955 (alpha[2]) | tail ESS min 1236 (rho) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.706 | stuck 0
- share0.2_rho60: R-hat max 1.0023 (alpha[9]) | bulk ESS min 2054 (rho) | tail ESS min 2411 (alpha[6]) | divergences 0/24000 (0.00%) | treedepth at max 14: 0/24000 (0.0%) | E-BFMI min 0.674 | stuck 0
- share0.9_rho5: R-hat max 1.0023 (rho) | bulk ESS min 2416 (rho) | tail ESS min 1444 (rho) | divergences 0/32000 (0.00%) | treedepth at max 15: 0/32000 (0.0%) | E-BFMI min 0.672 | stuck 0
- share0.9_rho60: R-hat max 1.0061 (tau[1]) | bulk ESS min 1157 (rho) | tail ESS min 2405 (lp__) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.727 | stuck 0

## Timing and gates

| fit | rung | seconds | gate | reasons | cached |
|---|---|---|---|---|---|
| a0_basin | 1 | 116 | pass |  | yes |
| a0_prior_predictive | 1 | 3 | pass |  | yes |
| a0_sens_scale05 | 1 | 56 | pass |  | yes |
| a0_sens_rho_uniform | 1 | 127 | pass |  | yes |
| a0_recovery_share0.2_rho5 | 2 | 186 | pass |  | no |
| a0_recovery_share0.2_rho60 | 3 | 733 | pass |  | no |
| a0_recovery_share0.9_rho5 | 4 | 486 | pass |  | no |
| a0_recovery_share0.9_rho60 | 1 | 59 | pass |  | no |

Fit records (JSON, with every diagnostic, hash and version) are beside the cached fits in `output/gp_cache/`.
