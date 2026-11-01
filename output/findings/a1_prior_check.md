# A1 under the corrected intercept prior (D25 follow-up)

Produced by `analyses/probe_latent_groups/a1_prior_check.R` on 2026-10-19. Three A1 fits on the D25 records (prior_alpha_sd = 10; seed 20261017, the seed of validate_a1.R), against the Normal(0, 2) posteriors stored in `a1_validation_summary.json`. Median [5%, 95%]. The gate and its reasons are validate_a1.R's. Exploratory design validation, not evidence about the record; a descriptive comparison, no threshold.

| record | true sigma_g | sigma_g, Normal(0, 2) | sigma_g, Normal(0, 10) | Delta, Normal(0, 2) | Delta, Normal(0, 10) | Delta_w, Normal(0, 2) | Delta_w, Normal(0, 10) | coherence, Normal(0, 2) | coherence, Normal(0, 10) | gate (new) | rung | seconds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| zero | 0.00 | 0.517 [0.048, 1.386] | 0.575 [0.060, 1.466] | 0.975 [0.093, 2.190] | 1.078 [0.117, 2.349] | 0.910 [0.057, 1.971] | 1.026 [0.086, 2.160] | 0.003 [-0.062, 0.095] | 0.006 [-0.056, 0.097] | fail: rhat,ess_bulk | 2 | 773 |
| sg0.5_r1 | 0.50 | 0.272 [0.029, 0.924] | 0.257 [0.023, 0.829] | 0.522 [0.055, 1.374] | 0.493 [0.045, 1.304] | 0.479 [0.025, 1.146] | 0.458 [0.033, 1.190] | 0.000 [-0.069, 0.079] | 0.000 [-0.065, 0.082] | fail: rhat | 2 | 994 |
| sg0.5_r5 | 0.50 | 0.268 [0.028, 1.098] | 0.269 [0.027, 0.873] | 0.506 [0.054, 1.451] | 0.515 [0.053, 1.318] | 0.448 [0.015, 1.135] | 0.476 [0.027, 1.157] | 0.000 [-0.066, 0.084] | 0.000 [-0.068, 0.080] | fail: ess_tail | 2 | 793 |

## Diagnostics (rule 16)

- zero: R-hat max 1.0297 (comembership[22,19]) | bulk ESS min 231 (comembership[21,19]) | tail ESS min 423 (comembership[9,7]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.697 | stuck 0
- sg0.5_r1: R-hat max 1.0128 (comembership[22,21]) | bulk ESS min 656 (comembership[20,3]) | tail ESS min 922 (comembership[26,3]) | divergences 0/12000 (0.00%) | treedepth at max 13: 1/12000 (0.0%) | E-BFMI min 0.677 | stuck 0
- sg0.5_r5: R-hat max 1.0086 (comembership[22,20]) | bulk ESS min 575 (comembership[27,25]) | tail ESS min 372 (comembership[20,10]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.701 | stuck 0

Visual set for the null fit: `figures/probe_latent_groups/a1prior_null_{rank,trace,pairs,energy}.png`.

